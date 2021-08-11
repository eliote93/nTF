MODULE TYPEDEF

!Define

TYPE ASYINFO_TYPE                !ASSEMBLY INFORMATION
  logical :: lempty, lfuel, lUse, lgeom
  logical :: lCentX, lCentY, lCentXY
  INTEGER :: nx, ny, nxy
  INTEGER :: nFuelPin = 0
  INTEGER :: GapType = 0
  INTEGER,pointer :: pin(:), pin2DIdx(:,:)
  INTEGER :: EdgeAsyIdx(0:3) = 0                 !1: X-dir edge, 2:Y-dir edge, 3:both
  INTEGER :: nThChGrp                            !n of T-H Chennel Group
  INTEGER, POINTER :: ThChGrp(:)                 !
END type

TYPE ASYGAP_TYPE
  LOGICAL :: lEmpty = .TRUE.
  INTEGER :: nx, ny
  INTEGER, POINTER :: GapPin(:), GapPin2D(:, :)
END TYPE

TYPE Asy_Type
  logical :: lCentX, lCentY, lCentXY, ldum
  INTEGER :: ixa, iya, color                       !--- CNJ Edit : Red-Black Domain Decomposition
  INTEGER :: AsyType
  INTEGER :: PartialAsyFlag                        ! Partial Assembly Flag
  INTEGER,pointer :: NeighIdx(:)                   ! Neighborhood Assembly Index
  INTEGER, POINTER :: GlobalPinIdx(:)              !
  INTEGER, POINTER :: PinIdx_Ext(:, :)             !
  REAL :: CX, CY                                   ! Assemlby Center Location
  REAL :: wt = 1.0
  !Control Rod
  LOGICAL :: lCrAsy = .FALSE.                       !
  INTEGER :: nCrPin = 0                            !
  INTEGER, POINTER :: CrPinIdx(:)                  ! Global Pin Idx
  INTEGER, POINTER :: CrPinBankId(:)               !
END TYPE

type coreinfo_type
  INTEGER :: nxya,nxa,nya   !Assembly Number
  INTEGER :: nxy,nx,ny      !Actual Cell numbe
  INTEGER :: nxc0, nyc0, nxyc0
  INTEGER :: nxc, nyc, nxyc
  INTEGER :: nxy0,nx0,ny0   !
  INTEGER :: nCoreFXR, nCoreFSR  ! # of Fxr and Fsr in the Core
  INTEGER :: nFuelPin = 0   !Fuel Pin Number
  INTEGER :: nAsyCell, nAsyGT         !Number of GT per Asy
  INTEGER :: nz, nzfm, nsubplane
  LOGICAL :: lGap = .FALSE.
  LOGICAL :: lEdge = .TRUE.
  LOGICAL :: lRot = .FALSE.
  LOGICAL :: lCbd = .FALSE.            !Checker Board
  LOGICAL :: RadSym(6)      !Radial Symmetry Flag
  INTEGER :: RadBC(6)       !Radial Boundary Condition
  INTEGER :: AxBC(2)        !Axial

  REAL :: PowerCore = 0._8      !Core Power
  REAL :: Hm_Mass0 = 0._8
  REAL :: FuelVol, FuelVolFm, TotVol       !Fuel Volume and TotalVolume
  REAL :: DelT = 0.0                         !Time Step Size
  Logical, POINTER :: lFuelPlane(:), lCladPlane(:), lAICPlane(:)
  INTEGER :: nBasePlane
  INTEGER, pointer :: CoreIdx(:,:)
  INTEGER, pointer :: CoreMap(:)
  INTEGER, POINTER :: THChMap(:)
  INTEGER, POINTER :: SubPlaneMap(:), SubPlaneRange(:, :)
  INTEGER, POINTER :: nSubPlanePtr(:)
  INTEGER :: nPinType, nCellType, nAsyType, nCellType0
  REAL, POINTER :: hz(:), hzfm(:), HzFmInv(:), HzInv(:)
  REAL,POINTER :: AsyCentX(:), AsyCentY(:)
  REAL, POINTER :: PinVol(:, :), PinVolFm(:, :)
  REAL, POINTER :: AsyVol(:, :)
  !--- CNJ Edit : CASMO Linear Source
  REAL, POINTER :: CoreFsrCentroid(:, :, :)
  REAL, POINTER :: CoreFsrMxx(:, :), CoreFsrMyy(:, :), CoreFsrMxy(:, :)
  !--- CNJ Edit : Group-wise Albedo (BYS Request)
  REAL, POINTER :: groupAlbedo(:, :)
  TYPE(Pin_Type), POINTER :: Pin(:)
  TYPE(ResVarPin_Type), POINTER :: ResVarPin(:,:)
  TYPE(Asy_Type), POINTER :: Asy(:)
  TYPE(AsyInfo_type), POINTER :: AsyInfo(:)
  TYPE(pininfo_type), POINTER :: PinInfo(:)
  TYPE(cell_type), POINTER :: CellInfo(:), BaseCellInfo(:) ! --- 180724 JSR
  TYPE(AsyGap_Type), POINTER :: AsyGap(:)
  LOGICAL :: lDcpl = .FALSE.
  !INTEGER :
  INTEGER :: iRefPln
  INTEGER :: nMappingPln, MappingPln(100)
  !
  !LOGICAL :: lCrIn = .FALSE.
  !INTEGER :: nCrConf, nCrBank
  !TYPE(CrConf_Type), POINTER :: CrConf(:)
  !TYPE(CrBank_Type), POINTER :: CrBank(:)
  LOGICAL :: lCrInfo = .FALSE.
  TYPE(CrInfo_Type), POINTER :: CrInfo
  TYPE(MiscStruct_Type), POINTER :: MiscStruct
END type

type pininfo_type
  logical :: lempty, luse, lfuel
  logical :: lCentX, lCentY, lCentXY
  logical :: lgap
  INTEGER :: ncell                      !number of axial cell
  INTEGER :: nFsrMax                    !Maximum number of FSR in a pincell
  INTEGER :: nFxrMax                    !Maxium number of FXR in a Pincell
  INTEGER :: PartialAsyFlag             !Partial Assembly Ray Index 0:None, 1:x, 2: y, 3:x-y
  INTEGER,pointer :: cell(:)            !Cell Type Info
  INTEGER :: EdgePinIdx(0:3)
  INTEGER :: GapType

END type

type basicgeom
  INTEGER :: nbd
  logical :: lcircle,lrect
  logical :: lCCent = .FALSE.
  INTEGER :: ncircle,nline
  INTEGER :: CCentType
  INTEGER :: nx,ny
  REAL,pointer :: bdline(:,:)  !Boundary Line
  REAL,pointer :: bdlineRange(:,:)
  REAL,pointer :: line(:,:)
  REAL,pointer :: circle(:,:)
  REAL :: cx,cy, lx,ly, x(2), y(2)
  REAL,POINTER :: delx(:), dely(:)
  INTEGER :: inpnx, inpny
  INTEGER :: inpdivx(10), inpdivy(10)
  INTEGER :: inpn       ! --- 180723 JSR
  INTEGER :: inpdiv(20) ! --- 180723 JSR
  REAL :: inpdelx(10), inpdely(10)
END type

type pin_type
  !logical :: lcad                            !Cad Input type
  LOGICAL :: lfuel = .FALSE.
  LOGICAL :: lGT = .FALSE.
  LOGICAL :: lRadRef = .FALSE.
  LOGICAL :: lBaffle = .FALSE.
  LOGICAL :: lBarrel = .FALSE.
  LOGICAL :: lGd = .FALSE.
  INTEGER :: ix, iy                           !Global Cartesian Coordinate
  INTEGER :: PinType                          !Pintype
  INTEGER :: AsyType                          !
  INTEGER :: ncell                            !number of axial cell
  INTEGER :: nFsrMax                          !Maximum number of FSR in a pincell
  INTEGER :: nFxrMax                          !Maxium number of FXR in a Pincell
  INTEGER :: iasy, ipin                       !Assembly Index and Local Pin Index
  INTEGER :: isuperPin                        !SuperPin Index
  INTEGER, POINTER :: cell(:)
  LOGICAL, POINTER :: lMox(:)
  LOGICAL, POINTER :: lAIC(:)
  INTEGER :: nBd                              !Number of Boundary
  INTEGER :: FsrIdxSt                         !Flat source Region(FSR) Index Starting Point
  INTEGER :: FxrIdxSt                         !Flat XS region Index(FXR) Starting Point
  !INTEGER :: nFxr
  REAL, POINTER :: BdLength(:)
  REAL, POINTER :: Center2SurfaceL(:)         !Pin Center to surface Length
  INTEGER, POINTER :: NeighIdx(:)             !Neighborhood Index
  INTEGER, POINTER :: NeighSurfIdx(:)         !Neigh Surf
  INTEGER, POINTER :: hCelGeo(:)              !Hexagonal Cell Geometry Base - KSC 180904
  !CR Info
  LOGICAL :: lCrPin = .FALSE.
  LOGICAL :: lInsert = .FALSE.
  INTEGER :: CrBankId = 0
  REAL, POINTER :: FracCrIn(:)
END type


type ResVarPin_type
  LOGICAL :: lres=.false.,lresA=.false.,lresC=.false.,lfuel = .FALSE.
  LOGICAL :: lcrres = .FALSE.
  REAL(4), POINTER :: avgxseq(:,:,:),avgxseq_1g(:),avgxseq_mg(:,:)
  REAL(4), POINTER :: rifa(:,:),rifs(:,:),riff(:,:)
  REAL(4), POINTER :: FnAdj(:)
  INTEGER :: niso,igt, niso_Res
  INTEGER, POINTER :: idiso(:), idx_Res(:), idiso_Res(:)
  REAL,POINTER :: pnum(:), pnum_Res(:)!, pnumrat(:,:) ! # of nuclides, # of nuclides removing duplcated in resonance, ND ratio
  ! Equivalent Dancoff Cell method
  ! Dancoff : Dancoff factor
  ! sigpF   :
  ! siglpM  :
  ! lbar    : mean chord length
  ! EquiRad : Dancoff factor equivalent radius of 1D cell
  ! invsurf4:
  REAL :: Dancoff=0._8, sigpF=0._8, siglpM=0._8, lbar=0._8, EquiRad=0._8, invsurf4, XSEsc = 0._8
  REAL,POINTER :: rad_cp(:),X(:,:,:)
  ! Exist but not used in code
  REAL, POINTER :: delr(:), Qsurfvol(:),vol(:)
  ! Variables for PSM
  INTEGER :: icelpsm
  REAL, POINTER :: Pji(:,:,:), eta(:) ! First-flight CP j<-i / Shadowing correction factor
  REAL :: fuelvol, nonfuelvol
  INTEGER :: nfueldiv, nfuelring
  REAL ::  siglpMcat(3), mavgm(3), alpham(3), invsubalpham(3), sigpm_inval(3) ! for grouped scattering nuclides
  LOGICAL :: lmsrc(3)
! REAL :: ALPHA1, ALPHA2, BETA1, BETA2
END type

type cell_type
  logical :: lempty ,luse, lrect, lres, lfuel, lGd   !Empty, Used, Rectangular cell, Resonance , Fuel
  LOGICAL :: lCR= .FALSE.
  LOGICAL :: lMox = .FALSE.
  LOGICAL :: lAIC = .FALSE.
  LOGICAL :: lsSPH = .FALSE.
  logical :: lCCell = .FALSE.                        !Corner Cell
  logical :: lCentX, lCentY, lCentXY
  logical :: lgap
  logical :: lcad = .FALSE.
  LOGICAL :: lhole = .FALSE.  ! EDIT JSU 20190626
  INTEGER :: ibFuel, ieFuel   ! Beginning and End of fuel cell...
  INTEGER :: nDivAzi                              !Azimuthal Division
  INTEGER :: nBd
  INTEGER :: nFSR
  INTEGER :: nFXR
  INTEGER :: nCellRay                     !
  INTEGER :: EdgeCellIdx(0:3) = 0
  INTEGER :: GapType
  INTEGER :: icel0
  INTEGER :: iReg0(0:300)
  INTEGER :: prcell
   REAL :: rr(300) ! --- 180914 JSR
  INTEGER :: iSSPHcell !--- 210414 JSU Added | Spectral SPH geom index
  !Volume & Geometry Index

  REAL,pointer :: vol(:)
  REAL, POINTER :: MaxSeg(:, :)              !Maximum Segment Length
  type(basicgeom):: Geom
  !XS Index Info
  INTEGER,POINTER :: iReg(:)                       !Composition or Mixture number of FSR
  INTEGER,POINTER :: FxrIdxSt(:)                   !First region number beloning to this xs regio
  INTEGER,POINTER :: nFsrInFxr(:)                  !Number of flat source region which is include in the FXR
  INTEGER,POINTER :: MapFxr2FsrIdx(:,:)            !Flat XS region Index Starting Point-> FSR Starting Index
  TYPE(CellRayInfo_type),POINTER :: CellRay(:)     !Cell Ray Information
  TYPE(THCell_Type), POINTER :: THCell
  TYPE(CadGeom_Type), POINTER :: CadGeom
  CHARACTER(246) :: CellInput
  !CR                                              !
  LOGICAL :: lCrCell = .FALSE.                     !
  TYPE(CrCell_Type), POINTER :: CrCell             ! Control Rod Cell
  INTEGER :: basecellstr                           !--- BYS edit : Base Cell Structure
  !Spectral SPH Factor related
  INTEGER :: icelSSPH
  REAL(4),POINTER :: SPHfactor(:,:)
  INTEGER,POINTER :: matidx(:),matfxridx(:)
  REAL,POINTER :: matrad(:),rad_cp(:),fxrvol(:),q_cp(:),matvol(:)
  INTEGER :: nreg_cp,nmat,nmodfxr,cldfxridx,nfueldiv,ngapdiv,ncladdiv,srdidx
  REAL :: fuelgapcldvol,invnmodfxr,FuelRefTEMP0,FuelRad0,U238ND0,MODXSV0
  ! PSM
  REAL    :: nonfuelvol, fuelvol
  INTEGER :: nnonfuel
  INTEGER :: icelPSM = 0 ! EDUT JSU 20200611
  LOGICAL :: lPSMcel = .FALSE.
END type

!Rad Decomposition
TYPE RadDcmp_TYPE
  LOGICAL :: lRadDcmp = .TRUE.
  INTEGER :: nDom                          ! Number of Domain
  INTEGER :: nxylocal(100)                 ! Number of local pin cell
  !INTEGER :: nxbeg(100), nybeg(100), nxend(100), nyend(100)
  INTEGER, POINTER :: nxbeg(:,:),nxend(:,:)
  INTEGER, POINTER :: nybeg(:), nyend(:)
  INTEGER, POINTER :: PinIdx(:, :)         ! Pin Index
  INTEGER, POINTER :: PinDom(:)
END TYPE

!Ray Related Type

TYPE CellRayInfo_type
  LOGICAL :: luse                            !Used or Non-Used Option
  INTEGER :: idxst, idxend
  INTEGER :: nSeg                            !Number of Segments
  INTEGER, POINTER :: LocalFsrIdx(:)         !Local Fsr(within pincell structure) Index
  INTEGER, POINTER :: LenSegInt(:)           !Segment Length Information(micron Unit)
  !INTEGER :: InOutSurf(2)                   !In-Out Surface Input
  REAL, POINTER :: LenSeg(:)                 !Segment Length

  REAL, POINTER :: pts(:,:)
END TYPE

TYPE Element_Type
  LOGICAL :: lboundary = .FALSE.
  INTEGER :: itype = 2
  INTEGER :: nnode = 3
  INTEGER :: nsurf = 3
  INTEGER :: GeomType, PhyType
  INTEGER :: NodeIdx(3)
  INTEGER :: NeighElement(3)=0
  INTEGER :: NeighSurf(3)=0
  INTEGER :: SurfNode(2,3)
  REAL :: x(4), y(4)
  REAL :: vol
  REAL :: line(3,3)
END TYPE

TYPE CadGeom_Type
  INTEGER :: nElement
  INTEGER :: nnode
  REAL, POINTER :: x(:), y(:)
  CHARACTER(50) :: MeshFormat, idum
  TYPE(Element_Type), POINTER :: Element(:)
END TYPE

!TYPE CellRayBase_Type
!  INTEGER :: itype
!  INTEGER :: AziAngleStIdx(100)              !Azimuthal Angle Index
!  TYPE(CellRayInfo_type), POINTER :: CellRay
!END TYPE

TYPE AsyRayInfo_type                         !Assembly Ray Information
  !LOGICAL :: lCentX, lCentY, lCentXY         !Partial Modular Ray Flag(1/2 or 1/4 in case of CENT option card)

  INTEGER :: AziAngIdx                       !Azimuthal Angle Index Info
  INTEGER :: ModRayIdx                       !Modular Ray Index
  INTEGER :: NextRayIdx(2)                   !NextAsyRayIdx
  INTEGER :: InOutSurf(2)                   !
  INTEGER :: nCellRay
  INTEGER, POINTER :: PinRaySurf(:,:)        !Pin Ray In & Out Surface Number
  INTEGER, POINTER :: PinIdx(:)              !Pin Index of Assembly Ray
  INTEGER, POINTER :: PinRayIdx(:)           !Pin Ray Index
  INTEGER :: PartialAsyRayFlag =0            !Partial Assembly Ray Index 0:None, 1:x, 2: y, 3:x-y
  INTEGER :: PartialAsyRayIdx(0:3)             !1: x, 2: y, 3: x-y
  REAL :: InOutPoint(2,2)
END TYPE

TYPE ModRayInfo_type
  INTEGER :: NextRayIdx(2)                   !In And Out Ray Index
  INTEGER :: InOutSurf(2)                    !
  INTEGER :: iAziAngIdx                      !
  INTEGER :: idum1                           !
  REAL :: InOutPoint(2,2)                    !
END TYPE

TYPE CoreRayInfo_type                        !Set of Modular Ray(Core Boundary to Boundary)
  INTEGER :: nRay                            !Number of Modular Ray
  INTEGER :: iAng
  INTEGER, POINTER :: AsyRayIdx(:)           !Modular Ray Index Numbers info
  INTEGER, POINTER :: AsyIdx(:)              !Assembly Idx
  INTEGER, POINTER :: ixa(:), iya(:)
  INTEGER :: InOutSurf(2)
END TYPE

TYPE RotRayInfo_Type
  INTEGER :: nRay                            !number of CORE RAY
  INTEGER :: nSeg                            !Number of total segment belong to
  INTEGER :: Ang1, Ang2
  INTEGER :: OmpAng1, OmpAng2
  INTEGER, POINTER :: RayIdx(:)              !Core Ray Index
  INTEGER, POINTER :: Dir(:)                 !Searching
END TYPE

!--- CNJ Edit : Domain Decomposition
TYPE DcmpAsyRayInfo_Type
  INTEGER :: iRotRay
  INTEGER :: iAsy        ! Global Index of Assembly
  INTEGER :: iRay        ! Numeric # for each Assembly
  INTEGER :: nMaxCellRay
  INTEGER :: nAsyRay     ! # of Reflections

  INTEGER, POINTER, DIMENSION(:) :: AsyRayList ! (iRef), Index of Modular Ray
  INTEGER, POINTER, DIMENSION(:) :: AziList    ! (iRef), Index of Azimuthal Angle
  INTEGER, POINTER, DIMENSION(:) :: DirList    ! (iRef). Index of Direction

  LOGICAL :: lRotRayBeg(2) = .FALSE. ! T : End Point is Starting Point of Rotational Ray
END TYPE

!--- CNJ Edit : Angular Multigrid Ray Tracing
TYPE MultigridInfo_Type
  INTEGER, POINTER :: AziList(:)
  INTEGER :: nAzi, nPolar
  REAL, POINTER :: wtang(:, :), wtsurf(:, :, :)
  REAL, POINTER :: mwt(:, :, :), mwt2(:, :, :), Comp(:, :, :)
  REAL, POINTER :: EXPA(:, :), EXPB(:, :)
END TYPE

TYPE FastCoreRayDat_Type
   INTEGER, POINTER :: PinIdx(:, :)
   INTEGER, POINTER :: CellIdx(:, :)
   INTEGER, POINTER :: CellRayIdxSt(:, :, :)
   INTEGER, POINTER :: SurfIdx(:, :, :)
   INTEGER, POINTER :: Ray1DIdx(:, :, :)
   INTEGER, POINTER :: nTotRaySeg(:)
   INTEGER, POINTER :: nTotCellRay(:)

   TYPE(FastRaySegDat_Type), POINTER :: RaySeg(:)
END TYPE

TYPE FastRaySegDat_Type
  INTEGER, POINTER :: FsrIdx(:)
  REAL, POINTER :: Lenseg(:)
  REAL, POINTER :: OptLen(:, :)
  REAL, POINTER :: ExpApp(:, :)
END TYPE


TYPE AziAngleInfo_Type
  INTEGER :: nModRay                        !Number of Modular Ray
  INTEGER :: nAsyRay                        !Number of Assembly Ray
  INTEGER :: nCoreRay                       !Number of Core Ray
  INTEGER :: ModRayIdxSt                    !Starting Index of Modular Ray
  INTEGER :: AsyRayIdxSt                    !Starting Index of Assembly Ray
  INTEGER :: CoreRayIdxSt                   !Starting Index of Core Ray
  REAL :: ang, sinv, cosv, tanv, del, weight
  TYPE(ModRayInfo_type), POINTER :: ModRay(:)
  TYPE(AsyRayInfo_type), POINTER :: AsyRay(:)
END TYPE

TYPE PolarAngle_Type
  REAL :: ang, sinv, cosv, weight
END TYPE

TYPE RayInfo_Type
  INTEGER :: nAziAngle = 8                                !Azimuthal Angle
  INTEGER :: nPolarAngle = 2                              !Polar Angle
  INTEGER :: nPolarAngleHemi
  INTEGER :: nModRay, nAsyRay, nCoreRay, nRotRay
  INTEGER :: nPhiAngSv
  INTEGER, POINTER :: PhiangInSvIdx(:, :), PhiangOutSvIdx(:, :)  !Angular FLux Save Array Adress(Index)
  !--- CNJ Edit : Domain Decomposition
  INTEGER, POINTER :: DcmpAsyAziList(:, :, :)
  INTEGER, POINTER :: DcmpAsyLinkInfo(:, :, :, :)
  INTEGER, POINTER :: DcmpAsyRayCount(:)
  !--- CNJ Edit : Node Majors, GPU Acceleration
  INTEGER, POINTER :: RotRayAziList(:, :)
  REAL :: del = 0.05                                             !Ray Spacing
  TYPE(AziAngleInfo_Type), POINTER :: AziAngle(:)
  TYPE(PolarAngle_Type), POINTER :: PolarAngle(:)
  TYPE(ModRayInfo_type), POINTER :: ModRay(:)
  TYPE(AsyRayInfo_type), POINTER :: AsyRay(:)
  TYPE(CoreRayInfo_Type), POINTER :: CoreRay(:)
  TYPE(RotRayInfo_Type), POINTER :: RotRay(:)
  !TYPE(RayInfo4CMFD_TYPE), POINTER :: RayInfo4CMFD
  TYPE(DcmpAsyRayInfo_Type), POINTER :: DcmpAsyRay(:, :)   !--- CNJ Edit : Domain Decomposition
  TYPE(CellRayInfo_type), POINTER :: CellRay1D
  TYPE(FastCoreRayDat_Type), POINTER :: FastCoreRayDat(:, :)
  TYPE(MultigridInfo_Type), POINTER :: MultigridInfo(:)
END TYPE


TYPE TrackingDat_Type
  LOGICAL :: lAlloc       = .FALSE.
  LOGICAL :: lAllocP1     = .FALSE.
  LOGICAL :: lAllocLinSrc = .FALSE.
  LOGICAL :: lAllocNM     = .FALSE.
  
  INTEGER, POINTER, DIMENSION(:,:)   :: FsrIdx, ExpAppIdx, AziMap
  INTEGER, POINTER, DIMENSION(:,:,:) :: ExpAppIdxnm
  
  ! Basic
  REAL, POINTER, DIMENSION(:,:)   :: wtang, hwt
  REAL, POINTER, DIMENSION(:,:,:) :: wtsurf, comp, mwt, mwt2
  
  ! GM TRC
  REAL, POINTER, DIMENSION(:)     :: phis, src, xst, PhiAngOut
  REAL, POINTER, DIMENSION(:,:)   :: OptLenList, EXPA, EXPB, ExpApp, PhiAngOutPolar, PhiAngIn
  REAL, POINTER, DIMENSION(:,:,:) :: ExpAppPolar, jout
  
  ! GM Pn
  REAL, POINTER, DIMENSION(:,:)   :: phim
  REAL, POINTER, DIMENSION(:,:,:) :: SrcAng1, SrcAng2
  
  ! NM Trc
  REAL, POINTER, DIMENSION(:,:)     :: phisnm, srcnm, xstnm
  REAL, POINTER, DIMENSION(:,:,:)   :: PhiAngOutnm, PhiAngInnm, OptLenListnm
  REAL, POINTER, DIMENSION(:,:,:,:) :: joutnm, ExpAppnm
  
  ! NM Pn
  REAL, POINTER, DIMENSION(:,:,:)   :: phimnm
  REAL, POINTER, DIMENSION(:,:,:,:) :: SrcAngnm1, SrcAngnm2
  
  ! AFSS
  REAL, POINTER, DIMENSION(:,:,:,:)   :: phia1g
  REAL, POINTER, DIMENSION(:,:,:,:,:) :: phiaNg
  
  ! CASMO LS
  REAL, POINTER, DIMENSION(:)       :: FsrMxx, FsrMyy, FsrMxy
  REAL, POINTER, DIMENSION(:,:)     :: FsrCentroid
  REAL, POINTER, DIMENSION(:,:,:)   :: dMxx, dMyy, dMxy, srcSlope, phimx, phimy, q0, x0, y0
  REAL, POINTER, DIMENSION(:,:,:,:) :: dCentroid, E1, E3, R1, R3, cmOptLen, cmOptLenInv, q1
  
  ! Dcmp.
  REAL, POINTER, DIMENSION(:,:,:,:)   :: DcmpPhiAngIn1g, DcmpPhiAngOut1g
  REAL, POINTER, DIMENSION(:,:,:,:,:) :: DcmpPhiAngInNg, DcmpPhiAngOutNg
END TYPE

TYPE RayInfo4CMFD_TYPE
  !Data Information for

  INTEGER :: nRotRay, nCoreRay, nPolAngle
  INTEGER, POINTER :: RotRayInOutCell(:, :)
  INTEGER, POINTER :: PhiangInSvIdx(:, :)
  REAL, POINTER :: PhiAngIn(:, :, :, :)
  TYPE(PolarAngle_Type), POINTER :: PolarAngle(:)
  !--- CNJ Edit : Domain Decomposition
  REAL, POINTER :: AsyPhiAngIn(:, :, :, :, :, :)
  INTEGER, POINTER :: DcmpAsyRayInCell(:, :, :), DcmpAsyRayInSurf(:, :, :)
  INTEGER, POINTER :: DcmpAsyRayCount(:)
END TYPE

TYPE CellRayBase_SearchTree

  REAL :: Mult = 100000
  INTEGER :: nline, idum1, ndiv
  INTEGER,POINTER :: nPtsInSruf(:)                     !Number of points which is belong to same surf
  INTEGER,POINTER :: InSurfIdxSt(:)       !In-Surface
  INTEGER,POINTER :: InOutSurf(:,:)
  INTEGER, POINTER :: pts(:,:,:)
END TYPE

TYPE CellRayBase_Type

  INTEGER :: CellType
  INTEGER :: nline                          !Number of Cell Ray
  INTEGER, POINTER :: AziAngIdx(:)          !Array Starting index of given angle index
  INTEGER, POINTER :: nlines(:)             !Number of Cell Rays which are belong to a certain angle
  TYPE(CellRayInfo_type), POINTER :: CellRay(:)
  TYPE(CellRayBase_SearchTree), POINTER :: SearchTree(:)
ENDTYPE


!CMFD Related Typeb
type scatmat

    INTEGER ib,ie
    REAL,pointer :: from(:), dummy
    REAL :: WithInGroupScat
    REAL :: self   !--- CNJ Edit : Preserving Self-scattering
END type

TYPE PinXs_Type

  TYPE(SCATMAT),POINTER,DIMENSION(:) :: XSS
  REAL,POINTER, DIMENSION(:,:) :: DTIL, DHAT, PDHAT, DTIL2      !Radial Dhat and Dtil
  REAL, POINTER, DIMENSION(:, :, :) :: partialDhat   !--- CNJ Edit : p-CMFD Acceleration
  REAL, POINTER, DIMENSION(:, :) :: atil, ahat       !--- CNJ Edit : Domain Decomposition
  REAL,POINTER, DIMENSION(:,:) :: AxDtil, AxDhat
  REAL,POINTER, DIMENSION(:) :: XSD, XSD2, XST, XSTR, XSR, XSNF, XSKF, CHI, XSA, XSDA !,PHI,FAC,XSRD  !--- BYS edit : XSA added
  REAL, POINTER, DIMENSION(:) :: PHI, FAC, XSRD
  REAL, POINTER, DIMENSION(:,:,:) :: Dcpl_DHAT
  REAL :: FuelTemp, ModTemp, PinTemp

  !Kinetic Parameter
  REAL :: BETAT, BETAP, OMEGA
  REAL :: xstr1g, xstr1g0                  !Absoprtion XS 1group
  !REAL :: BETAT,OMEGA,BETAP
  REAL,POINTER,DIMENSION(:) :: BETA, CHIP, VELO,RVDELT
END TYPE

TYPE MicXsFtn_Type

  TYPE(SCATMAT),POINTER :: XSS(:, :, :)
  REAL, POINTER :: sig_T(:, :, :), sig_Nf(:, :, :), sig_kf(:, :, :)
  REAL :: Tfuel, Tmod, Tclad, Tst
END TYPE
!Axial Solver TYPE

TYPE AXFLX_TYPE

  LOGICAL :: lUse = .FALSE.
  REAL,POINTER :: PHI(:,:,:) !FLUX, Modlar Space Flux
  REAL,POINTER :: PSI(:)  !Fission Source
  REAL,POINTER :: TLKG(:,:,:), LkgShape(:, :)  !, SubNodePhi(:, :), SubNodeLkg(:, :)
  REAL,POINTER :: JOUT(:,:,:), DHAT(:,:), DTIL(:, :), PDHAT(:, :)
  REAL, POINTER :: SP3DHAT(:, :, :), SP3DTIL(:, :, :), SP3PDHAT(:, :, :)
  REAL,POINTER :: S(:,:), KSQ(:,:), QT(:,:) !EIGENVECTOR,KAFFA SQUARE,QHAT TRANS
  REAL, POINTER :: partialDhat(:, :, :)   !--- CNJ Edit : p-CMFD Acceleration
  REAL, POINTER :: TranPsi(:), TranPsid(:), Prec(:, :), TranSrc(:, :)
END TYPE

TYPE AxGEOM_TYPE
  INTEGER :: IX, IY, NG, idum
  INTEGER :: myzb, myze, myzbf, myzef
  INTEGER :: NMESH, NCOMP
  INTEGER :: BC(2)
  REAL :: Area
  REAL, POINTER :: H(:), HINV(:)
  INTEGER, POINTER :: COMP(:), idum2
  LOGICAL :: lTransient = .FALSE.
END TYPE

TYPE BiLU_TYPE

  REAL, POINTER :: DelInv(:,:), Deliau(:, :), Al(:, :)
  REAL, POINTER :: DelInv2G(:, :,:), Deliau2G(:, :, :), Al2G(:, :, :)
  LOGICAL :: lAlloc = .FALSE.
  LOGICAL :: lReset = .TRUE.
  LOGICAL :: l2G = .FALSE.
END TYPE

!TYPE BiLU2G_TYPE
!
!  REAL, POINTER :: DelInv(:, :,:), Deliau(:, :, :), Al(:, :, :)
!  LOGICAL :: lAlloc = .FALSE.
!  LOGICAL :: lReset = .TRUE.
!END TYPE

TYPE CMFDLS_TYPE

  REAL, POINTER :: Diag(:, :)
  REAL, POINTER :: RadOffDiag(:, :, :), AxOffDiag(:, :, :)
  REAL, POINTER :: Diag2g(:, :, :)
  REAL, POINTER :: RadOffDiag2g(:, :, :, :), AxOffDiag2g(:, :, :, :)
  TYPE(BiLU_TYPE) :: BiLU
  !TYPE(BiLU2G_TYPE) :: BiLU2G
  INTEGER, POINTER :: NeighIdx(:, :)
  INTEGER, POINTER :: AxialPlaneMap(:)
  INTEGER :: myzbf, myzef, nxy, nbd
  INTEGER :: nproc, comm, myrank   !For MPI solver : number of Process, communicator, rank
  INTEGER :: nThread
  LOGICAL :: l2G = .FALSE.
END TYPE


TYPE DcplCmfdLs_Type

  TYPE(CMFDLS_Type), POINTER :: CMFDLS(:)
END TYPE

TYPE FMInfo_TYPE

  REAL, POINTER :: Phis(:, :, :) !nFsr, iz, ng
  REAL, POINTER :: Phi1a(:, :, :, :, :) !Polar, nFsr, Azimuthal, iz, ng
  REAL, POINTER :: Phi2a(:, :, :, :, :) !Polar, nFsr, Azimuthal, iz, ng
  REAL, POINTER :: PhiAngIn(:, :, :, :)
  REAL, POINTER :: AsyPhiAngIn(:, :, :, :, :, :)   !--- CNJ Edit : Domain Decomposition
  REAL, POINTER :: Psi(:, :), PsiD(:, :)
  REAL, POINTER :: PsiC(:, :), PsicD(:, :)
  REAL, POINTER :: Power(:, :)
  REAL, POINTER :: RadJout(:, :, :, :, :)
  REAL, POINTER :: AxSrc(:, :, :), AxPXS(:, :, :)
  REAL, POINTER :: phim(:, :, :, :)
  REAL, POINTER :: LinSrcSlope(:, :, :, :)
  REAL, POINTER :: PhiCrit(:), SpecConv(:)

  REAL, POINTER :: gPhis(:, :, :) !-- JSU EDIT 2017/08/10

  REAL,POINTER :: w(:)   !Under Relaxation Factor
  TYPE(FXRInfo_TYPE), POINTER :: Fxr(:, :)

  !Transient Variables----------------------------------------------------------------------------------
  REAL, POINTER :: TranPower(:, :)
  REAL, POINTER :: TranPhi(:, :, :)
  REAL, POINTER :: Prec(:, :, :)                     !Precursor Number density
  REAL, POINTER :: PrecSrc(:, :)
  REAL, POINTER :: PrecSrcK(:, :, :)
  REAL, POINTER :: TranPsi(:, :), TranPsid(:, :)
  REAL, POINTER :: ResSrc(:, :, :)                   !Connected with CmInfo%ResSrc using Pointers
  REAL, POINTER :: neighPhis(:, :, :)                !Neighbor phis for AFW decusping

  !BDF
  REAL, POINTER :: TranPhi2(:, :, :), TranPhi3(:, :, :), TranPhi4(:, :, :), TranPhi5(:, :, :)

  !SCM
  REAL, POINTER :: TranPrec(:, :, :)

  REAL, POINTER :: PhiSCM(:, :, :)
  REAL, POINTER :: PsiSCM(:, :)
  REAL, POINTER :: xsnfSCM(:, :, :)

  !AfSrc
  REAL, POINTER :: TranPhi1a(:, :, :, :, :), TranPhi2a(:, :, :, :, :)

  !AM3
  REAL, POINTER :: ResSrcD(:, :, :)
  !-----------------------------------------------------------------------------------------------------
ENDTYPE


TYPE CMInfo_Type

  REAL, POINTER :: PhiC(:, :, :), PhiFM(:, :, :)
  REAL, POINTER :: PsiC(:, :), PsicD(:, :)
  REAL, POINTER :: PsiFM(:, :), PsiFMD(:, :)
  REAL, POINTER :: RadJout(:, :, :, :, :)
  REAL, POINTER :: AxSrc(:, :, :), AxPXS(:, :, :)
  REAL, POINTER :: AxDtil(:, :, :, :), AxDhat(:, :, :, :), AxPDhat(:, :, :, :)
  TYPE(PinXS_Type), POINTER :: PinXS(:, :)
  TYPE(CMFDLS_TYPE), POINTER :: CoreCMFDLS(:)
  !TYPE(AXFLX_TYPE), POINTER :: AxFlx(:, :)
  TYPE(RayInfo4CMFD_Type), POINTER :: RayInfo4Cmfd

  TYPE(CMFDLS_TYPE), POINTER :: GcCMFDLS(:)
  !TYPE(CMFDLS_TYPE), POINTER :: CMFD2GLS(:)
  TYPE(PinXS_Type), POINTER :: GcPinXS(:, :)
  REAL, POINTER :: GcPsiC(:, :), GcPsicD(:, :)
  REAL, POINTER :: GcPhiC(:, :, :)
  REAL, POINTER :: phic_adj(:,:,:)

  !Transient Variables----------------------------------------------------------------------------------
  REAL, POINTER :: TranPhiCm(:, :, :), TranPhiFm(:, :, :)
  REAL, POINTER :: PrecCm(:, :, :)
  REAL, POINTER :: TranPsiCm(:, :), TranPsiCmd(:, :)
  REAL, POINTER :: PrecFm(:, :, :)
  REAL, POINTER :: TranPsiFm(:, :), TranPsiFmd(:, :)
  REAL, POINTER :: ResSrcCm(:, :, :)
  REAL, POINTER :: ResSrcFm(:, :, :)

  REAL, POINTER :: GcTranSrc(:, :, :)

  !Adaptive Theta
  REAL, POINTER :: ThetaCM(:, :, :)

  !BDF
  REAL, POINTER :: TranPhiCm2(:, :, :), TranPhiFm2(:, :, :)
  REAL, POINTER :: TranPhiCm3(:, :, :), TranPhiFm3(:, :, :)
  REAL, POINTER :: TranPhiCm4(:, :, :), TranPhiFm4(:, :, :)
  REAL, POINTER :: TranPhiCm5(:, :, :), TranPhiFm5(:, :, :)

  !SCM Variable
  REAL, POINTER :: PhiCSCM(:, :, :)
  REAL, POINTER :: PsiCSCM(:, :)

  !AM3
  REAL, POINTER :: ResSrcCmD(:, :, :)
  REAL, POINTER :: ResSrcFmD(:, :, :)
  !-----------------------------------------------------------------------------------------------------
END TYPE
!XS related TYPE

TYPE DcplInfo_Type

  LOGICAL :: lfeedback = .FALSE.
  !LOGICAL :: lReftemp = .FALSE.
  INTEGER :: nRefPln,nPln
  INTEGER :: nRefTemp = 1
  INTEGER, POINTER :: RefPln(:)
  INTEGER, POINTER :: Pln_Map(:)
  INTEGER, POINTER :: Pln_Map0(:)
  INTEGER, POINTER :: RefPln_Map(:, :)
  REAL, POINTER :: Ref_Temp(:)
  REAL :: eigv(10, 100) = 1._8
  TYPE(FmInfo_Type), POINTER :: DcplFmInfo(:, :)
  TYPE(CmInfo_Type), POINTER :: DcplCmInfo(:, :)
  TYPE(ThInfo_Type), POINTER :: DcplThInfo(:)
  TYPE(FxrInfo_Type), POINTER :: FXR(:, :)
  TYPE(MicXsFtn_Type), POINTER :: MicXsFtn(:, :)
  LOGICAL :: ldum1
  INTEGER :: idum1
ENDTYPE


TYPE DhatData_Type

  INTEGER :: ntemp, idum
  REAL :: dtil(0:30, 500)
  REAL :: dhat(0:30, 500)
END TYPE

TYPE CellXSData_Type

  LOGICAL :: lAlloc = .FALSE.
  INTEGER :: nFsr, nFxr
  INTEGER :: ntemp
  INTEGER :: nFsrMax, nTempMax

  TYPE(FxrInfo_Type), POINTER :: Fxr(:)
  REAL, POINTER :: phis(:, :, :)
  TYPE(CEll_Type), POINTER :: CellInfo
  TYPE(DhatData_Type) :: DhatData(4)

  LOGICAL :: ldum1
  INTEGER:: idum1
END TYPE

TYPE GroupInfo_Type

  LOGICAL :: lUpScat, ldum
  INTEGER :: ng, nofg, norg, nchi, ntiso, ntiso_depl, iresoGrp1, iresoGrp2
  INTEGER :: UpscatRange(2)
  INTEGER, POINTER :: InScatRange(:, :), OutScatRange(:, :)
  INTEGER :: nGC                 !Number of Condensed Group
  INTEGER :: GCStruct(2,300)     !Coarse Group => Fine Group
  INTEGER :: InvGCStruct(300) !FineGroup => Coarse Group

  INTEGER :: nprec = 6
! PHOTON Information...  |-- JSU EDIT 2019.05.08
  INTEGER :: ngg, nele
  LOGICAL :: lUpScat_Ph
  INTEGER :: UpScatRange_Ph(2)
  INTEGER, POINTER :: InScatRange_Ph(:, :), OutScatRange_Ph(:, :)
  ! AVERAGE (GEOMATRICAL AVERAGE) OF ENERGY IN THE GROUP                  |-- JSU EDIT 2017.09.13. |
  REAL, POINTER :: GamAvgE(:), NeuAvgE(:)   ! in Joule (1.602*1.e-19eV)
END TYPE


TYPE BenchXS_type

  LOGICAL :: lempty, lfuel, lcr
  REAL, pointer :: xst(:), xstr(:), xskf(:), xsnf(:), chi(:), xsa(:)      !Total, Transport, kaffa-fission nu-fission. CHI, absorptiom
  REAL, pointer :: xss(:,:), xss0(:,:), xss1(:,:), xss0t(:,:)
  REAL, pointer :: xss2(:,:), xss3(:,:)
  REAL, pointer :: xs0sum(:), xs1sum(:)
  TYPE(scatmat), pointer :: PnSM(:,:)
  REAL, POINTER, DIMENSION(:,:) :: rt1, rt2, w1, w2

  !Transient
  REAL :: wt = 0._8
  INTEGER :: iso0 = 0
  INTEGER :: iso1 = 0
  LOGICAL :: lCusping = .FALSE.

  !NEACRP
  REAL :: boronppm = 0.
END type

TYPE DynBenchXS_type

  LOGICAL :: lempty, lfuel, lcr
  REAL, pointer :: xst(:,:), xstr(:,:), xskf(:,:),xsnf(:,:),chi(:,:),xsa(:,:)      !Total, Transport, kaffa-fission nu-fission. CHI, absorptiom
  REAL, pointer :: xss(:,:,:), xss0(:,:,:), xss1(:,:,:), xss0t(:,:,:)
  REAL, pointer :: xss2(:,:,:), xss3(:,:,:)
  REAL, pointer :: xs0sum(:,:), xs1sum(:,:)
  TYPE(scatmat), pointer :: PnSM(:,:,:)
  REAL, POINTER, DIMENSION(:,:,:) :: rt1, rt2, w1, w2

  REAL, POINTER :: Beta(:,:)
  REAL, POINTER :: ChiD(:,:)
  REAL, POINTER :: Velo(:,:)
  REAL, POINTER :: Lambda(:,:)

  !Transient
  REAL :: wt = 0._8
  INTEGER :: iso0 = 0
  INTEGER :: iso1 = 0
  LOGICAL :: lCusping = .FALSE.
  INTEGER :: ntemp
  REAL, POINTER :: temp(:)

END TYPE

TYPE NeaCrpBenchXs_type

  INTEGER :: basexsl
  LOGICAL :: lCA = .FALSE.
  REAL :: c0, tm0, rho0, td0
  ! Follows the xs key in the neacrp benchmark
  REAL :: xs0(9)
  REAL :: gradxs_c(9), gradxs_tm(9), gradxs_rho(9), gradxs_td(9)
END TYPE

TYPE BenchKinParam_type

  REAL, pointer :: BETA(:)
  REAL, pointer :: Velo(:)
  REAL, pointer :: ChiD(:)
  REAL, pointer :: ChiDg(:,:)
END TYPE

TYPE XsMac_Type

  INTEGER :: id = 0
  logical :: lalloc = .FALSE.
  logical :: lfuel
  INTEGER :: ng = 0
  REAL, POINTER :: xsmact(:)               !Total Xs
  REAL, POINTER :: xsmactr(:)              !Transport
  REAL, POINTER :: xsmaca(:)               !Absorption
  REAL, POINTER :: xsmacf(:)               !Fission
  REAL, POINTER :: xsmacnf(:)              !Nu-Fission
  REAL, POINTER :: xsmackf(:)              !kappa-Fission
  REAL, POINTER :: xsmacs(:)               !Total Scattering Xs
  REAL, POINTER :: xsmacstr(:)             !Total Scattering Xs
  REAL, POINTER :: xsmacsm(:,:)            !Scattering Matrices
  LOGICAL :: lallocsm = .FALSE.
  REAL, POINTER :: xsmacp1sm(:,:)          !P1Scattering Matrices
  REAL, POINTER :: xsmacp2sm(:,:)          !P2Scattering Matrices
  REAL, POINTER :: xsmacp3sm(:,:)          !P3Scattering Matrices
  REAL, POINTER :: CHI(:)
  INTEGER :: niso
  Logical :: lIsoAlloc = .FALSE.
  REAL, POINTER :: IsoXsMacA(:, :), IsoXsMacTr(:, :), IsoXsMacT(:, :) !BYS edit 14/06/01
  REAL, POINTER :: IsoXsMacS0(:, :), IsoXsMacS1(:, :), IsoXsMacSS(:, :)
  REAL, POINTER :: ISoXsMacnf(:, :), IsoXsMacf(:, :), IsoXsMackf(:, :)
  REAL, POINTER :: IsoXsRadCap(:, :)                                                   !-- JSU EDIT 20170727
  REAL, POINTER :: phiIso(:,:)                                                         !PHS add  15/10/14
  REAL, POINTER :: isoxsmacsm(:,:,:)            !PHS add  28/12/15 isotope-wise scattering Matrices
! FOR EXPLICIT KAPPA CALCULATION                                        |-- JSU EDIT 2019/08/14
  LOGICAL :: lKERMAAlloc = .FALSE. !
  LOGICAL :: lISOKERMAAlloc = .FALSE.
  REAL, POINTER, DIMENSION(:)   :: MacKERMA_t, MacKERMA_s, MacKERMA_d, MacKERMA_p, MacKERMA_f
  REAL, POINTER, DIMENSION(:,:) :: IsoMacKERMA_t, IsoMacKERMA_s, IsoMacKERMA_d, IsoMacKERMA_p, IsoMacKERMA_f
  REAL, POINTER, DIMENSION(:)   :: MacDelkf
  REAL, POINTER, DIMENSION(:,:) :: IsoMacDelkf
! For FXR-wise Point-wise calculation
  LOGICAL :: lAllocPXS = .FALSE., lAllocPXSS = .FALSE.
  REAL, POINTER, DIMENSION(:)   :: PXST_MIX, PXSA_MIX, PXSF_MIX, PXSS_MIX, UFGFLX_MIX, PXSLS_MIX
  REAL, POINTER, DIMENSION(:)   :: PXST_ISO, PXSA_ISO, PXSF_ISO, PXSS_ISO, UFGFLX_ISO, PXSLS_ISO
END TYPE

TYPE DeplGd_TYPE

  INTEGER :: niso
  INTEGER, POINTER :: idiso(:)
  REAL :: phi1g, n155
  REAL, POINTER :: xsa(:), xsf(:), xsn2n(:), xsn3n(:)
END TYPE

TYPE DeplPCP_TYPE

  REAL ::pnum(2, -1:1, 64152:64160) = 0
  REAL :: f(64152:64160) = 1

END TYPE

TYPE Fxrinfo_type

  LOGICAL :: lfuel = .FALSE.
  LOGICAL :: lAIC = .FALSE.
  LOGICAL :: lCLD = .FALSE.
  LOGICAL :: ldepl = .FALSE.
  LOGICAL :: lres = .FALSE.
  LOGICAL :: lGD = .FALSE.
  LOGICAL :: lh2o = .FALSE.                                   !Fuel, depletion, resonnance
  LOGICAL :: lUse = .FALSE.
  LOGICAL :: lVoid = .FALSE.                                  ! Void
  INTEGER :: imix                                             ! Mixture
  INTEGER :: niso = 0                                         ! No. Isotope for Neutronics calculation
  INTEGER :: niso_Res = 0                                     ! No. Isotope for Neutronics calculation (Removing duplicated nuclides)
  INTEGER :: niso_depl = 0                                    ! No. Isotope for Depletion Calcualtion
  INTEGER :: niso_past = 0                                    ! No. Isotope of previous step for Depletion Calcualtion
  INTEGER :: FsrIdxSt, nFsrInFxr                              !# Global FSr Index # region
  REAL :: Burnup = 0, Burnup_past = 0, Hmkg0 = 0              !# Local Burnup(Mwd/kgU), Local Burnup past(Mwd/kgU), Initially loaded Heavy Metal(MgU)
  REAL :: temp, area, xstilde                                 !Temperature, Area
  INTEGER :: ndim = 0                                         ! Size of Array for Isotope Data
  INTEGER, POINTER :: idiso(:), idiso_past(:), idiso_pastpsm(:), idx_Res(:), idiso_Res(:) !Isotope Id List
  REAL, POINTER :: pnum(:), pnum_past(:), chi(:)              !Number density, Equip. XS,
  REAL, POINTER :: pnum_all(:), pnum_past_all(:)              !Depletion Bug Fix ! 16/02/11 Depletion timestep bug fixed
  REAL, POINTER :: pnum_Res(:)!, pnumrat(:,:)
  LOGICAL :: L_PNUM_ALL = .FALSE.                             ! 16/02/11 Depletion timestep bug fixed
  ! Xenon Under Relaxation Variables
  REAL :: pnXe = 0, pnI = 0, absXe = 0
  ! MLG SGFSP Variables
  REAL(4), POINTER :: xseq_f_1g(:), xseq_f_mg(:,:), xseq_c_1g(:) ! Escape(Equivalent) XS with MLG
  REAL(4), POINTER :: FnAdj(:), FtAdj(:,:)                       ! NDCF and TCF
  ! Categorization or Isotope-wise SGFSP
  REAL(4), POINTER :: xseq(:, :, :), NDAF(:, :, :)               ! Escape(Equivalent) XS and NDAF with categorization or isotope-wise
  ! Effective XS Ratio
  REAL(4), POINTER :: fresoa(:), fresof(:), fresos(:), fresostr(:), fresonf(:), fresokf(:) !Effective XS
  REAL(4), POINTER :: fresoAIso(:,:),fresoFIso(:,:),fresoSIso(:,:),fresoSSIso(:,:),fresoS1Iso(:,:)
  !
  REAL, POINTER :: DelInflow(:)
  REAL, POINTER :: Dcpl_Temp(:), Dcpl_pnum(:, :)
  REAL, POINTER :: Dcpl_fresoA(:, :), Dcpl_fresoF(:, :), Dcpl_fresoS(:, :), Dcpl_fresoStr(:, :)
  REAL(4), POINTER :: fresocapIso(:,:)                            !-- JSU EDIT 20170727  : Resonance Treatment for Isotopewise Radioactive Capture
  TYPE(DeplGd_TYPE), POINTER :: DeplXs1G(:)                   !Quadratic Interpolation of XS for GD pin
  TYPE(DeplPCP_TYPE), POINTER :: DeplPCP
  !Xenon Dynamics
  TYPE(XeDynFxr_Type), POINTER :: XeDyn
  !Kinetic Parameter
  LOGICAL :: lMixtureMixing = .FALSE.
  INTEGER :: imix0
  REAL ::  h2ofrac = 1.0_8
  REAL :: betat
  REAL :: siga, siga0, siga4g(4), siga4g0(4)
  REAL, POINTER :: beta(:), velo(:), veloh(:), chid(:), chip(:) !, chidk(:,:)
  !Control Rod Cusping

  LOGICAL :: lCrFxr = .FALSE.
  LOGICAL :: lCrCspFtn = .FALSE.
  LOGICAL :: lCrRes = .FALSE.
  !REAL, POINTER :: fCusping(:)
  TYPE(CspFXR_TYPE), POINTER :: CspFXR
  !Spectral SPH Factor related
  LOGICAL :: lSSPHalloc
  REAL(4),POINTER :: SPHfactor(:)
  INTEGER :: ipin

  !AFW Cusping
  LOGICAL :: lCusping = .FALSE.
  INTEGER :: iso0, iso1
  REAL :: wt = 0._8

  ! NeaCrp
  REAL :: DopTemp, rho
END TYPE

TYPE XeDynFxr_Type

  REAL ::  ARate_Xe135(2)
  REAL ::  Prate_Xe135(2)
  REAL ::  Prate_I135(2)
  REAL ::  Pn_Xe135(2)
  REAL ::  Pn_I135(2)
END TYPE

TYPE DcplFxrInfo_Type

  TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
  INTEGER :: nTemp
  REAL :: SubGrpRefTemp

  !INTEGER :: idum1
END TYPE


!Parallel Enviorment Variables

TYPE PE_TYPE

  LOGICAL :: master = .TRUE.
  LOGICAL :: slave = .FALSE.
  LOGICAL :: lidle = .FALSE.
  LOGICAL :: lDcplParallel = .FALSE.
  LOGICAL :: RadParallel = .FALSE.
  LOGICAL :: AxialParallel = .FALSE.         !Axial Parallel, Radial Parallel
  LOGICAL :: lAxSolParallel = .FALSE.        !Axial Solver Parallel
  LOGICAL :: RTMaster = .TRUE.
  LOGICAL :: RTSlave = .FALSE.
  LOGICAL :: CMFDMaster = .TRUE.
  LOGICAL :: CMFDSlave = .FALSE.
  LOGICAL :: lCmfdGrp = .TRUE.               !Participating a CMFD calculation or not
  LOGICAL :: lRTGrp = .TRUE.                 !Participating a RT calculation or not
  LOGICAL :: lThread = .FALSE.
  LOGICAL :: NonBlockPinXsComm = .FALSE.
  LOGICAL :: lUsrAxDcp = .FALSE.             !User defined Axial Decomposition
  INTEGER :: nproc, nproc0, myrank, myrank0
  INTEGER :: nThread = 1
  INTEGER :: nAxThread = 1                         !Threadid
  INTEGER :: nCmfdThread = 1
  INTEGER :: nDeplThread = 1
  INTEGER :: nCmfdGrp = 0
  INTEGER :: myzb, myze, myzbf, myzef, nxy, nz, nzfm, ng
  INTEGER :: nRay, idum
  INTEGER, POINTER :: RayIdx(:), idum1(:)
  INTEGER :: AxDomRange(2, 0:500), RadDomRange(2, 0:500)  !
  INTEGER :: SubPlnDomRange(2, 0:500)
  INTEGER :: UsrAxDecomp(0:500)
  INTEGER, POINTER :: AxDomList(:), RadDomList(:)
  INTEGER :: myRayBeg, myRayEnd
  INTEGER :: nMyRay
  INTEGER :: myNxyBeg, myNxyEnd
  INTEGER :: myOmpNxyBeg(100), myOmpNxyEnd(100)
  INTEGER :: myOmpFsrBeg(100), myOmpFsrEnd(100)
  INTEGER, POINTER :: RayList(:)
  INTEGER :: MPI_NTRACER_COMM, MPI_RT_COMM, MPI_CMFD_COMM, MPI_AXN_COMM, MPI_COMM
  INTEGER :: MPI_RTMASTER_COMM, MPI_NULL
  INTEGER :: nRtProc, nCmfdProc, nAxnProc
  INTEGER :: myRTrank, myCMFDrank

  INTEGER :: myRefPlnBeg, myRefPlnEnd
  INTEGER :: MPI_DCPLMASTER_COMM
  INTEGER :: RefPlnRange(2, 0:500)
  LOGICAL :: lmyRefPln(500) = .FALSE.

  INTEGER :: WorldGroup, LocalGroup
  LOGICAL :: lSubCh_proc=.FALSE.

  !--- CNJ Edit : Intel Math Kernel Library
  LOGICAL :: lMKL = .FALSE.

  !--- CNJ Edit : GPU Acceleration
  LOGICAL :: lCUDA = .FALSE., lCUDACMFD = .FALSE., lCUDATH = .FALSE., lCUDADepl = .FALSE.

  !--- CNJ Edit : Domain Decomposition + MPI
  INTEGER :: myAsyBeg, myAsyEnd
  INTEGER :: myPinBeg, myPinEnd
  INTEGER :: myFsrBeg, myFsrEnd
  INTEGER :: myFxrBeg, myFxrEnd
  INTEGER, POINTER :: nFsr(:), nPin(:), nAsy(:)
  INTEGER, POINTER :: Fsr_displs(:), Pin_displs(:), Asy_displs(:)
END TYPE

!Power Distribution

TYPE PowerDist_Type

  REAL :: fxyzn, fxyzp, fxyn, fxyp, frn, frp, fz, fxavg, fxmax
  REAL :: Pin3DNormalizer, Pin2DNormalizer, Fm2DNormalizer, Fm3DNormalizer
  REAL :: Asy3DNormalizer, Asy2DNormalizer
  REAL :: Axial1DNormalizer
  REAL, POINTER :: PinPower2D(:, :), PinPower3D(:, :, :)
  REAL, POINTER :: AsyPower2D(:), AsyPower3D(:, :)
  REAL, POINTER :: Axial1DPower(:)
  REAL :: pwsum
END TYPE

TYPE DancoffDist_Type

  REAL, POINTER :: Dancoff3D(:, :, :)
END TYPE
!Mixture Information

TYPE Mixture_Type

  LOGICAL :: lempty = .TRUE.
  LOGICAL :: lres = .FALSE.
  LOGICAL :: lfuel = .FALSE.
  LOGICAL :: ldepl = .FALSE.
  LOGICAL :: lh2o = .FALSE.
  LOGICAL :: lCLD = .FALSE.
  LOGICAL :: lGd = .FALSE.
  LOGICAL :: lMox = .FALSE.
  LOGICAL :: lAIC = .FALSE.
  LOGICAL :: lcompact = .FALSE.
  LOGICAL, POINTER :: lSSPH(:)
  CHARACTER(4) :: name
  INTEGER :: niso, ntriso, imatrix
  REAL :: temp, dens, fv0
  REAL :: h2ofrac0 = 1.0
  INTEGER, POINTER :: idiso(:), itriso(:)
  REAL, POINTER :: fweig(:), pnum(:), ftriso(:)
  INTEGER :: deplopt
END TYPE

!T/H Condition
TYPE THInfo_TYPE

  REAL, POINTER :: RefFuelTemp(:)
  REAL, POINTER :: FxrTemp(:, :)
  REAL, POINTER :: RelPower(:, :)
  REAL, POINTER :: Tdop(:, :)
  REAL, POINTER :: Tcool(:, :)
  REAL, POINTER :: TcoolInOut(:, :, :)
  REAL, POINTER :: Tfvol(:, :, :)
  REAL, POINTER :: DenCool(:, :)
  REAL, POINTER :: qvol(:, :)
  REAL, POINTER :: UserDefFuelT(:, :)
  REAL, POINTER :: UserDefModT(:, :)
  REAL, POINTER :: CBMCool(:,:)

  REAL :: PowFa                     ! Nomial Asy Power at the Full Power Condi, (watt)
  REAL :: Pexit                     ! Core Exit Coolant Pressure (Pa)
  REAL :: PowLin                    ! Linear Power Density
  REAL :: PowLv = 1.0_8             ! initial power level
  REAL :: Tin                       ! inlet Temp           : tin
  REAL :: Denin                     ! Inlet Density        : din
  REAL :: TDopIn                    ! Inlet Doppler Temp   : tdopin
  REAL :: NomCoolT                  ! Nomial Coolant Temp  : tcnom
  REAL :: NomFuelT                  ! Nomial Fuel Temp     : tfnom
  REAL :: MdotFA                  !
  REAL :: tfmax                     ! Maxium Fuel Center Line Temperature
  REAL :: tfcl_max, tfavg_max
  REAL :: TdopChg                   ! Doppler Temperature Change Max.
  REAL :: TModoutAvg                ! Average Moderator Temperature
  REAL :: pwpeakf                   ! powr peaking factor
  TYPE(FuelTh_Type), POINTER :: FuelTh(:)
  TYPE(CoolantTH_Type), POINTER :: CoolantTh(:)

  REAL :: Rhou
  REAL, POINTER :: TdopBase(:,:)
END TYPE

TYPE FuelTh_Type

  LOGICAL :: lFuel = .FALSE.
  LOGICAL, POINTER :: lMox(:)
  REAL, POINTER :: hflux(:)
  REAL, POINTER :: qvol(:)
  REAL, POINTER :: Tcool(:)
  REAL, POINTER :: htcoef(:)             ! Heat Tranfer Coefficient
  REAL, POINTER :: tfuel(:, :)           ! Tfuel(nr,nz)
  REAL, POINTER :: tfvol(:, :)
  REAL, POINTER :: teff(:)
  REAL :: TdopMax !

  REAL, POINTER :: tfueld(:, :), qvold(:), qshaped(:)
END TYPE

TYPE CoolantTH_Type

  LOGICAL :: lFuel = .FALSE.
  REAL, POINTER :: hcool(:)              ! volume enthalpy in J/Kg
  REAL, POINTER :: rhou(:)               ! junction mass flux in Kg/m^2-sec
  REAL, POINTER :: rhohu(:)              ! junction enthaply flux in J/m^2-s
  REAL, POINTER :: u(:), ud(:)           ! junction velocity in m/s
  REAL, POINTER :: qeff(:)               ! effective volumetric HGR in cool.
  REAL, POINTER :: qvol(:)               ! Pointing Varialbe
  REAL, POINTER :: Tcool(:)              ! Pointing Varialbe
  REAL, POINTER :: TcoolInOut(:, :)      !
  REAL, POINTER :: DenCool(:)

  !Variable for transient
  REAL, POINTER :: TCoold(:)
  REAL, POINTER :: DenCoold(:)           ! Coolant Density at the previous Time Step
  REAL, POINTER :: hcoold(:)             ! Volume enthalpy at the previous Time Step in J/kg
  REAL, POINTER :: rhoud(:)              ! Junnction mass flux
  REAL, POINTER :: rhohud(:)
  REAL, POINTER :: qeffd(:)              !

  !Variable for Assembly Dependent acf
  REAL :: acf, zetap, Deq
END TYPE


TYPE THOpt_Type

  REAL :: KFuelCorrelation(0:5)     !Fuel Conduction Correlation
  REAL :: KCladCorrelation(0:3)
  REAL :: CpFuelCorrelation(0:3)
  REAL :: CpCladCorrelation(0:3)
  REAL :: Hgap = 10000._8
  REAL :: PinDim(5)                 !  1 : PinRadius, Caldding Outer, Clad Thickness, Guide Tube Radius
  INTEGER :: nrpellet  = 20
  INTEGER :: nCondRing = 9

  INTEGER :: KFuelModel = 1            !  1: NEACRP-l-335  2: FRAPCON
  INTEGER :: KCladModel = 1            !  1: Zirc          2: SS
  INTEGER :: CpFuelModel = 1           !  1: NEACRP-l-335  2: FRAPCON
  INTEGER :: CpCladModel = 1           !  1: Zirc-NEACRP-l-335  2: SS 3: Zirc-FRAPCON
  INTEGER :: hGapModel = 1             !  1: Fixed Value   2: Interpolation  3: Modified Ross Stoute (Campbell)
  REAL :: RhoFuel = 10.282
  REAL :: RhoClad = 6.6
  LOGICAL :: AA_STH = .FALSE.
END TYPE

TYPE THVar_Type

  REAL :: rs, rw, tw, rgap, rgto, rgti           !  Pallet R, Clad R, Thickness of Cladding, GT Radius, m unit
  REAL :: ChannelPitch              !  Pitch of Channel, m unit
  REAL :: AsyPitch                  !  Pitch of Asy, m unit
  REAL :: acf                       !  Fuel Channel flow area
  REAL :: afp                       !  Fuel pellet area
  REAL :: acfgt                     !  GT channel Flow Area
  REAL :: xi                        !  Wetted Perimeter
  REAL :: zeta                      !  Heated Perimeter
  REAL :: zetap                     !  Heated Perimeter Density
  REAL :: deq                       !  Equivalence Diameter
  REAL :: hact                      !  Active Height
  REAL :: FuelDelR, CldDelR         !  Radial mesh size for fuel and clad
  REAL :: FracDC = 0.               !  fraction of heat deposited directly in coolant
  REAL :: Fracdf                    !

  REAL :: Tin                       ! inlet Temp           : tin
  REAL :: MdotFA                    ! mass flow rate per channel
  REAL :: Pexit

  REAL :: DelT                      !Current Time Step Size

  REAL, POINTER :: Hz(:), r(:), ChanVol(:)
  REAL :: BoilingTemp


  !INTEGER :: nAsyCh                 !  # channel per asy
  !INTEGER :: nAsyGT                 !  # GT per Asy
  REAL :: nAsyCh                 !  # channel per asy
  REAL :: nAsyGT                 !  # GT per Asy
  INTEGER :: nzth                   !
  INTEGER :: npr, npr1, npr2, npr3, npr4, npr5
  LOGICAL :: lhact = .FALSE.

  INTEGER :: nChType
  TYPE(THCh_Type), POINTER :: THCh(:)

END TYPE

TYPE THCh_Type
  REAL :: nAsyCh, nAsyGT
  REAL :: ChannelPitch, AsyPitch, hact
  REAL :: acf, xi, zetap, Deq
END TYPE

TYPE THCell_Type

  INTEGER :: CldReg(2), FuelReg(2), CoolReg(2)
  INTEGER, POINTER :: FuelMapping(:, :)
  REAL, POINTER :: Frac(:, :)
  INTEGER, POINTER :: CondMapping(:, :)
  REAL, POINTER :: InvFrac(:, :)
END TYPE

TYPE CrCell_TYPE

  LOGICAL :: LCrIn = .FALSE.
  LOGICAL :: LCrOut = .FALSE.
  INTEGER :: CellCrIdx
  INTEGER :: inIdx, outIdx
  INTEGER :: nCrFxr, nCrFsr                 ! # of Control Rod Fxr, Fsr
  INTEGER, POINTER :: CrFxrIdx(:)           ! Control Rod Fxr Index
  INTEGER, POINTER :: CrFsrIdx(:)           ! Control ROd Fsr Index
  INTEGER, POINTER :: imix_crin(:)          ! Control Rod Material
  INTEGER, POINTER :: imix_crout(:)         ! Control Rod Material
END TYPE

TYPE CrAsyConf_TYPE

  INTEGER :: nCrPin
  INTEGER, POINTER :: CrLoc(:, :)
  INTEGER, POINTER :: CrPinIdx(:)
ENDTYPE

TYPE CrBank_TYPE
  LOGICAL :: lCrDecusp = .FALSE.
  LOGICAL :: lCrDir = .TRUE.
  LOGICAL :: lCrIn = .FALSE.
  LOGICAL :: lCrFullIn = .FALSE.
  LOGICAL :: lPosInp = .FALSE.      !Position is specified or not
  LOGICAL :: lCrSearch = .FALSE.    !Control Rod Moving or Not
  INTEGER :: PosType = 2            !1 : specify step input, 2: specify control rod location directly
  CHARACTER(10) :: BankName
  REAL :: RodInPos = 0              !Fully Inserted Position
  REAL :: RodStepSize = 0
  REAL :: RodPos = 0                !Rod Position
  INTEGER :: RodStep = 0
  INTEGER :: nCrAsy = 0
  INTEGER, POINTER :: CrAsyLoc(:)
  REAL, POINTER :: FracCrIn(:)       !Ratio of the Inserted region for each plane(1: fully inserted, 0 No insertion)
  INTEGER :: isonew = 0              !Mixture index for Benchmark XS
END TYPE

TYPE CrPosDat_Type

  INTEGER :: nDat = 0
  LOGICAL :: LExist(100) = .FALSE.
  REAL :: PosInpDat(100, 100)
  REAL :: RodPos(100, 100)
  INTEGER :: RodStep(100, 100)

END TYPE

TYPE CrSearch_Type

  LOGICAL :: lInit = .FALSE.
  LOGICAL :: lFullIn = .FALSE.
  LOGICAL :: lFullOut = .FALSE.
  LOGICAL :: lSearchConv = .FALSE.
  INTEGER :: CRMV_DOM(2)=0
  INTEGER :: iter = 0
  INTEGER :: eigv_iter = 0
  REAL :: CrMv = 0       !Control Rod Movement
  REAL :: CrMvd = 0      !Control Rod Movement of previous iteration
  REAL :: CrWorth = 10.     !Cr Worth
  REAL :: CrWorthd = 10.    !Cr Worth of previous movement
  REAL :: eigvd             !Eigenvalue of Previous
END TYPE

TYPE CrInfo_Type

  LOGICAL :: lCrInfo = .TRUE.
  LOGICAL :: lCrChg = .FALSE.
  INTEGER :: nCrCell = 0
  INTEGER :: nCrAsyConf = 0
  INTEGER :: nCrBank = 0

  INTEGER :: nCrPin = 0
  INTEGER :: nCrAsy = 0

  INTEGER, POINTER :: CRCellMAP(:, :)
  TYPE(CrAsyConf_Type), POINTER :: CrAsyConf(:)
  TYPE(CrBank_Type), POINTER :: CrBank(:)

  TYPE(CrSearch_Type) :: CrSearch
  TYPE(CrPosDat_Type), POINTER :: CrPosDat
END TYPE

!LP_Suffling Information
TYPE LpShf_Type

  INTEGER :: nrstfile = 0
  LOGICAL :: lDataExist(200) = .FALSE.
  CHARACTER(256) :: RstFiles(200)
  REAL :: PUL(200) = 0                   ! Cooling Time
  LOGICAL :: lPUL = .FALSE.
  REAL :: HM_Mass0(200) = 0
  INTEGER :: cycleid                  ! Current Cycle No
  INTEGER :: RstCycleId
  INTEGER :: IdxSt(2) = 1
  LOGICAL :: lPLN_MAP = .FALSE.
  INTEGER :: PlN_MAP(200)
  TYPE(Shf_Type), POINTER :: Shf(:)   ! Shuffling Information

END TYPE

TYPE Shf_Type

  CHARACTER(80) :: fields
  LOGICAL :: LFreshFuel = .TRUE.
  LOGICAL :: LNoneFuel = .FALSE.
  LOGICAL :: lRead =.FALSE.
  INTEGER :: cycleid           !For fresh fuel
  INTEGER :: ix, iy, irot      !Assembly Location at the
END TYPE

TYPE XsChange_TYPE

  LOGICAL :: lUse = .FALSE.
  INTEGER :: Iso0, iso1
  INTEGER :: isonew = 0
  INTEGER :: iStepBeg, iStepend
  REAL :: tbeg, tend
  LOGICAL :: lStart = .FALSE.
  LOGICAL :: lComplete = .FALSE.
  LOGICAL :: lsupl = .FALSE.   !Supplementary option
  CHARACTER(256) :: field1, field2
  INTEGER :: izbeg, izend, nasy
  INTEGER, POINTER :: AsyList(:)
  REAL :: wt = 0
  LOGICAL :: lCusping = .FALSE.
  LOGICAL :: lCuspingDirection = .TRUE. !TRUE for Insertion, False for Withdrawal
  LOGICAL :: lStepFunc = .FALSE.
END TYPE

TYPE XsCntlRod_type

  INTEGER :: iso0, iso1
  INTEGER :: isonew = 0
  INTEGER :: izbeg, izend, nasy
  INTEGER, POINTER :: AsyList(:)
  CHARACTER(256) :: field1, field2
  REAL :: wt = 0
  LOGICAL :: lCusping = .FALSE.
  LOGICAL :: lCuspingDirection = .TRUE. !TRUE for Insertion, False for Withdrawal
  LOGICAL :: lStepFunc = .FALSE.
END TYPE

TYPE XSNoise_TYPE

  INTEGER :: itype
  REAL :: amp, freq, phase
  INTEGER :: ixa, iya, ixya
  INTEGER :: ix, iy, ixy
  INTEGER :: izbeg, izend
  LOGICAL :: lfirst = .TRUE.
  INTEGER :: iso0
  INTEGER :: isonew = 0
END TYPE

TYPE TranCntl_TYPE

  INTEGER :: nchange = 0
  TYPE(XsChange_TYPE), POINTER :: XsChange(:)
  TYPE(XsNoise_TYPE), POINTER :: XsNoise(:)
  TYPE(XsCntlRod_TYPE), POINTER :: XsCntlRod(:)
  REAL :: Tend, DelT0
  REAL :: Tstep_inp(0:100), Tdiv_inp(0:100)
  REAL :: Tdiv_corrector
  REAL :: dtth
  INTEGER :: nthstep

  REAL :: freq_inp(1:100)
  REAL :: freq_now
  INTEGER :: nfreq
  INTEGER :: nNoise = 0

  LOGICAL :: lfixtmprw = .FALSE.               ! rod worth calculation with fixed temperature condition

  LOGICAL :: lchidk = .FALSE.
  LOGICAL :: lExpTrsf = .FALSE.      !Exponential Transform
  LOGICAL :: lExpMOC = .FALSE.       !Exponential Source Term Control
  LOGICAL :: lStepFunc = .FALSE.
  LOGICAL :: lStepApprox = .FALSE.
  LOGICAL :: lCusping = .FALSE.
  LOGICAL :: lStepImplicit = .FALSE.
  LOGICAL :: ImplicitSwitch = .FALSE.
  LOGICAL :: lCorrector = .FALSE.
  LOGICAL :: lGuess = .FALSE.
  LOGICAL :: lPCQSIter = .FALSE.
  LOGICAL :: lIQS = .FALSE.
  LOGICAL :: lAdptT = .FALSE.
  INTEGER :: PCRtype = 1
  REAL :: cmfdres
  REAL :: delpsifm

  LOGICAL :: lTheta
  REAL :: Theta = 0.5
  REAL :: Theta0

  LOGICAL :: lMocUpdt = .TRUE.
  LOGICAL :: lAxNUpdt = .TRUE.
  LOGICAL :: lSGFSPUpdt = .FALSE.
  LOGICAL :: lCondiMOC = .FALSE.
  LOGICAL :: lXsPerturb = .FALSE.
  LOGICAL :: lKineticBen = .FALSE.
  LOGICAL :: lDynamicBen = .FALSE.
  LOGICAL :: lMethod = .FALSE.
  INTEGER :: it_woSG = 0
  INTEGER :: it_woMOC = 0

  INTEGER :: nowstep = 0
  INTEGER :: nstep
  REAL :: T(50000), DelT(50000)

  INTEGER :: nTWriteOut = 0
  REAL :: TWriteOut(4000)
  INTEGER :: StepWriteOut(4000)
  INTEGER :: NowWriteStep = 0
  INTEGER :: nMaxOuter = 15
  INTEGER :: nMaxCMFD = 50
  REAL :: res_conv = 1.e-5
  REAL :: res_conv2 = 1.e-6
  REAL :: psi_conv = 1.e-6
  REAL :: cmfd_res_conv = 1.e-7

  !BDF
  LOGICAL :: Cmfd_Bdf = .FALSE.
  LOGICAL :: MOC_Bdf = .FALSE.
  INTEGER :: BDF_ORDER = 5
  LOGICAL :: lCN_Step = .FALSE.
  INTEGER :: CN_Step
  REAL :: Coeff(0:5)

  !AT
  LOGICAL :: lAdpTheta
  LOGICAL :: lAdpThetaSt = .FALSE.
  LOGICAL :: AdpThetaMethod = .FALSE.
  INTEGER :: AdpThetaStep
  REAL, POINTER :: ThetaCM(:, :, :)

  !SCM
  LOGICAL :: lSCM
  LOGICAL :: lSCM_Prec = .FALSE.
  LOGICAL :: lAmpFrqFirst = .TRUE.
  REAL :: AmpFrq = 0._8, AmpFrqd = 0._8, AMpFrqdd = 0._8
  REAL :: AmpFrq1 = 0._8, AmpFrq2 = 0._8
  REAL :: AvgAmpFrq = 0._8
  REAL :: EigD = 1._8, EigD1 = 1._8, EigD2 = 1._8
  REAL, POINTER :: ShpFrqCM(:, :, :), ShpFrqCMd(:, :, :)
  REAL, POINTER :: ShpFrqFM(:, :, :), ShpFrqFMd(:, :, :)
  REAL, POINTER :: AvgShpFrqCM(:, :, :), AvgShpFrqFM(:, :, :)
  REAL, POINTER :: PrecFrqCM(:, :, :), PrecFrqFM(:, :, :)
  REAL, POINTER :: ShpFrqCM_2nd(:, :, :), ShpFrqFm_2nd(:, :, :)
  REAL, POINTER :: ShpFrqCMd_2nd(:, :, :), ShpFrqFmd_2nd(:, :, :)
  REAL, POINTER :: AvgShpFrqCM_2nd(:, :, :), AvgShpFrqFM_2nd(:, :, :)

  !AM3
  LOGICAL :: lAM = .FALSE.
  REAL :: Beta0, Gam1, Gam2

  !Method
  INTEGER :: TD !Temporal Discretization
  LOGICAL :: lblockGS =.TRUE.

  LOGICAL :: lIQSAA = .FALSE.
  REAL, POINTER :: IQSAA_x(:,:), IQSAA_g(:,:)
  INTEGER :: IQSAA_m

  !Noise Sampling
  REAL :: Speriod, Sbeg, Send !Sampling Period, Sampling Begin, Sampling End
  LOGICAL :: lNNSampling
  REAL, POINTER :: Ssteps(:)
  INTEGER :: nSstep
  INTEGER :: iSstep = 1
  INTEGER :: Sio = 81
END TYPE

TYPE TranInfo_TYPE


  LOGICAL :: lBenchXS = .TRUE.
  LOGICAL :: lLibXs = .FALSE.

  INTEGER :: nPrec =6
  REAL, POINTER :: Chid(:)           ! Delayed Neutron
  REAL, POINTER :: chidk(:,:)
  !REAL, POINTER :: Chia(:)           ! Average Fission Spectrum
  !REAL, POINTER :: Chip(:)           ! Prompt Fission Spectrum
  REAL, POINTER :: lambda(:)         ! Decay constant of Precursor
  REAL, POINTER :: Invlambda(:)
  !REAL, POINTER :: beta(:)           ! Delayed Neutron Fraction
  !REAL, POINTER :: neut_velo(:)         ! Neutron Velocity

  !REAL, POINTER :: Fsrbetat(:, :) !Total Delayed Neutron
  REAL, POINTER :: CellOmegam(:, :, :), Cellomega0(:, :, :), Cellomegap(:, :, :)
  REAL, POINTER :: FxrOmegam(:, :, :), FxrOmega0(:, :, :), Fxromegap(:, :, :)

  REAL, POINTER :: Expo(:, :, :), Expo_alpha(:, :, :)
  REAL, POINTER :: FmExpo(:, :, :), FmExpo_alpha(:, :, :)
  REAL, POINTER :: CorePower_History(:)

  REAL, POINTER :: PhiShape0(:, :), PhiShape(:, :)

  !REAL, POINTER :: Kappa(:)   !Exp(-lambda * DelT)

  REAL :: eigv0 = 1.0_8               !Initial Eigenvalue
  REAL :: PowerLevel0 = 1.0_8
  REAL :: PowerLevel = 1.0
  REAL :: fnorm =1.0_8                ! normalize flux such that average group flux be unity
  REAL :: UnitPowerLevel0 = 1.0       ! Power Level for normalized flux level at the initial condition of transeint
  REAL :: reactivity = 0.             ! Current Time step Reacitivy
  REAL :: reactivity_dynamic = 0.
  REAL :: corebeta = 0.               ! Current Time step coreavgBeta
  REAL :: corebeta_dynamic = 0.
  REAL :: corebeta_dynamic_MOC = 0._8
  REAL :: TranEig = 0.                ! Current Time step Estimated Eig
  REAL :: TranEig_dynamic = 0.
  REAL :: delrho = 0.                 ! Current Time step Estimated delrho
  REAL :: delrho_dynamic = 0.
  REAL :: lifetime = 0.               ! Current Time step neutron LifeTime
  REAL :: lifetime_Dynamic = 0.
  REAL :: lifetime_Dynamic_MOC = 0.
  REAL :: factor_F = 0.

  REAL :: Inv_Factor_F0
  REAL :: Inv_Factor_K0
  REAL :: Inv_lifetime0
  REAL :: Prev_delrho
  REAL :: Prev_corebetat = 0.
  REAL :: Prev_lifetime = 0.
  REAL :: Prev_Factor_F = 0.
  REAL, POINTER :: Prev_corePrec(:)
  REAL, POINTER :: Prev_coreBeta(:), coreBetak(:)
  REAL :: PwSum0
  INTEGER :: nfuelcell

  REAL :: AmpPredictor = 1.
  REAL :: Amp = 1.
  REAL :: TranAmp = 1.
  REAL :: TranAmpd = 1.
  REAL :: AmpRatio = 1.
  REAL :: AmpTilt = 0.
  REAL :: AmpTiltd = 0.
  REAL :: PrecRatio(6) = 1.

  REAL :: PowerLeveld = 1.0_8
  REAL :: rtdblr = 0                  !Inverse of Doubling Time

  REAL, POINTER :: RefTemp0(:), RefTemp(:)           !Reference Temperature
  REAL :: Tfcl_ref_SG, Tfcl_ref_MOC

  !-- Dynamic Benchmark (C5G7-TD Phase II)
  REAL :: InitTemp
  REAL :: InitPow = 1.e-6
  REAL, POINTER :: fuelTemp(:,:)
END TYPE
!
TYPE CspFXR_TYPE

  INTEGER :: niso(2)
  REAL, POINTER :: pnum(:, :)
  INTEGER, POINTER :: isolist(:, :)
  REAL, POINTER :: fcsp(:, :)
  INTEGER, POINTER :: isomap(:, :)
  REAL :: vwt(2) = 0.0
  LOGICAL :: lH2o(2) = .FALSE.
  REAL :: H2OFrac(2) = 1.0
END TYPE

TYPE MiscStruct_Type

  !Core Ring Structure
  LOGICAL :: lRing = .FALSE.
  LOGICAL :: lBarrel = .FALSE.
  INTEGER :: nring = 0
  INTEGER :: max_nring = 100
  REAL :: rad_ring(2, 0:100) = 0
  REAL :: ring_cent(2)
  REAL :: ring_plnbeg(0:100) = 0
  REAL :: ring_plnend(0:100) = 0
  INTEGER :: mix_ring(0:100) = 0
  REAL :: VolCor_ring(0:100)
  LOGICAL :: lVolCor_ring(0:100)
  !Baffle
  LOGICAL :: lBaffle = .FALSE.
  INTEGER :: Mix_Baffle               !
  REAL :: Thick_Baffle                !Baffle
END TYPE

TYPE XeDynState_Type

  INTEGER :: istep = 0
  REAL :: T, delT
  REAL :: PowLv = 1._8
  REAL :: FlowLv = 1._8
  LOGICAL :: lCntlRod = .FALSE.
  INTEGER :: CntlRodPos = 0
  LOGICAL :: lXeEq = .FALSE.
  LOGICAL :: lBoronSearch = .FALSE.
  LOGICAL :: lBoronVal = .FALSE.
  LOGICAL :: lPrevBoron = .FALSE.
  REAL :: Target_keff = 1.0_8
  REAL :: BoronPPM = 0
  CHARACTER(256) :: Inpline=''
END TYPE

TYPE XeDynInfo_Type

  LOGICAL :: lCalculation = .FALSE.
  INTEGER :: nTimeStep = 0
  INTEGER :: nState = 0
  INTEGER :: nMaxTimeStep = 1000
  REAL :: Tsec(0:1000), Tmin(0:1000), Th(0:1000), Tday(0:1000)
  INTEGER :: StateMapping(0:1000)
  LOGICAL :: lPerturb(0:1000)
  CHARACTER(10) :: Unit
  TYPE(XeDynState_Type) :: CoreState(0:1000)
END TYPE

END module

