MODULE PARAM
! CONSTANTS
INTEGER, PARAMETER :: XDIR=1,YDIR=2,ZDIR=3
INTEGER, parameter :: FULL = 360, QUARTER = 90
INTEGER, parameter :: NG2 = 2      ! 2. the number of energy groups for coarse group.
INTEGER, parameter :: NDIRMAX = 3  ! the number of maximum dimension, so this is
INTEGER, parameter :: FAST=1, THERMAL=2
INTEGER, parameter :: PLUS=1, MINUS=-1
logical, parameter :: TRUE=.true.,FALSE=.false.
REAL   , parameter :: CINF=1.0E30,ZERO=0.,ONE=1.,               &
                    RTHREE=1.0/3.0,HALF=0.5,BIG=1.0E30,         &
                    RFOUR=1.0/4.0,RFIVE=1.0/5.0,RSIX=1.0/6.0,   &
                    RSEVEN=1.0/7.0,R10=1.0/10.0,R14=1.0/14.0
REAL   , parameter :: CKELVIN=273.15_8, PI=3.14159265358979_8, AVOGADRO=0.6022137_8, INVPI=0.3183098861837907_8, hpi=1.570796326794897_8
INTEGER, PARAMETER :: CENTER=0, LEFT=1, RIGHT=2, BOTTOM=1, TOP=2
REAL, PARAMETER :: awh2o=18.01228_8, awboron=10.8120002746582_8 !, b10frac=0.198_8   !water and boron atomic weight
REAL, PARAMETER :: sigpH=20.4780_8, sigpO=3.8883_8, sigpB10=2.1424_8, sigpB11=4.8400_8
!INTEGER, PARAMETER :: WEST=1,EAST=2,NORTH=3,SOUTH=4

INTEGER, PARAMETER :: forward = 1, backward = 2   !--- CNJ Edit : Domain Decomposition
INTEGER, PARAMETER :: red = 1, black = 2   !--- CNJ Edit : Red-Black Domain Decomposition
INTEGER, PARAMETER :: SOUTH=1, WEST=2, NORTH=3, EAST = 4
INTEGER, PARAMETER :: PREV=-1, CURR=0, NEXT=1
INTEGER, PARAMETER :: VoidCell = 0, RefCell = -1, RotCell = -2, CbdCell = -3


REAL, PARAMETER :: epsm1=1.e-1_8, epsm2=1.e-2_8, epsm3=1.e-3_8, epsm4=1.e-4_8
REAL, PARAMETER :: epsm5=1.e-5_8, epsm6=1.e-6_8, epsm7=1.e-7_8, epsm8=1.e-8_8
REAL, PARAMETER :: epsm10=1.e-10_8, epsm20=1.e-20_8, epsm30=1.e-30_8

INTEGER, PARAMETER :: lsseigv=1, ldcplsseigv=2, ldepletion=3, lTransient = 4 !lshape=2, lreactcoef=3, ladjoint=4, ltransient=5,
INTEGER, PARAMETER :: lCrCspGen = 5, lXenonDynamics = 6, lBranch = 7, lEFTsearch =  8
INTEGER, PARAMETER :: lP1SENM = 1, lP3SENM = 2

INTEGER, PARAMETER :: MaxAxialPlane = 500
INTEGER, PARAMETER :: MaxPrec = 100
INTEGER, PARAMETER :: ngmax = 2000
INTEGER, PARAMETER :: nzmax = 500
INTEGER, PARAMETER :: nMaxFsr = 500
INTEGER, PARAMETER :: nCellXMax = 200, nCellMax = 2500
INTEGER, PARAMETER :: nThreadMax = 48
character(20), parameter :: AxSolverName(2) = (/'P1 SENM ',  'SP3 SENM'/)
character(1),parameter :: DOT='.',BANG='!',BLANK=' ', SLASH='/',AST='*', POUND='#'
character(512), parameter :: BLANK0=' '
character(126),parameter :: hbar1 = &
'------------------------------------------------------------------------------------------------------------------------------'
character(126),parameter :: hbar2 = &
'=============================================================================================================================='
CHARACTER(132) MESG
!DATA AxSolverName /'P1 SENM',  'SP3 SENM'/

REAL, PARAMETER :: eVtoJ=1.602*1.e-19

END MODULE