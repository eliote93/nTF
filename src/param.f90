MODULE PARAM
  
! CONSTANTS
INTEGER, PARAMETER :: XDIR = 1, YDIR = 2, ZDIR = 3
INTEGER, PARAMETER :: NG2 = 2
INTEGER, PARAMETER :: FAST = 1, THERMAL = 2
INTEGER, PARAMETER :: PLUS = 1, MINUS = -1
INTEGER, PARAMETER :: CENTER = 0, LEFT = 1, RIGHT = 2, BOTTOM = 1, TOP = 2
INTEGER, PARAMETER :: FORWARD = 1, BACKWARD = 2
INTEGER, PARAMETER :: RED = 1, BLACK = 2, BLUE = 3
INTEGER, PARAMETER :: SOUTH = 1, WEST = 2, NORTH = 3, EAST = 4
INTEGER, PARAMETER :: PREV = -1, CURR = 0, NEXT = 1
INTEGER, PARAMETER :: VoidCell = 0, RefCell = -1, RotCell = -2, CbdCell = -3
INTEGER, PARAMETER :: lsseigv = 1, ldcplsseigv = 2, ldepletion = 3, lTransient = 4, lCrCspGen = 5, lXenonDynamics = 6, lBranch = 7, lEFTsearch =  8, lNNFSP = 9
INTEGER, PARAMETER :: lP1SENM = 1, lP3SENM = 2
INTEGER, PARAMETER :: MaxPrec = 100, ngmax = 2000, nzmax = 500, nMaxFsr = 500, nCellXMax = 200, nCellMax = 2500, nThreadMax = 48

! LOGICAL
LOGICAL, PARAMETER :: TRUE = .true., FALSE = .false.

! REAL
REAL, PARAMETER :: ZERO = 0., ONE = 1.
REAL, PARAMETER :: RTHREE = 1.0/3.0, HALF = 0.5, BIG = 1.0E30
REAL, PARAMETER :: RFOUR  = 1.0/4.0, RFIVE =1.0/5.0, RSIX = 1.0/6.0
REAL, PARAMETER :: RSEVEN = 1.0/7.0, R10 = 1.0/10.0
REAL, PARAMETER :: CKELVIN = 273.15_8, PI = 3.14159265358979_8, AVOGADRO = 0.6022137_8, INVPI = 0.3183098861837907_8, HPI = 1.570796326794897_8
REAL, PARAMETER :: awh2o = 18.01228_8, awboron = 10.8120002746582_8
REAL, PARAMETER :: sigpH = 20.4780_8, sigpO = 3.8883_8, sigpB10 = 2.1424_8, sigpB11 = 4.8400_8
REAL, PARAMETER :: epsm1 = 1.e-1_8, epsm2 = 1.e-2_8, epsm3 = 1.e-3_8, epsm4 = 1.e-4_8
REAL, PARAMETER :: epsm5 = 1.e-5_8, epsm6 = 1.e-6_8, epsm7 = 1.e-7_8, epsm8 = 1.e-8_8
REAL, PARAMETER :: epsm10 = 1.e-10_8, epsm20 = 1.e-20_8, epsm30 = 1.e-30_8

CHARACTER(20),  PARAMETER :: AxSolverName(2) = (/'P1 SENM ',  'SP3 SENM'/)
CHARACTER(1),   PARAMETER :: DOT = '.', BANG = '!', BLANK = ' ', SLASH = '/', AST = '*', POUND = '#'
CHARACTER(512), PARAMETER :: BLANK0 = ' '
CHARACTER(126), PARAMETER :: hbar1 = &
'------------------------------------------------------------------------------------------------------------------------------'
CHARACTER(126),PARAMETER :: hbar2 = &
'=============================================================================================================================='
CHARACTER(132) MESG

END MODULE PARAM