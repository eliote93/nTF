#include <defines.h>
!--- CNJ Edit : 3D Calculation Initialization Module with Intel MKL
#ifdef __INTEL_MKL

MODULE MKL_INIT
    
USE MKL_3D
IMPLICIT NONE

LOGICAL :: lMKLInit = FALSE

CONTAINS

SUBROUTINE SetMKLEnv(CoreInfo, FmInfo, RayInfo)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,       FmInfo_Type,        RayInfo_Type
USE CORE_MOD,       ONLY : GroupInfo,           GcGroupInfo
USE CNTL,           ONLY : nTracerCntl
USE CMFD_COMMON,    ONLY : AllocHomoXSVar
USE MKL_CHEBYSHEV,  ONLY : ChebyshevInit
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(FmInfo_Type) :: FmInfo
TYPE(RayInfo_Type) :: RayInfo

mklCntl%lAxRefFDM = nTracerCntl%lAxRefFDM
mklCntl%lRefPinFDM = nTracerCntl%lNonFuelPinFDM

CALL SetGeometry(CoreInfo, FmInfo, RayInfo, GroupInfo, GcGroupInfo, nTracerCntl%l3dim)

mklCMFD%ng = mklGeom%ng
mklGcCMFD%ng = mklGeom%ngc

CALL AllocHomoXSVar(CoreInfo, mklGeom%ng)
CALL AllocCMFDVar(mklCMFD, GroupInfo, mklGeom%myzb, mklGeom%myze, mklGeom%l3dim)
IF (mklCntl%lGcCMFD) THEN
  CALL SetGcInfo(GroupInfo, GcGroupInfo)
  CALL AllocCMFDVar(mklGcCMFD, GcGroupInfo, 1, mklGeom%nzCMFD, mklGeom%l3dim)
ENDIF
IF (mklGeom%l3dim) CALL AllocAxialVar(RayInfo%PolarAngle)
IF (mklCntl%lChebyshev) CALL ChebyshevInit()

END SUBROUTINE

SUBROUTINE SetGeometry(CoreInfo, FmInfo, RayInfo, GroupInfo, GcGroupInfo, l3dim)
USE TYPEDEF,        ONLY : CoreInfo_Type,       FmInfo_Type,        RayInfo_Type,       GroupInfo_Type
USE CNTL,           ONLY : nTracerCntl
USE CMFD_COMMON,    ONLY : SetSuperPin
USE PE_MOD,         ONLY : PE
USE HexCmfdConst,   ONLY : HexCPSuperPin, HexSetGlobalPin
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(FmInfo_Type) :: FmInfo
TYPE(RayInfo_Type) :: RayInfo
TYPE(GroupInfo_Type) :: GroupInfo, GcGroupInfo
LOGICAL :: l3dim

mklGeom%ng = GroupInfo%ng
mklGeom%nFsr = CoreInfo%nCoreFsr
mklGeom%nxmax = CoreInfo%nx
mklGeom%ny = CoreInfo%ny
mklGeom%nxy = CoreInfo%nxy
mklGeom%nz = CoreInfo%nz
mklGeom%AxBC = CoreInfo%AxBC
mklGeom%myzb = PE%myzb
mklGeom%myze = PE%myze
mklGeom%nPolar2D = RayInfo%nPolarAngle
mklGeom%l3dim = l3dim

ALLOCATE(mklGeom%InScatRange(2, mklGeom%ng)); mklGeom%OutScatRange = GroupInfo%OutScatRange 
ALLOCATE(mklGeom%OutScatRange(2, mklGeom%ng)); mklGeom%InScatRange = GroupInfo%InScatRange

mklCntl%lSuperpin = mklCntl%lSuperpin .AND. CoreInfo%lGap

CALL SetCMFDPlane(CoreInfo, l3dim)

IF (nTracerCntl%lHex) THEN
  CALL HexCPSuperPin(CoreInfo, mklGeom%superPin, mklGeom%nxy, mklCntl%lSuperpin)
  CALL HexSetGlobalPin(CoreInfo, FmInfo)
ELSE
  CALL SetSuperPin(CoreInfo, mklGeom%superPin, mklGeom%nxy, mklGeom%myzb, mklGeom%myze, mklCntl%lSuperpin)
  CALL SetGlobalPin(CoreInfo, FmInfo)
END IF

END SUBROUTINE

SUBROUTINE SetCMFDPlane(CoreInfo, l3dim)
USE TYPEDEF,        ONLY : CoreInfo_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
LOGICAL :: l3dim

REAL, POINTER :: hz(:), hzfm(:)
INTEGER, POINTER :: planeMap(:)
INTEGER, POINTER :: fmRange(:, :), nzfm(:)
INTEGER :: i, ipin, ipin_map, iz, myzb, myze
INTEGER :: nz, nzCMFD
INTEGER, POINTER :: nSubplane(:)

hz => CoreInfo%hz
nz = CoreInfo%nz
myzb = mklGeom%myzb
myze = mklGeom%myze

IF (myzb .EQ. 1) mklGeom%lBottom = TRUE
IF (myze .EQ. nz) mklGeom%lTop = TRUE

ALLOCATE(nSubplane(nz))

DO iz = 1, nz
  nSubplane(iz) = INT(hz(iz) / mklCntl%CMFDHeight) + 1
ENDDO

IF (.NOT. mklCntl%lSubplane) nSubplane = 1

nzCMFD = 0
DO iz = myzb, myze
  nzCMFD = nzCMFD + nSubplane(iz)
ENDDO

ALLOCATE(mklGeom%hz(nz))
ALLOCATE(mklGeom%hzfm(0 : nzCMFD + 1))
ALLOCATE(mklGeom%nzfm(nz))
ALLOCATE(mklGeom%planeMap(nzCMFD))
ALLOCATE(mklGeom%fmRange(myzb : myze, 2))

IF (mklCntl%lGcCMFD) THEN
  ALLOCATE(mklGcCMFD%planeMap(nzCMFD))
  DO iz = 1, nzCMFD
    mklGcCMFD%planeMap(iz) = iz
  ENDDO
ENDIF

hzfm => mklGeom%hzfm
planeMap => mklGeom%planeMap
fmRange => mklGeom%fmRange
nzfm => mklGeom%nzfm
mklGeom%nzCMFD = nzCMFD
mklGeom%hz = hz

iz = 1
DO i = myzb, myze
  nzCMFD = nSubplane(i)
  iz = iz + nzCMFD
  planeMap(iz - nzCMFD : iz - 1) = i
  hzfm(iz - nzCMFD : iz - 1) = hz(i) / nzCMFD
  fmRange(i, 1) = iz - nzCMFD
  fmRange(i, 2) = iz - 1
  nzfm(i) = fmRange(i, 2) - fmRange(i, 1) + 1
ENDDO

IF (.NOT. mklGeom%lBottom) THEN
  hzfm(0) = hz(myzb - 1) / nSubplane(myzb - 1)
ELSE
  hzfm(0) = hzfm(1)
ENDIF

IF (.NOT. mklGeom%lTop) THEN
  nzCMFD = mklGeom%nzCMFD
  hzfm(nzCMFD + 1) = hz(myze + 1) / nSubplane(myze + 1)
ELSE
  nzCMFD = mklGeom%nzCMFD
  hzfm(nzCMFD + 1) = hzfm(nzCMFD)
ENDIF

mklCMFD%planeMap => mklGeom%planeMap

ALLOCATE(mklGeom%lRefPlane(myzb : myze))

DO iz = myzb, myze
  mklGeom%lRefPlane(iz) = .NOT. CoreInfo%lFuelPlane(iz)
ENDDO

DEALLOCATE(nSubplane)

END SUBROUTINE

SUBROUTINE SetGlobalPin(CoreInfo, FmInfo)
USE TYPEDEF,        ONLY : CoreInfo_Type,       FmInfo_Type,        FxrInfo_Type,       Pin_Type,                   &
                           Cell_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(FmInfo_Type) :: FmInfo

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
TYPE(superPin_Type), POINTER :: superPin(:)
INTEGER :: node(mklGeom%nxmax, mklGeom%ny)
INTEGER :: nxy, nx, ny, nzCMFD
INTEGER :: i, j, ibd, ifxr, icel, ix, iy, ixy, ixy_map, iz, izf, ipin, ipin_map, ineighpin
REAL, POINTER :: hzfm(:)

Fxr => FmInfo%Fxr
Pin => CoreInfo%Pin
Cell => CoreInfo%CellInfo
nx = mklGeom%nxmax
ny = mklGeom%ny
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
superPin => mklGeom%superPin
hzfm => mklGeom%hzfm

ALLOCATE(mklGeom%pinMap(nxy), mklGeom%pinMapRev(-1 : nxy))
ALLOCATE(mklGeom%nx(ny), mklGeom%pinRange(2, ny), mklGeom%ixRange(2, ny))

node = 0
DO ipin = 1, nxy
  ix = superPin(ipin)%ix
  iy = superPin(ipin)%iy
  node(ix, iy) = ipin
ENDDO

DO iy = 1, ny
  DO ix = 1, nx
    IF (node(ix, iy) .GT. 0) EXIT
  ENDDO
  mklGeom%ixRange(1, iy) = ix
  DO ix = nx, 1, -1
    IF (node(ix, iy) .GT. 0) EXIT
  ENDDO
  mklGeom%ixRange(2, iy) = ix
  mklGeom%nx(iy) = mklGeom%ixRange(2, iy) - mklGeom%ixRange(1, iy) + 1
ENDDO

ipin = 0
DO iy = 1, ny
  DO ix = mklGeom%ixRange(1, iy), mklGeom%ixRange(2, iy)
    ipin = ipin + 1
    mklGeom%pinMap(ipin) = node(ix, iy)
    mklGeom%pinMapRev(node(ix, iy)) = ipin
  ENDDO
  mklGeom%pinRange(1, iy) = ipin - mklGeom%nx(iy) + 1
  mklGeom%pinRange(2, iy) = ipin
ENDDO
mklGeom%pinMapRev(-1) = -1
mklGeom%pinMapRev(0) = 0

ALLOCATE(mklGeom%PinVolFm(nxy, nzCMFD))

DO izf = 1, nzCMFD
  iz = mklGeom%planeMap(izf)
  DO ipin = 1, nxy
    ipin_map = mklGeom%pinMap(ipin)
    mklGeom%PinVolFm(ipin, izf) = superPin(ipin_map)%Area * hzfm(izf)
  ENDDO
ENDDO

ALLOCATE(mklGeom%lRefPin(nxy))
ALLOCATE(mklGeom%lRefCell(nzCMFD, nxy))

DO ipin = 1, nxy
  ipin_map = mklGeom%pinMap(ipin)
  mklGeom%lRefPin(ipin) = .NOT. ANY(superPin(ipin_map)%lFuel)
!  DO ibd = 1, 4
!    ineighpin = superPin(ipin_map)%NeighIdx(ibd); IF (ineighpin .LE. 0.0) CYCLE
!    mklGeom%lRefPin(ipin) = mklGeom%lRefPin(ipin) .AND. .NOT. ANY(superPin(ineighpin)%lFuel)
!  ENDDO
  DO izf = 1, nzCMFD
    iz = mklGeom%planeMap(izf)
    mklGeom%lRefCell(izf, ipin) = .NOT. superPin(ipin_map)%lFuel(iz)
!    DO ibd = 1, 4
!      ineighpin = superPin(ipin_map)%NeighIdx(ibd); IF (ineighpin .LE. 0.0) CYCLE
!      mklGeom%lRefCell(izf, ipin) = mklGeom%lRefCell(izf, ipin) .AND. .NOT. superPin(ineighpin)%lFuel(iz)
!    ENDDO
  ENDDO
ENDDO

ALLOCATE(mklGeom%lH2OCell(nzCMFD, nxy))

DO ixy = 1, nxy
  ixy_map = mklGeom%pinMap(ixy)
  DO izf = 1, nzCMFD
    iz = mklGeom%planeMap(izf)
    mklGeom%lH2OCell(izf, ixy) = TRUE
    DO i = 1, superPin(ixy_map)%nxy
      ipin = superPin(ixy_map)%pin(i)
      icel = Pin(ipin)%Cell(iz)
      DO j = 1, Cell(icel)%nFxr
        ifxr = Pin(ipin)%FxrIdxSt + j - 1
        IF (.NOT. Fxr(ifxr, iz)%lH2O) mklGeom%lH2OCell(izf, ixy) = FALSE
      ENDDO
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE SetGcInfo(GroupInfo, GcGroupInfo)
USE TYPEDEF,        ONLY : GroupInfo_Type
IMPLICIT NONE

TYPE(GroupInfo_Type) :: GroupInfo, GcGroupInfo

INTEGER :: ng, ngc
INTEGER :: ig, igc, ib, ie, ibc, iec

ng = mklGeom%ng
ngc = mklGeom%ngc

IF (.NOT. ASSOCIATED(mklGeom%GcStruct)) THEN
  ALLOCATE(mklGeom%GcStruct(2, 2))
  mklGeom%GcStruct(1, 1) = 1
  mklGeom%GcStruct(2, 1) = GroupInfo%UpScatRange(1) - 1
  mklGeom%GcStruct(1, 2) = GroupInfo%UpScatRange(1)
  mklGeom%GcStruct(2, 2) = GroupInfo%UpScatRange(2)
ENDIF
ALLOCATE(mklGeom%GcStructInv(ng))
DO igc = 1, ngc
  mklGeom%GcStructInv(mklGeom%GcStruct(1, igc) : mklGeom%GcStruct(2, igc)) = igc
ENDDO
ALLOCATE(GcGroupInfo%InScatRange(2, ngc))
GcGroupInfo%InScatRange(1, :) = ngc
GcGroupInfo%InScatRange(2, :) = 1
GcGroupInfo%ng = ngc
DO ig = 1, ng
  igc = mklGeom%GcStructInv(ig)
  ib = GroupInfo%InScatRange(1, ig); ie = GroupInfo%InScatRange(2, ig);  
  ibc = mklGeom%GcStructInv(ib); iec = mklGeom%GcStructInv(ie);
  GcGroupInfo%InScatRange(1, igc) = MIN(GcGroupInfo%InScatRange(1, igc), ibc)
  GcGroupInfo%InScatRange(2, igc) = MAX(GcGroupInfo%InScatRange(2, igc), iec)
ENDDO
GcGroupInfo%lUpScat = FALSE
DO igc = 1, ngc
  IF (GcGroupInfo%InScatRange(2, igc) .GT. igc) THEN
    GcGroupInfo%lUpScat = TRUE
    EXIT
  ENDIF
ENDDO
IF (GcGroupInfo%lUpScat) THEN
  GcGroupInfo%UpScatRange(1) = ngc; GcGroupInfo%UpScatRange(2) = 1
  DO igc = 1, ngc
    IF (GcGroupInfo%InScatRange(2, igc) .GT. igc) THEN
      GcGroupInfo%UpScatRange(1) = MIN(GcGroupInfo%UpScatRange(1), igc)
      GcGroupInfo%UpScatRange(2) = MAX(GcGroupInfo%UpScatRange(2), igc)
    ENDIF
  ENDDO
  GcGroupInfo%UpScatRange(2) = ngc
ENDIF

END SUBROUTINE

SUBROUTINE AllocCMFDVar(CMFD, GroupInfo, myzb, myze, l3dim)
USE TYPEDEF,        ONLY : GroupInfo_Type
USE CMFD_COMMON,    ONLY : AllocPinXS
IMPLICIT NONE

TYPE(mklCMFD_Type) :: CMFD
TYPE(GroupInfo_Type) :: GroupInfo
INTEGER :: myzb, myze
LOGICAL :: l3dim

INTEGER :: ig, ng, nxy, nzCMFD

ng = CMFD%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

CMFD%bottomRange = (/ 1, nxy /)
CMFD%topRange = (/ nxy * (nzCMFD - 1) + 1, nxy * nzCMFD /)
ALLOCATE(CMFD%InScatRange(2, ng)); CMFD%InScatRange = GroupInfo%InScatRange

ALLOCATE(CMFD%M(ng))
ALLOCATE(CMFD%S(nxy * nzCMFD, ng, ng))
ALLOCATE(CMFD%F(nxy * nzCMFD, ng))
ALLOCATE(CMFD%Chi(nxy * nzCMFD, ng))
  
IF (mklCntl%lJacobi) THEN
  ALLOCATE(CMFD%Jacobi(ng))
  DO ig = 1, ng
    ALLOCATE(CMFD%Jacobi(ig)%invDiag(nxy * nzCMFD))
    IF (l3dim) ALLOCATE(CMFD%Jacobi(ig)%offDiag(nxy, 2))
  ENDDO
ENDIF

IF (l3dim) THEN
  IF (mklCntl%lDavidson) THEN
    ALLOCATE(CMFD%AxOffDiag(ng, nxy, 2))
    ALLOCATE(CMFD%trAxOffDiag(ng, nxy, 2))
    CMFD%AxOffDiag = 0.0
    CMFD%trAxOffDiag = 0.0
  ELSE
    ALLOCATE(CMFD%AxOffDiag(nxy, 2, ng))
    CMFD%AxOffDiag = 0.0
  ENDIF
  IF (mklCntl%DcplLv .EQ. 2) THEN
    ALLOCATE(CMFD%dcplAxOffDiag(nxy, nzCMFD, 2, ng)); CMFD%dcplAxOffDiag = 0.0
  ENDIF
ENDIF
  
ALLOCATE(CMFD%phis(nxy, nzCMFD, ng))
ALLOCATE(CMFD%phic(nxy, myzb : myze, ng))
ALLOCATE(CMFD%neighphis(nxy, ng, 2))
ALLOCATE(CMFD%src(nxy * nzCMFD, ng))
ALLOCATE(CMFD%psi(nxy, nzCMFD))
ALLOCATE(CMFD%psid(nxy, nzCMFD))

IF (mklCntl%lChebyshev) ALLOCATE(CMFD%phisd(nxy, nzCMFD, ng, 2))

IF (mklCntl%lSPAI) THEN
  ALLOCATE(CMFD%SPAI(ng))
ELSE
  ALLOCATE(CMFD%ILU(ng))
ENDIF

IF (mklCntl%lDavidson) THEN
  ALLOCATE(CMFD%Davidson)
  ALLOCATE(CMFD%Davidson%u(ng * nxy * nzCMFD))
  ALLOCATE(CMFD%Davidson%t(ng * nxy * nzCMFD))
  ALLOCATE(CMFD%Davidson%r(ng * nxy * nzCMFD))
ENDIF

ALLOCATE(CMFD%AxDtil(2, nxy, nzCMFD, ng))
ALLOCATE(CMFD%AxDhat(2, nxy, nzCMFD, ng))
CMFD%AxDhat = 0.0

ALLOCATE(CMFD%theta(ng, nxy, nzCMFD)); CMFD%theta = 0.0

CALL AllocPinXS(CMFD%PinXS, GroupInfo, nxy, myzb, myze)

END SUBROUTINE

SUBROUTINE AllocAxialVar(PolarAngle)
USE TYPEDEF,        ONLY : PolarAngle_Type
USE MKL_LINMOC,     ONLY : AllocLinearMOC
USE MKL_FLATMOC,    ONLY : AllocFlatMOC
USE MKL_NODAL,      ONLY : AllocNodal
USE MKL_FDM,        ONLY : AllocFDM
IMPLICIT NONE

TYPE(PolarAngle_Type), POINTER :: PolarAngle(:)
INTEGER :: ng, nFsr, nxy, nzCMFD, nPolar1D, ScatOd
INTEGER :: myzb, myze

ng = mklGeom%ng
nFsr = mklGeom%nFsr
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
nPolar1D = mklGeom%nPolar1D
myzb = mklGeom%myzb
myze = mklGeom%myze
ScatOd = mklCntl%scatOrder

ALLOCATE(mklAxial%phic(ng, nxy, nzCMFD)); mklAxial%phic = 0.0
ALLOCATE(mklAxial%Jout(3, ng, 2, nzCMFD, nxy)); mklAxial%Jout = 0.0

IF (mklCntl%lMOC) THEN

  ALLOCATE(mklAxial%Angle(nPolar1D))

  IF (nPolar1D .GT. 4) THEN
    CALL GaussLegendre()
  ELSE
    CALL TabuchiYamamoto()
  ENDIF
  
ENDIF

SELECT CASE (mklCntl%AxSolver)
CASE (NODAL)
  CALL AllocNodal()
CASE (FDM)
  CALL AllocFDM()
CASE (MOC)
  ALLOCATE(mklAxial%PhiAngIn(nPolar1D, ng, nxy, 2))
  ALLOCATE(mklAxial%PhiAngOut(nPolar1D, ng, nxy, 2))
  ALLOCATE(mklAxial%SmP1(ng, ng, nzCMFD, nxy)); mklAxial%SmP1 = 0.0
  ALLOCATE(mklAxial%SmP2(ng, ng, nzCMFD, nxy)); mklAxial%SmP2 = 0.0
  ALLOCATE(mklAxial%SmP3(ng, ng, nzCMFD, nxy)); mklAxial%SmP3 = 0.0
  IF (mklCntl%lCASMO) THEN
    CALL AllocLinearMOC()
  ELSE
    CALL AllocFlatMOC()
  ENDIF
END SELECT

ALLOCATE(mklAxial%atil(2, nxy, nzCMFD, ng))
   
END SUBROUTINE

SUBROUTINE GaussLegendre()

IMPLICIT NONE

INTEGER :: nPolar1D
INTEGER :: ipol
REAL :: mu
REAL, POINTER :: abscissa(:), weight(:)

nPolar1D = mklGeom%nPolar1D

ALLOCATE(abscissa(nPolar1D * 2), weight(nPolar1D * 2))

CALL gauleg(nPolar1D * 2, abscissa, weight)

DO ipol = 1, nPolar1D
  mu = - abscissa(ipol)
  mklAxial%Angle(ipol)%cosv = mu
  mklAxial%Angle(ipol)%rcosv = 1.0 / mu
  mklAxial%Angle(ipol)%wt = weight(ipol) / 2.0
  mklAxial%Angle(ipol)%wtsurf = weight(ipol) * mu / 2.0
ENDDO

DEALLOCATE(abscissa, weight)
  
CONTAINS

!********************************************************************************
!* Calculation of GAUSS-LEGENDRE abscissas and weights for Gaussian Quadrature
!* integration of polynomial functions.
!*      For normalized lower and upper limits of integration -1.0 & 1.0, and
!* given n, this routine calculates, arrays xabsc(1:n) and  weig(1:n) of length n,
!* containing the abscissas and weights of the Gauss-Legendre n-point quadrature
!* formula.  For detailed explanations finding weights & abscissas, see
!* "Numerical Recipes in Fortran */
!********************************************************************************
	SUBROUTINE  gauleg(ngp, xabsc, weig)

      implicit none
      INTEGER, parameter :: dbp = SELECTED_REAL_KIND (15,307)
      INTEGER  i, j, m
      REAL(dbp)  p1, p2, p3, pp, z, z1
      INTEGER, INTENT(IN) :: ngp            ! # of Gauss Points
      REAL(dbp), INTENT(OUT) :: xabsc(ngp), weig(ngp)
      REAL(dbp)  :: EPS, M_PI
      PARAMETER (EPS=3.0d-15)       	!EPS is the relative precision
      PARAMETER (M_PI=3.141592654d0)      ! Pi value

	   m = (ngp + 1) / 2
!* Roots are symmetric in the interval - so only need to find half of them  */

	   do i = 1, m				! Loop over the desired roots */

     		z = cos( M_PI * (i-0.25d0) / (ngp+0.5d0) )
!*   Starting with the above approximation to the ith root,
!*          we enter the main loop of refinement by NEWTON'S method   */
100     	p1 = 1.0d0
        	p2 = 0.0d0
!*  Loop up the recurrence relation to get the Legendre
!*  polynomial evaluated at z                 */

        	do j = 1, ngp
           	p3 = p2
           	p2 = p1
           	p1 = ((2.0d0*j-1.0d0) * z * p2 - (j-1.0d0)*p3) / j
        	enddo

!* p1 is now the desired Legendre polynomial. We next compute pp,
!* its derivative, by a standard relation involving also p2, the
!* polynomial of one lower order.      */
        	pp = ngp*(z*p1-p2)/(z*z-1.0d0)
        	z1 = z
        	z = z1 - p1/pp             ! Newton's Method  */

        	if (dabs(z-z1) .gt. EPS) GOTO  100

      	xabsc(i) =  - z                    	! Roots will be bewteen -1.0 & 1.0 */
      	xabsc(ngp+1-i) =  + z                	! and symmetric about the origin  */
      	weig(i) = 2.0d0/((1.0d0-z*z)*pp*pp) ! Compute the weight and its       */
      	weig(ngp+1-i) = weig(i)               ! symmetric counterpart         */

      end do     ! i loop

   End subroutine gauleg
   
END SUBROUTINE

SUBROUTINE TabuchiYamamoto()

IMPLICIT NONE

INTEGER :: nPolar1D
INTEGER :: ipol
REAL :: mu
REAL :: weight(4, 4), sinv(4, 4)

DATA weight / 1.0000000, 0.0000000, 0.0000000, 0.0000000, &
              0.2128540, 0.7871460, 0.0000000, 0.0000000, &
              0.0462330, 0.2836190, 0.6701480, 0.0000000, &
              0.0116656, 0.0855345, 0.3079014, 0.5948985  /
DATA sinv   / 0.7071068, 0.0000000, 0.0000000, 0.0000000, &
              0.3639000, 0.8999000, 0.0000000, 0.0000000, &
              0.1666480, 0.5377070, 0.9329540, 0.0000000, &
              0.0834435, 0.2908650, 0.6328485, 0.9492286  /

nPolar1D = mklGeom%nPolar1D

DO ipol = 1, nPolar1D
  mu = SQRT(1.0 - sinv(ipol, nPolar1D) ** 2)
  mklAxial%Angle(ipol)%cosv = mu
  mklAxial%Angle(ipol)%rcosv = 1.0 / mu
  mklAxial%Angle(ipol)%wt = weight(ipol, nPolar1D) / 2.0
  mklAxial%Angle(ipol)%wtsurf = weight(ipol, nPolar1D) * mu / 2.0
ENDDO

END SUBROUTINE

SUBROUTINE SetPolarMap(PolarAngle, Angle)
USE TYPEDEF,        ONLY : PolarAngle_Type
IMPLICIT NONE

TYPE(PolarAngle_Type), POINTER :: PolarAngle(:)
TYPE(mklAngle_Type) :: Angle

INTEGER :: nPolar2D
INTEGER :: ipol
REAL :: cosv

nPolar2D = mklGeom%nPolar2D
cosv = Angle%cosv

IF (PolarAngle(1)%cosv .LT. cosv) THEN
  Angle%wt2D(2) = (cosv - PolarAngle(1)%cosv) / (PolarAngle(2)%cosv - PolarAngle(1)%cosv)
  Angle%wt2D(1) = 1.0 - Angle%wt2D(2)
  Angle%polar2D(1) = 1
  Angle%polar2D(2) = 2
  RETURN
ENDIF

DO ipol = 1, nPolar2D - 1
  IF (PolarAngle(ipol)%cosv .GT. cosv .AND. PolarAngle(ipol + 1)%cosv .LT. cosv) THEN
    Angle%wt2D(2) = (cosv - PolarAngle(ipol)%cosv) / (PolarAngle(ipol + 1)%cosv - PolarAngle(ipol)%cosv)
    Angle%wt2D(1) = 1.0 - Angle%wt2D(2)
    Angle%polar2D(1) = ipol
    Angle%polar2D(2) = ipol + 1
    RETURN
  ENDIF
ENDDO

Angle%wt2D(2) = (cosv - PolarAngle(nPolar2D - 1)%cosv) / (PolarAngle(nPolar2D)%cosv - PolarAngle(nPolar2D - 1)%cosv)
Angle%wt2D(1) = 1.0 - Angle%wt2D(2)
Angle%polar2D(1) = nPolar2D - 1
Angle%polar2D(2) = nPolar2D
RETURN

END SUBROUTINE

! SUBROUTINE AllocBILUVar(CMFD)
! 
! IMPLICIT NONE
! 
! TYPE(mklCMFD_Type) :: CMFD
! 
! TYPE(mklBILU_Type), POINTER :: BILU
! TYPE(blockMat_Type), POINTER :: blockMat
! TYPE(triLU_Type), POINTER :: blockDiag, Del
! INTEGER :: ng, nx, nxmax, ny, nxy, nzCMFD
! INTEGER :: ig, iz, iy
! 
! ng = CMFD%ng
! nxmax = mklGeom%nxmax
! ny = mklGeom%ny
! nxy = mklGeom%nxy
! nzCMFD = mklGeom%nzCMFD
! 
! ALLOCATE(CMFD%BILU(ng))
! 
! DO ig = 1, ng
!   BILU => CMFD%BILU(ig)
!   ALLOCATE(BILU%lower(nxy, nzCMFD), BILU%upper(nxy, nzCMFD))
!   ALLOCATE(BILU%blockMat(nzCMFD))
!   DO iz = 1, nzCMFD
!     blockMat => BILU%blockMat(iz)
!     blockMat%ny = ny
!     ALLOCATE(blockMat%lower(nxmax, ny), blockMat%upper(nxmax, ny))
!     ALLOCATE(blockMat%blockDiag(ny), blockMat%Del(ny))
!     DO iy = 1, ny
!       blockDiag => blockMat%blockDiag(iy); Del => blockMat%Del(iy)
!       nx = mklGeom%nx(iy); blockDiag%nx = nx; Del%nx = nx
!       ALLOCATE(blockDiag%lower(2 : nx), blockDiag%diag(1 : nx), blockDiag%upper(1 : nx - 1))
!       ALLOCATE(Del%lower(2 : nx), Del%diag(1 : nx), Del%upper(1 : nx - 1))
!     ENDDO
!   ENDDO
! ENDDO
! 
! END SUBROUTINE

END MODULE
    
#endif