#include <defines.h>
MODULE Tlkg_mod
USE PARAM
USE TYPEDEF,   ONLY : AxFlx_TYPE,            PinXS_TYPE,            &
                      CoreInfo_TYPE,         Pin_TYPE,    CmInfo_Type
IMPLICIT NONE

INTEGER, PARAMETER :: ConvLkg2ndMod = 0
INTEGER, PARAMETER :: IntraLkg2ndMod = 1

INTERFACE RadTlkgUpdate
  MODULE PROCEDURE RadTlkgUpdate_type0
  MODULE PROCEDURE RadTlkgUpdate_type1
  MODULE PROCEDURE RadTlkgUpdate_type2
END INTERFACE

CONTAINS


SUBROUTINE RadTlkgUpdate_type0(Core, CmInfo, TLKG, ixy1, ixy2, iz1, iz2, ng)
IMPLICIT NONE
TYPE(CoreInfo_TYPE) :: Core
TYPE(CmInfo_TYPE) :: CmInfo
REAL, POINTER :: Tlkg(:, :, :)
INTEGER :: ixy1, ixy2, ng, iz1, iz2


TYPE(PinXS_TYPE), POINTER :: PinXS(:, :)
REAL, POINTER :: PhiC(:, :, :)

REAL, POINTER :: PinVol(:, :), Hz(:)
TYPE(Pin_TYPE), POINTER :: Pin(:)
INTEGER :: iz, iz0, ixy, ig, i, j, k
INTEGER :: nxy, nbd
REAL :: RadLkg, Jnet, Dtil, Dhat, neighphi, MyPhi

PinXS => CmInfo%PinXS
PhiC => CmInfo%PhiFm

Pin => Core%Pin
Hz => Core%Hzfm
PinVol => Core%PinVolFm
nxy = Core%nxy
nbd = 4
DO ig = 1, ng
  DO iz = iz1, iz2
    iz0  =Core%SubPlaneMap(iz)
    DO ixy = 1, nxy
      RadLkg = 0
      MyPhi = PhiC(ixy, iz, ig)
      DO i = 1, nbd
        DHat = PinXS(ixy, iz0)%DHat(i, ig)
        DTil = PinXS(ixy, iz0)%Dtil(i, ig)
        j = Pin(ixy)%NeighIdx(i)
        neighphi = 0.
        IF(j .GT. 0) THEN
          neighphi = PhiC(j, iz, ig)
        ENDIF
        Jnet = (Dtil - Dhat) * MyPhi - (Dtil + Dhat) * NeighPhi
        RadLkg = RadLkg + Jnet
      ENDDO
      RadLkg = RadLkg * Hz(iz) / PinVol(ixy, iz0)
      TLKG(ixy, iz, ig) = RadLkg
    ENDDO
  ENDDO
ENDDO
NULLIFY(PhiC, PinXS)

NULLIFY(Pin)   
NULLIFY(Hz); NULLIFY(PinVol)
END SUBROUTINE

SUBROUTINE RadTlkgUpdate_type1(Core, PhiC, PinXS, TLKG, ixy1, ixy2, iz1, iz2, ng)
IMPLICIT NONE
TYPE(CoreInfo_TYPE) :: Core
TYPE(PinXS_TYPE) :: PinXS(ixy1:ixy2, iz1:iz2)
REAL :: Tlkg(ixy1:ixy2, iz1:iz2, ng), PhiC(ixy1:ixy2, iz1:iz2, ng)
INTEGER :: ixy1, ixy2, ng, iz1, iz2

REAL, POINTER :: PinVol(:, :), Hz(:)
TYPE(Pin_TYPE), POINTER :: Pin(:)
INTEGER :: iz, iz0, ixy, ig, i, j, k
INTEGER :: nxy, nbd
REAL :: RadLkg, Jnet, Dtil, Dhat, neighphi, MyPhi


Pin => Core%Pin
Hz => Core%Hzfm
PinVol => Core%PinVolFm
nxy = Core%nxy
nbd = 4
DO ig = 1, ng
  DO iz = iz1, iz2
    iz0  =Core%SubPlaneMap(iz)
    DO ixy = 1, nxy
      RadLkg = 0
      MyPhi = PhiC(ixy, iz, ig)
      DO i = 1, nbd
        DHat = PinXS(ixy, iz0)%DHat(i, ig)
        DTil = PinXS(ixy, iz0)%Dtil(i, ig)
        j = Pin(ixy)%NeighIdx(i)
        neighphi = 0.
        IF(j .GT. 0) THEN
          neighphi = PhiC(j, iz, ig)
        ENDIF
        Jnet = (Dtil - Dhat) * MyPhi - (Dtil + Dhat) * NeighPhi
        RadLkg = RadLkg + Jnet
      ENDDO
      RadLkg = RadLkg * Hz(iz) / PinVol(ixy, iz0)
      TLKG(ixy, iz, ig) = RadLkg
    ENDDO
  ENDDO
ENDDO

NULLIFY(Pin)   
NULLIFY(Hz); NULLIFY(PinVol)
END SUBROUTINE

SUBROUTINE RadTlkgUpdate_type2(Core, PinXS, AxFlx, ixy1, ixy2, iz1, iz2, ng)
IMPLICIT NONE
TYPE(CoreInfo_TYPE) :: Core
INTEGER :: ixy1, ixy2, ng, iz1, iz2
TYPE(PinXS_Type), POINTER :: PinXS(:, :)
TYPE(AxFlx_Type), POINTER :: AxFlx(:, :)
!REAL :: PHIC(ixy1 : ixy2, iz1 : iz2, ng)


REAL, POINTER :: PinVol(:, :), Hz(:)
TYPE(Pin_TYPE), POINTER :: Pin(:)
INTEGER :: iz, iz0, ixy, ig, i, j, k
INTEGER :: nxy, nbd
REAL :: RadLkg, Jnet, Dtil, Dhat, neighphi, MyPhi


Pin => Core%Pin
Hz => Core%Hzfm
PinVol => Core%PinVolFm
nxy = Core%nxy
nbd = 4
DO ig = 1, ng
  DO iz = iz1, iz2
    iz0  =Core%SubPlaneMap(iz)
    DO ixy = 1, nxy
      RadLkg = 0
      !MyPhi = PhiC(ixy, iz, ig)
      MyPhi = AxFlx(iz, ixy)%Phi(0, 1, ig)
      DO i = 1, nbd
        DHat = PinXS(ixy, iz0)%DHat(i, ig)
        DTil = PinXS(ixy, iz0)%Dtil(i, ig)
        j = Pin(ixy)%NeighIdx(i)
        neighphi = 0.
        IF(j .GT. 0) THEN
          !neighphi = PhiC(j, iz, ig)
          neighphi = AxFlx(iz, j)%Phi(0, 1, ig)
        ENDIF
        Jnet = (Dtil - Dhat) * MyPhi - (Dtil + Dhat) * NeighPhi
        RadLkg = RadLkg + Jnet
      ENDDO
      RadLkg = RadLkg * Hz(iz) / PinVol(ixy, iz0)
      AxFlx(iz, ixy)%Tlkg(0, 1, ig) = RadLkg
    ENDDO
  ENDDO
ENDDO
!DO iz = iz1, iz2
!  WRITE(*, '(100e15.5)') (AxFlx(i, 1)%Tlkg(0, 1, 1), i = iz1, iz2)
!ENDDO
NULLIFY(Pin)   
NULLIFY(Hz); NULLIFY(PinVol)
END SUBROUTINE

SUBROUTINE IntraLkgShape(Core, AxFlx, PinXS, ixy1, ixy2, iz1, iz2, ng)
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(AxFlx_Type) :: AxFlx(iz1:iz2, ixy1:ixy2)
TYPE(PinXS_TYPE) :: PinXS(ixy1:ixy2, iz1:iz2)
INTEGER :: ixy1, ixy2, iz1, iz2, ng

TYPE(pin_type), POINTER :: Pin(:)
REAL, POINTER :: PinVol(:, :), Hz(:)

INTEGER :: iz, iz0, ixy, ig, ibd, nbd, ineigh
REAL :: Dhat, Dtil
REAL :: LkgShape(1:4), myshape(1:4), neighshape(1:4), jshape(1:4)
REAL :: SubNodeLkg(-1:1), MySubPhi(-1:1), NeighSubPhi(-1:1), SubNodeJ(-1:1)
Pin => Core%Pin
Hz => Core%Hzfm
PinVol => Core%PinVolFm

nbd = 4
DO iz = iz1, iz2
  iz0  =Core%SubPlaneMap(iz)
  DO ixy = ixy1, ixy2
    DO ig = 1, ng
      LkgShape = 0
      DO ibd = 1, nbd
        Dhat = PinXS(ixy, iz0)%Dhat(ibd, ig)
        Dtil = PinXS(ixy, iz0)%Dtil(ibd, ig)
        myshape = AxFlx(iz, ixy)%phi(1:4, 1, ig)
        neighShape = 0
        ineigh = Pin(ixy)%NeighIdx(ibd)
        IF(ineigh .GT. 0) THEN
          NeighShape = AxFlx(iz, ineigh)%phi(1:4, 1, ig)
        ENDIF
        Jshape = (Dtil - Dhat) * Myshape- (Dtil + Dhat) * NeighShape
        LkgShape = LkgShape + JShape
      ENDDO
      LkgShape = LkgShape * Hz(iz) / PinVol(ixy, iz0)
      AxFlx(iz, ixy)%LkgShape(1:4, ig) = LkgShape
    ENDDO
  ENDDO
ENDDO
#ifdef ExactSubLkg
DO iz = iz1, iz2
  iz0  =Core%SubPlaneMap(iz)
  DO ixy = ixy1, ixy2
    DO ig = 1, ng
      SubNodeLkg = 0
      DO ibd = 1, nbd
        Dhat = PinXS(ixy, iz0)%Dhat(ibd, ig)
        Dtil = PinXS(ixy, iz0)%Dtil(ibd, ig)
        MySubPhi = AxFlx(iz, ixy)%SubNodePhi(-1:1, ig)
        ineigh = Pin(ixy)%NeighIdx(ibd)
        IF(ineigh .GT. 0) THEN
          NeighSubPhi = AxFlx(iz, ineigh)%SubNodePhi(-1:1, ig)
        ENDIF
        SubNodeJ = (Dtil - Dhat) * MySubPhi- (Dtil + Dhat) * NeighSubPhi
        SubNodeLkg = SubNodeLkg + SubNodeJ
      ENDDO
      SubNodeLkg = SubNodeLkg * Hz(iz) / PinVol(ixy, iz0)
      AxFlx(iz, ixy)%SubNodeLkg(-1:1, ig) = SubNodeLkg
    ENDDO
  ENDDO
ENDDO
#endif
NULLIFY(Pin, Hz, PinVol)
END SUBROUTINE

SUBROUTINE SetTLkgShape0(AxFlx, bc, nz, ng, lsolver, lkgmod)
USE PARAM
USE TYPEDEF, ONLY : AxFlx_Type
IMPLICIT NONE
TYPE(AxFlx_Type) :: AxFlx(nz)
INTEGER ::  bc(2), lsolver, nz, ng
REAL :: tlkg(0 : 2, 0 : nz + 1), avgtlkg(-1 : 1)
REAL :: coeff(0:2,2)
INTEGER :: lkgmod

INTEGER :: iz, ig
!lkgmod =0
DO ig = 1, ng
  tlkg(0, 0) = 0; tlkg(0, nz + 1) = 0;
  DO iz = 1, nz
    tlkg(0, iz) = AxFlx(iz)%Tlkg(0, 1, ig)
  ENDDO
  if(BC(1) .EQ. RefCell) tlkg(0, 0) = tlkg(0, 1)
  if(BC(2) .EQ. RefCell) tlkg(0, nz + 1) = tlkg(0, nz)
  DO iz = 1, nz
    avgtlkg = tlkg(0, iz - 1 : iz +1)
    tlkg(0:2, iz) = Convtlkg2nd(avgtlkg)
  ENDDO
  Do iz = 1, nz
    AxFlx(iz)%Tlkg(:, 1, ig) = tlkg(0:2, iz) 
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE SetTLkgShape(AxFlx, bc, hz, nz, ng, lsolver, lkgmod)
USE PARAM
USE TYPEDEF, ONLY : AxFlx_Type
IMPLICIT NONE
TYPE(AxFlx_Type) :: AxFlx(nz)
REAL :: hz(nz)
INTEGER ::  bc(2), lsolver, nz, ng
REAL :: tlkg(0 : 2, 0 : nz + 1), avgtlkg(-1 : 1), hz0(0:nz+1), hr(-1:1)
REAL :: coeff(0:2,2)
INTEGER :: lkgmod

INTEGER :: iz, ig
!lkgmod =0
DO iz = 1, nz
  hz0(iz) = hz(iz)  
ENDDO
hz0(0) = hz0(1); hz0(nz+1) = hz0(nz)
DO ig = 1, ng
  tlkg(0, 0) = 0; tlkg(0, nz + 1) = 0;
  DO iz = 1, nz
    tlkg(0, iz) = AxFlx(iz)%Tlkg(0, 1, ig) !(iz)%(avg, ?, g)
  ENDDO
  if(BC(1) .EQ. RefCell) tlkg(0, 0) = tlkg(0, 1)
  if(BC(2) .EQ. RefCell) tlkg(0, nz + 1) = tlkg(0, nz)
  DO iz = 1, nz
    avgtlkg = tlkg(0, iz - 1 : iz +1)  !original
    !avgtlkg(-1) = tlkg(0, iz)
    !avgtlkg(0) = tlkg(0, iz )
    !avgtlkg(+1) = tlkg(0, iz)
    hr(0) = 1; hr(-1) = hz0(iz-1) / hz0(iz);hr(1) = hz0(iz+1) / hz0(iz)
    tlkg(0:2, iz) = Convtlkg2nd(avgtlkg)    ! ORIGINAL
    !tlkg(0:2, iz) = Convtlkg2nd0(avgtlkg,hr) ! Hr
  ENDDO
  Do iz = 1, nz
    AxFlx(iz)%Tlkg(:, 1, ig) = tlkg(0:2, iz) 
  ENDDO
ENDDO

END SUBROUTINE



FUNCTION IntraLkg2nd(avgtlkg, shape4th)
REAL :: IntraLkg2nd(0:2)
REAL :: avgtlkg, shape4th(1:4)

IntraLkg2nd(0) = avgtlkg
IntraLkg2nd(1) = Shape4th(1) - Shape4th(3) / 9._8
IntraLkg2nd(2) = Shape4th(2) - 5._8 * Shape4th(4) / 9._8

END FUNCTION

FUNCTION IntraLkg2nd_2(avgtlkg, subnodelkg)
REAL :: IntraLkg2nd_2(0:2)
REAL :: avgtlkg, subnodelkg(-1:1)

IntraLkg2nd_2(0) = avgtlkg
IntraLkg2nd_2(1) = -0.75_8*(subnodelkg(-1)-subnodelkg(1))
IntraLkg2nd_2(2) = -0.75_8*(2._8*subnodelkg(0)-subnodelkg(-1)-subnodelkg(1))
!IntraLkg2nd(1) = Shape4th(1) - Shape4th(3) / 9._8
!IntraLkg2nd(2) = Shape4th(2) - 5._8 * Shape4th(4) / 9._8

END FUNCTION

FUNCTION Convtlkg2nd(avgtlkg)
REAL :: avgtlkg(-1:1),Convtlkg2nd(0:2)
  Convtlkg2nd(0) = avgtlkg(0) 
  Convtlkg2nd(1) = 0.25_8 * (-avgtlkg(-1) + avgtlkg(1))
  Convtlkg2nd(2) = 0.0833333333333333333_8 * (-2._8 * avgtlkg(0) + avgtlkg(-1) + avgtlkg(1))
  !tlkg2nd_coeff(0:0) = 0
END FUNCTION

FUNCTION Convtlkg2nd0(avgtlkg, hr)
REAL :: avgtlkg(-1:1),Convtlkg2nd0(0:2), hr(-1:1)

  Convtlkg2nd0(0) = avgtlkg(0) 
  Convtlkg2nd0(1) = -(avgtlkg(-1)-avgtlkg(0))/((1._8+hr(-1))*(1._8+2*hr(-1)))
  Convtlkg2nd0(1) = Convtlkg2nd0(1) + (avgtlkg(1)-avgtlkg(0))/((1._8+hr(1))*(1._8+2*hr(1)))
  Convtlkg2nd0(1) = Convtlkg2nd0(1) /(1._8/(1._8+2*hr(-1))+1._8/(1._8+2*hr(1)))
  Convtlkg2nd0(2) = (avgtlkg(-1)-avgtlkg(0)) / (1._8+hr(-1)) + (avgtlkg(1)-avgtlkg(0)) / (1._8+hr(1))
  Convtlkg2nd0(2) = Convtlkg2nd0(2) / (2._8 + 2._8*hr(-1) + 2._8*hr(1))
  !Convtlkg2nd(1) = 0.25_8 * (-avgtlkg(-1) + avgtlkg(1))
  !Convtlkg2nd(2) = 0.0833333333333333333_8 * (-2._8 * avgtlkg(0) + avgtlkg(-1) + avgtlkg(1))
  !tlkg2nd_coeff(0:0) = 0
END FUNCTION

SUBROUTINE PrintLkg(Core,AxFlx, nz, ng)
TYPE(CoreInfo_Type) :: Core
TYPE(AxFlx_TYPE) :: AxFlx(nz)
INTEGER :: nz, ng
INTEGER :: ig, iz
INTEGER :: i
INTEGER, SAVE :: iter
DATA iter /0/
REAL, POINTER :: hz(:)
REAL :: x, xh, hsum, lkg(0:10), coeff(0:4)
character(20) :: filename
iter = iter + 1
IF(iter < 10) THEN
  write(filename,'(A7,A2,I1,A4)') 'Lkginfo','00',iter,'.out'
ELSEIF(iter <100) THEN
  write(filename,'(A7,A1,I2,A4)') 'Lkginfo','0',iter,'.out'
ELSE
  write(filename,'(A7,I3,A4)') 'Lkginfo', iter,'.out'
ENDIF
hz => Core%hz
OPEN(unit=99, file = filename, STATUS = 'REPLACE')
ig = 7
hsum =0 
DO iz = 1, nz
  x=-1.0_8 - 0.05 
  DO i = 1, 41
    x = x + 0.05
    xh = (x + 1._8) * 0.5 * hz(iz) + hsum
    coeff(0:2) = AxFlx(iz)%tlkg(0:2, 1, ig)
    Lkg(0) = AxFlx(iz)%tlkg(0, 1, ig)
    Lkg(1) = LpExpFtn(coeff(0:2), x, 2)
    coeff(1:4) = AxFlx(iz)%LkgShape(1:4, ig)
    Lkg(2) = LpExpFtn(coeff(0:4), x, 4)
    coeff(0:2) = IntraLkg2nd(AxFlx(iz)%tlkg(0, 1, ig), AxFlx(iz)%LkgShape(1:4, ig))
    Lkg(3) = LpExpFtn(coeff(0:2), x, 2)
    Lkg(4) = LpExpFtn(AxFlx(iz)%phi(0:4, 1, ig), x, 2)
    WRITE(99, '(2i5, f10.5, 10E15.5)') IG, IZ, Xh, lKG(0:4) 
  ENDDO
  hsum = hsum + hz(iz)
ENDDO
CLOSE(99)
NULLIFY(hz)
END SUBROUTINE

FUNCTION LpExpFtn(coeff, x, n)
REAL :: Coeff(0:n)
REAL :: x
REAL :: LpExpFtn
INTEGER :: n
INTEGER :: i
LpExpFtn = 0
DO i = 0, n
  LpExpFtn = LpExpFtn + LP(I, x) * COEFF(I) 
ENDDO

END FUNCTION

FUNCTION LP(I,X)
IMPLICIT NONE
INTEGER :: I,J
REAL :: X,Y,LP,LPCOEFF(0:4,0:4)
DATA LPCOEFF    /    1._8,    0._8,    0._8,   0._8,    0._8, &
                     0._8,    1._8,    0._8,   0._8,    0._8, &
                   -0.5_8,    0._8,   1.5_8,   0._8,    0._8, &
                     0._8,  -1.5_8,    0._8,  2.5_8,    0._8, &
                  0.375_8,    0._8, -3.75_8,   0._8, 4.375_8 /    
Y=0
DO J=0,I
  Y=Y+(X**J)*LPCOEFF(J,I)
ENDDO          
LP=Y                     
END FUNCTION


END MODULE