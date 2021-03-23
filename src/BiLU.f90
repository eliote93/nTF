#include <defines.h>
MODULE BiLU_MOD
USE PARAM
USE TYPEDEF,    ONLY : CoreInfo_TYPE,      CMFDLS_TYPE,      Pin_TYPE,      &
                       BiLU_TYPE,          PE_TYPE
USE ALLOCS
IMPLICIT NONE
INTEGER, PRIVATE :: nx, ny, nxy
INTEGER, PRIVATE :: nz, myzb, myze, nzfm, myzbf, myzef
INTEGER, POINTER, PRIVATE :: nxbeg(:), nxend(:)
INTEGER, POINTER, PRIVATE :: node(:, :)
LOGICAL, PRIVATE :: lSetGeom = .FALSE.

REAL, POINTER, PRIVATE :: Au(:), del(:)                  !Size(nx)
REAL, POINTER, PRIVATE :: ainvl(:), ainvu(:), ainvd(:)            !Size(nx)
REAL, POINTER, PRIVATE ::  DelInv(:,:), Deliau(:, :), Al(:, :)
LOGICAL, PRIVATE :: lPrivateVarAlloc = .FALSE.

REAL, POINTER, PRIVATE :: Au2G(:, :), del2G(:, :)                  !Size(nx)
REAL, POINTER, PRIVATE :: ainvl2G(:, :), ainvu2G(:, :), ainvd2G(:, :)            !Size(nx)
REAL, POINTER, PRIVATE ::  DelInv2G(:, :,:), Deliau2G(:, :, :), Al2G(:, :, :)
LOGICAL, PRIVATE :: lPrivateVar2GAlloc = .FALSE.

REAL, POINTER, PRIVATE :: buf(:)
REAL, POINTER, PRIVATE :: x0(:), s(:), b0(:)
REAL, POINTER, PRIVATE :: buf2g(:, :)
REAL, POINTER, PRIVATE :: x02g(:, :), s2g(:, :), b02g(:, :)

!Matrix Pointing
REAL, POINTER, PRIVATE :: Diag1g(:, :)
REAL, POINTER, PRIVATE :: RadOffDiag1g(:, :, :), AxOffDiag1g(:, :, :)
!Matrix Pointing
REAL, POINTER, PRIVATE :: Diag2g(:, :, :)
REAL, POINTER, PRIVATE :: RadOffDiag2g(:, :, :, :), AxOffDiag2g(:, :, :, :)

INTEGER, POINTER, PRIVATE :: NeighIdx(:, :)
INTEGER, POINTER, PRIVATE :: AxialPlaneMap(:)
CONTAINS

SUBROUTINE SetBiluSolverEnv(A)
TYPE(CMFDLS_TYPE) :: A

DelInv => A%BiLU%DelInv; Deliau => A%BiLU%DeliAu; Al => A%BiLU%Al
Diag1g => A%Diag; RadOffDiag1g => A%RadOffDiag; AxOffDiag1g => A%AxOffDiag
NeighIdx => A%NeighIdx; AxialPlaneMap => A%AxialPlaneMap

END SUBROUTINE

SUBROUTINE SetBilu2gSolverEnv(A)
TYPE(CMFDLS_TYPE) :: A

DelInv2g => A%BiLU%DelInv2g; Deliau2g => A%BiLU%DeliAu2g; Al2g => A%BiLU%Al2g
Diag2g => A%Diag2g; RadOffDiag2g => A%RadOffDiag2g; AxOffDiag2g => A%AxOffDiag2g
NeighIdx => A%NeighIdx; AxialPlaneMap => A%AxialPlaneMap

END SUBROUTINE

SUBROUTINE FreeBiluSolverEnv()
NULLIFY(DelInv, Deliau, Al)
NULLIFY(Diag1g, RadOffDiag1g, AxOffDiag1g)
NULLIFY(NeighIdx, AxialPlaneMap)
END SUBROUTINE


SUBROUTINE FreeBilu2gSolverEnv()
NULLIFY(DelInv2g, Deliau2g, Al2g)
NULLIFY(Diag2g, RadOffDiag2g, AxOffDiag2g)
NULLIFY(NeighIdx, AxialPlaneMap)
END SUBROUTINE

SUBROUTINE SolveBiLU(b, x)
USE Geom,             ONLY : Core
USE BasicOperation,   ONLY : CP_CA
REAL, POINTER :: B(:, :), X(:, :)
INTEGER :: l, k, k0, kp1
!REAL :: s(nxy), b0(nxy)

s(1:nxy) = 0.
x(1:nxy, myzbf:myzef) = 0.

DO k = myzbf, myzef
  k0 = AxialPlaneMap(k)
  DO l = 1, nxy
    b0(l) = b(l, k) - AxOffDiag1g(1, l, k) * s(l)
  ENDDO
  CALL SolveBiLU2d(k, b0, s)
  DO l = 1, nxy
    x(l, k) = s(l)
  ENDDO
ENDDO

DO k = myzef-1, myzbf, -1
  kp1 = k + 1
  k0 = AxialPlaneMap(k)
  DO l = 1, nxy
    b0(l) = x(l, kp1) * AxOffDiag1g(2, l, k)
  ENDDO
  CALL SolveBiLU2d(k, b0, s)
  DO l = 1, nxy
    x(l, k) = x(l, k) - s(l)
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE SolveBiLU2g(b, x)
USE Geom,             ONLY : Core
USE BasicOperation,   ONLY : CP_CA
USE MAT2x2OP
REAL, POINTER :: B(:, :, :), X(:, :, :)
INTEGER :: l, k, k0, kp1
!REAL :: s(nxy), b0(nxy)

s2g(1:2, 1:nxy) = 0.;
x(1:2, 1:nxy, myzbf:myzef) = 0.
DO k = myzbf, myzef
  k0 = AxialPlaneMap(k)
  DO l = 1, nxy
    b02g(:, l) = b(:, l, k) - SUBMATVECOP(AxOffDiag2g(:, 1, l, k), s2g(:, l))
  ENDDO
  CALL SolveBiLU2d2g(k, b02g, s2g)
  DO l = 1, nxy
    x(:, l, k) = s2g(:, l)
  ENDDO
ENDDO

DO k = myzef-1, myzbf, -1
  kp1 = k + 1
  k0 = AxialPlaneMap(k)
  DO l = 1, nxy
    b02g(:, l) = SUBMATVECOP(x(:, l, kp1), AxOffDiag2g(:, 2, l, k))
  ENDDO
  CALL SolveBiLU2d2g(k, b02g, s2g)
  DO l = 1, nxy
    x(:, l, k) = x(:, l, k) - s2g(:, l)
  ENDDO
ENDDO

!DO k = myzef-1, myzbf, -1
!  kp1 = k + 1
!  k0 = AxialPlaneMap(k)
!  DO l = 1, nxy
!    b0(l) = x(l, kp1) * AxOffDiag1g(2, l, k)
!  ENDDO
!  CALL SolveBiLU2d(k, b0, s)
!  DO l = 1, nxy
!    x(l, k) = x(l, k) - s(l)
!  ENDDO
!ENDDO
END SUBROUTINE

SUBROUTINE SolveBiLU2d(iz, b, x)
IMPLICIT NONE
INTEGER :: iz
REAL :: b(nxy), x(nxy)

REAL :: s1dl(nx), b01d(nx)
INTEGER :: i, j, k, k0, jp1, l, ls

k = iz; k0 = AxialPlaneMap(iz)
DO i = 1, nx
  s1dl(i) = 0
ENDDO
!forward solve
DO j = 1, ny
  DO i = nxbeg(j), nxend(j)
    l = node(i, j)
    b01d(i) = b(l) - RadOffDiag1g(3, l, k0) * s1dl(i)
  ENDDO
  CALL SolveBiLU1d(j, k, b01d, s1dl)
  DO i = nxbeg(j), nxend(j)
    l = node(i, j)
    x(l) = s1dl(i)
  ENDDO
ENDDO

!Backward solve
jp1 = ny
DO j = ny - 1, 1, -1
  DO i = nxend(j), nxbeg(j), -1
    l = node(i, j)
    ls = max(1, node(i,jp1))
    b01d(i) = x(ls) * RadOffDiag1g(1, l, k0)
  ENDDO
  CALL SolveBiLU1d(j, k, b01d, s1dl)
  DO i = nxend(j), nxbeg(j), -1
    l = node(i, j)
    x(l) = x(l) - s1dl(i)
  ENDDO
  jp1 = j
ENDDO
END SUBROUTINE

SUBROUTINE SolveBiLU2d2G(iz, b, x)
USE MAT2x2OP
IMPLICIT NONE
INTEGER :: iz
REAL :: b(2, nxy), x(2, nxy)

REAL :: s1dl(2, nx), b01d(2, nx)
INTEGER :: i, j, k, k0, jp1, l, ls

k = iz; k0 = AxialPlaneMap(iz)
DO i = 1, nx
  s1dl(:, i) = 0
ENDDO

!forward solve
DO j = 1, ny
  DO i = nxbeg(j), nxend(j)
    l = node(i, j)
    b01d(:, i) = b(:, l) - SUBMATVECOP(RadOffDiag2g(:, 3, l, k0), s1dl(:, i))
  ENDDO
  CALL SolveBiLU1d2g(j, k, b01d, s1dl)
  DO i = nxbeg(j), nxend(j)
    l = node(i, j)
    x(:, l) = s1dl(:, i)
  ENDDO
ENDDO

!Backward solve
jp1 = ny
DO j = ny - 1, 1, -1
  DO i = nxend(j), nxbeg(j), -1
    l = node(i, j)
    ls = max(1, node(i,jp1))
    b01d(:, i) = SUBMATVECOP(RadOffDiag2g(:, 1, l, k0), x(:, ls))
  ENDDO
  CALL SolveBiLU1d2g(j, k, b01d, s1dl)
  DO i = nxend(j), nxbeg(j), -1
    l = node(i, j)
    x(:, l) = x(:, l) - s1dl(:, i)
  ENDDO
  jp1 = j
ENDDO

END SUBROUTINE


SUBROUTINE SolveBiLU1d(irow, iz, b, x)
IMPLICIT NONE
INTEGER :: irow, iz
REAL :: b(nx), x(nx)

REAL :: y(nx), b1i
INTEGER :: ibeg, iend, im1, ip1
INTEGER :: i, j, k, l

ibeg = nxbeg(irow); iend = nxend(irow)
l = node(ibeg, irow); k = iz
i = ibeg
y(i) = DelInv(l, k) * b(i)
im1 = i
!Forward
DO i = ibeg+1, iend
  l = node(i, irow)
  b1i = b(i) - Al(l,k) * y(im1)
  y(i) = DelInv(l, k) * b1i
  im1 = i
ENDDO

!Backward
x(iend) = y(iend)
ip1 = iend
DO i = iend-1, ibeg, -1
  l = node(i, irow)
  x(i) = y(i) - Deliau(l, k) * x(ip1)
  ip1 = i
ENDDO
END SUBROUTINE

SUBROUTINE SolveBiLU1d2g(irow, iz, b, x)
USE MAT2x2OP
IMPLICIT NONE
INTEGER :: irow, iz
REAL :: b(2, nx), x(2, nx)

REAL :: y(2, nx), b1i(2)
INTEGER :: ibeg, iend, im1, ip1
INTEGER :: i, j, k, l

ibeg = nxbeg(irow); iend = nxend(irow)
l = node(ibeg, irow); k = iz
i = ibeg
y(:, i) = SUBMATVECOP(DelInv2g(:,l, k), b(:, i))
im1 = i
!Forward
DO i = ibeg+1, iend
  l = node(i, irow)
  b1i(:) = b(:,i) - SUBMATVECOP(Al2g(:, l,k), y(:, im1))
  y(:, i) = SUBMATVECOP(DelInv2g(:, l, k), b1i(:))
  im1 = i
ENDDO
!
!Backward
x(:, iend) = y(:, iend)
ip1 = iend

DO i = iend-1, ibeg, -1
  l = node(i, irow)
  x(:, i) = y(:, i) - SUBMATVECOP(Deliau2g(:, l, k), x(:, ip1))
  ip1 = i
ENDDO

END SUBROUTINE

SUBROUTINE SolveBiLU_OMP(b, x, idom)
USE LsRadDcmp_MOD,   ONLY : RadDcmp
USE BasicOperation,   ONLY : CP_CA
REAL, POINTER :: B(:, :), X(:, :)
INTEGER :: idom
INTEGER :: i, l, k, k0, kp1, nxylocal

INTEGER, POINTER :: List(:)
!REAL :: s(nxy), b0(nxy)

nxylocal = RadDcmp%nxylocal(idom)
List => RadDcmp%PinIdx(:, idom)

DO i = 1, nxylocal
  l = List(i)
  s(l) = 0
ENDDO
DO k = myzbf, myzef
  DO i = 1, nxylocal
    l = List(i)
    x(l, k) = 0
  ENDDO
ENDDO

DO k = myzbf, myzef
  k0 = AxialPlaneMap(k)
  DO i = 1, nxylocal
     l = List(i)
    b0(l) = b(l, k) - AxOffDiag1g(1, l, k) * s(l)
  ENDDO
!$OMP BARRIER
  CALL SolveBiLU2d_OMP(k, b0, s, idom)
  DO i = 1, nxylocal
    l = List(i)
    x(l, k) = s(l)
  ENDDO
!$OMP BARRIER
ENDDO

DO k = myzef-1, myzbf, -1
  kp1 = k + 1
  k0 = AxialPlaneMap(k)
  DO i = 1, nxylocal
     l = List(i)
    b0(l) = x(l, kp1) * AxOffDiag1g(2, l, k)
  ENDDO
!$OMP BARRIER
  CALL SolveBiLU2d_OMP(k, b0, s, idom)
  DO i = 1, nxylocal
     l = List(i)
    x(l, k) = x(l, k) - s(l)
  ENDDO
!$OMP BARRIER
ENDDO
END SUBROUTINE

SUBROUTINE SolveBiLU2d_OMP(iz, b, x, idom)
USE LsRAdDcmp_MOD, ONLY : RadDcmp
IMPLICIT NONE
INTEGER :: iz, idom
REAL :: b(nxy), x(nxy)


REAL :: s1dl(nx), b01d(nx)
!REAL :: b01d2(nx), s1d2(nx)
INTEGER :: i, j, k, k0, jp1, jp2, l, ls, ln


INTEGER :: jbeg, jend, lneigh
!INTEGER :: lneigh, idom0, idom1, idom
!INTEGER :: imod
!LOGICAL :: lbound

jbeg = RadDcmp%nybeg(idom); jend = RadDcmp%nyend(idom)

k = iz; k0 = AxialPlaneMap(iz)
DO i = 1, nx
  s1dl(i) = 0
ENDDO

IF(jbeg .NE. 1) THEN
  DO i = 1, nx
    s1dl(i) = 0
  ENDDO
  DO i = nxbeg(jbeg-1), nxend(jbeg-1)
    l=node(i,jbeg-1);
    ls = node(i, jbeg-2); ls = max(1, ls)
    b01d(i) = b(l)  -  RadOffDiag1g(3, l, k0) * b(ls) * DelInv(LS, K)
    s1dl(i) = b(l) / Diag1g(l, k)
  ENDDO
  CALL SolveBiLU1d(jbeg-1, k, b01d, s1dl)

ENDIF
!forward solve
DO j = jbeg, jend
  DO i = nxbeg(j), nxend(j)
    l = node(i, j)
    b01d(i) = b(l) - RadOffDiag1g(3, l, k0) * s1dl(i)
  ENDDO
  CALL SolveBiLU1d(j, k, b01d, s1dl)
  DO i = nxbeg(j), nxend(j)
    l = node(i, j)
    x(l) = s1dl(i)
  ENDDO
ENDDO

!$OMP BARRIER
 CALL DATACOPY_OMP(x, x0, idom)

!Backward solve
IF(jend .NE. ny) THEN
  jp1 = jend + 1; jp2 = jend + 2; j = jend
  DO i = nxend(jp1), nxbeg(jp1), -1
    l = node(i, jp1);  ls = max(1, node(i,jp2))
    ln = max(1, node(i, jp2+1))
    b01d(i) = RadOffDiag1g(1, l, k0) * (x0(ls))
  ENDDO
  CALL SolveBiLU1d(jp1, k, b01d, s1dl)
  DO i = nxend(jp1), nxbeg(jp1), -1
    l = node(i, jp1)
    s1dl(i) = x0(l) - s1dl(i)
  ENDDO
  DO i = nxend(j), nxbeg(j), -1
    l = node(i, j)
    lneigh = NeighIdx(1, l)
    b01d(i) = 0
    IF(lneigh .NE. -1) b01d(i) =   RadOffDiag1g(1, l, k0) * s1dl(i)
  ENDDO
  CALL SolveBiLU1d(j, k, b01d, s1dl)
  DO i = nxend(j), nxbeg(j), -1
    l = node(i, j)
    x(l) = x(l) - s1dl(i)
  ENDDO
ENDIF

jp1 = jend
DO j = jend - 1, jbeg, -1
  DO i = nxend(j), nxbeg(j), -1
    l = node(i, j)
    ls = max(1, node(i,jp1))
    b01d(i) = x(ls) * RadOffDiag1g(1, l, k0)
  ENDDO
  CALL SolveBiLU1d(j, k, b01d, s1dl)
  DO i = nxend(j), nxbeg(j), -1
    l = node(i, j)
    x(l) = x(l) - s1dl(i)
  ENDDO
  jp1 = j
ENDDO
!$OMP BARRIER
END SUBROUTINE


SUBROUTINE SolveBiLU1d_OMP(irow, iz, b, x, idom)
USE LsRAdDcmp_MOD, ONLY : RadDcmp
IMPLICIT NONE
INTEGER :: irow, iz, idom
REAL :: b(nx), x(nx)

REAL :: y(nx), b1i, x1, x2
INTEGER :: idom1, idom0, lneigh, lneigh2
INTEGER :: ibeg, iend, im1, im2, ip1, ip2
INTEGER :: i, j, k, l

ibeg = RadDcmp%nxbeg(irow, idom); iend = RadDcmp%nxend(irow, idom)
l = node(ibeg, irow); k = iz
i = ibeg
IF(ibeg .EQ. nxbeg(irow)) THEN
   y(i) = DelInv(l, k) * b(i)
ELSE
  im2 = ibeg - 2; im1 = ibeg - 1
  lneigh = NeighIdx(2, l); lneigh2 = NeighIdx(2, lneigh)
  x2=b(im2) * DelInv(lneigh2, k) * Al(lneigh, k)
  x1=(b(im1) -x2)* DelInv(lneigh, k)
  b1i = b(i) - Al(l,k) * x1
  y(i) = DelInv(l, k) * b1i
ENDIF
im1 = i
im2 = 0

!Forward
DO i = ibeg+1, iend
  l = node(i, irow);
  b1i = b(i) - Al(l,k) * y(im1)
  y(i) = DelInv(l, k) * b1i
  im2 = im1
  im1 = i
ENDDO

IF(iend .EQ. nxend(irow)) THEN
  x(iend) = y(iend)
ELSE
  l = node(iend, irow)
  ip1 = iend + 1; ip2 = iend + 2
  lneigh = NeighIdx(4, l); lneigh2 = NeighIdx(4, lneigh)
  x2= Deliau(lneigh, k) * y(ip2)
  x1 = y(ip1) -x2
  x(iend) = y(iend) - Deliau(l, k) * x1
ENDIF
ip1 = iend
DO i = iend-1, ibeg, -1
  l = node(i, irow)
  x(i) = y(i) - Deliau(l, k) * x(ip1)
  ip1 = i
ENDDO

END SUBROUTINE

SUBROUTINE MakeBiLU(Core, A, X, ng_,  PE)
TYPE(CoreInfo_Type) :: Core
TYPE(CMFDLS_TYPE) :: A(ng_)
TYPE(PE_TYPE) :: PE
REAL, POINTER :: X(:, :, :)
INTEGER :: ng_
INTEGER :: ig, iz

IF(.NOT. lSetGeom) CALL SetBiLUGeom(Core, ng_, PE)
IF(.NOT. lPrivateVarAlloc) CALL AllocBiLUVar()

!CALL SetMPIBiLU(A, X, ng_, 1, PE)
DO ig = 1, ng_
  !Allocate The BiLU matrix variable
  IF(.NOT. A(ig)%BiLU%lAlloc) CALL AllocBiLUMat(A(ig)%BiLU)
  !Variables Pointing
  DelInv => A(ig)%BiLU%DelInv; Deliau =>  A(ig)%BiLU%DeliAu; Al => A(ig)%BiLU%Al
  Diag1g => A(ig)%Diag; RadOffDiag1g => A(ig)%RadOffDiag; AxOffDiag1g => A(ig)%AxOffDiag
  NeighIdx => A(ig)%NeighIdx; AxialPlaneMap => A(ig)%AxialPlaneMap

  !2D incomplete LU
  DO iz = myzbf, myzef
    !Factorize 2-D(Plane)
    CALL FactorizeILU2D(iz)
  ENDDO
  !
  NULLIFY(DelInv, Deliau, Al)
  NULLIFY(Diag1g, RadOffDiag1g, AxOffDiag1g)
  NULLIFY(NeighIdx, AxialPlaneMap)
ENDDO
!CALL SetMPIBiLU(A, X, ng_, 2, PE)
END SUBROUTINE

SUBROUTINE MakeBiLU2G(Core, A, X, ng_, PE)
USE IOUTIL,       ONLY : TERMINATE
TYPE(CoreInfo_TYPE) :: Core
TYPE(CMFDLS_TYPE) :: A
TYPE(PE_TYPE) :: PE
REAL, POINTER :: X(:, :, :)
INTEGER:: ng_

INTEGER :: iz

IF(.NOT. A%l2G) CALL TERMINATE('MakeBiLU2G: Matrix does not have 2G structure')

IF(.NOT. lSetGeom) CALL SetBiLUGeom(Core, ng_, PE)
IF(.NOT. lPrivateVar2gAlloc) CALL AllocBiLUVar2G()
IF(.NOT. A%BiLU%lAlloc) CALL AllocBilUMat2G(A%BiLU)
A%BiLU%l2G = .TRUE.
!Variables Pointing
DelInv2g => A%BiLU%DelInv2g; Deliau2g => A%BiLU%DeliAu2g; Al2g => A%BiLU%Al2g
Diag2g => A%Diag2g; RadOffDiag2g => A%RadOffDiag2g; AxOffDiag2g => A%AxOffDiag2g
NeighIdx => A%NeighIdx; AxialPlaneMap => A%AxialPlaneMap
DO iz = myzbf, myzef
  CALL FactorizeILU2D2G(iz)
ENDDO
NULLIFY(DelInv2g, Deliau2g, Al2g)
NULLIFY(Diag2g, RadOffDiag2g, AxOffDiag2g)
NULLIFY(NeighIdx, AxialPlaneMap)

END SUBROUTINE

SUBROUTINE SetMPIBiLU(A, X, ng_, imod, PE)
#ifdef MPI_ENV
USE MPIComm_MOD, ONLY : GetNeighDat
#endif
IMPLICIT NONE
TYPE(CMFDLS_TYPE) :: A(ng_)
TYPE(PE_TYPE) :: PE
REAL, POINTER :: X(:, :, :)
REAL :: temp
INTEGER :: ng_, imod
INTEGER :: ig, iz, ixy
IF(imod .EQ. 1) THEN
  DO ig = 1, ng_
#ifdef MPI_ENV
    CALL GetNeighDat(X(1:nxy, myzbf-1:myzef+1, ig), nxy, myzbf, myzef,   &
                           PE%myCmfdRank, PE%nproc, PE%MPI_CMFD_COMM)
#endif
    DO ixy = 1, nxy
      iz = myzbf
      temp = 0
      IF(iz .NE. 1) temp = A(ig)%AxOffDiag(1, ixy, iz) * X(ixy, myzbf-1, ig) / X(ixy, myzbf, ig)
      A(ig)%Diag(ixy, iz) = A(ig)%Diag(ixy, iz) + temp
      iz = myzef
      temp = 0
      IF(iz .NE. nzfm) temp = A(ig)%AxOffDiag(2, ixy, iz) * X(ixy, myzef+1, ig) / X(ixy, myzef, ig)
      A(ig)%Diag(ixy, iz) = A(ig)%Diag(ixy, iz) + temp
    ENDDO
  ENDDO
ELSE
  DO ig = 1, ng_
#ifdef MPI_ENV
    CALL GetNeighDat(X(1:nxy, myzbf-1:myzef+1, ig), nxy, myzbf, myzef, PE%myCmfdRank, PE%nproc,   &
                    PE%MPI_CMFD_COMM)
#endif
    DO ixy = 1, nxy
      iz = myzbf
      temp = 0
      IF(iz .NE. 1) temp = A(ig)%AxOffDiag(1, ixy, iz) * X(ixy, myzbf-1, ig) / X(ixy, myzbf, ig)
      A(ig)%Diag(ixy, iz) = A(ig)%Diag(ixy, iz) - temp
      iz = myzef
      temp = 0
      IF(iz .NE. nzfm) temp = A(ig)%AxOffDiag(2, ixy, iz) * X(ixy, myzef+1, ig) / X(ixy, myzef, ig)
      A(ig)%Diag(ixy, iz) = A(ig)%Diag(ixy, iz) - temp
    ENDDO
  ENDDO
ENDIF

END SUBROUTINE

SUBROUTINE FactorizeILU2D(iz)   !FacILU
IMPLICIT NONE
INTEGER :: iz

INTEGER :: ixy
INTEGER :: i, j, k, l
INTEGER :: jm1, ln, lnm1, lnp1, k0

!REAL, POINTER, PRIVATE :: Au(:), del(:)                  !Size(nx)
k = iz; k0 = AxialPlaneMap(k)
j = 1
DO i = nxbeg(j), nxend(j)
  l = node(i, j)
  Del(i) = diag1g(l, k)
  Al(l, k) = RadOffDiag1g(2, l, k0)
  Au(i) = RadOffDiag1g(4, l, k0)
ENDDO
DO j = 2, ny
  jm1 = j - 1
  !Factorize Del_j
  CALL FactorizeILU1D(jm1, k)
  !Inverse of Factorize 1-D(line)
  CALL ApproxBlockInv1D(jm1, iz)
  !D_j+1 = a_j+1 - l_j+1 * inv(D_j) * u_j
  DO i = nxbeg(j), nxend(j)
    l = node(i, j)
    IF(i .GE. nxbeg(jm1) .AND. i .LE. nxend(jm1)) THEN
      ln=node(i,jm1)
      Del(i) = Diag1g(l, k) - RadOffDiag1g(3, l, k0) * AinvD(i) * RadOffDiag1g(1, ln, k0)
    ELSE
      Del(i) = Diag1g(l, k)
    ENDIF

    !Lower Diagonal Part
    Al(l, k) = RadOffDiag1g(2, l, k0)

    IF(i .NE. nxbeg(j) .AND. (i-1) .GE. nxbeg(jm1) .AND. (i-1) .LE. nxend(jm1)) THEN
      !IF(lnm1 .NE. 0)
      lnm1 = node(i-1, jm1)
      Al(l, k) = Al(l, k) - RadOffDiag1g(3, l, k0) * Ainvl(i) * RadOffDiag1g(1, lnm1, k0)
    ENDIF

    !Upper Diagoanl Part
    Au(i) = RadOffDiag1g(4, l, k0)
    IF(i .NE. nxend(j) .AND. (i+1) .GE. nxbeg(jm1) .AND. (i+1).le.nxend(jm1)) then
      !IF(lnp1 .NE. 0)
      lnp1 = node(i+1, jm1)
      Au(i) = Au(i) - RadOffDiag1g(3, l, k0) * Ainvu(i) * RadOffDiag1g(1, lnp1, k0)
    ENDIF
  ENDDO
  CONTINUE
ENDDO
CALL FactorizeILU1D(ny, k)
END SUBROUTINE

SUBROUTINE FactorizeILU2D2G(iz)   !FacILU
USE MAT2x2OP
IMPLICIT NONE
INTEGER :: iz

INTEGER :: ixy
INTEGER :: i, j, k, l
INTEGER :: jm1, ln, lnm1, lnp1, k0

!REAL, POINTER, PRIVATE :: Au(:), del(:)                  !Size(nx)
k = iz; k0 = AxialPlaneMap(k)
j = 1

DO i = nxbeg(j), nxend(j)
  l = node(i, j)
  Del2G(:, i) =Diag2g(:, l, k)
  Al2g(:, l, k) = RadOffDiag2g(:, 2, l, k0)
  Au2g(:, i) = RadOffDiag2g(:, 4, l, k0)
ENDDO

DO j = 2, ny
  jm1 = j - 1
  !Factorize Del_j
  CALL FactorizeILU1D2G(jm1, k)
  !Inverse of Factorize 1-D(line)
  CALL ApproxBlockInv1d2g(jm1, iz)
  !D_j+1 = a_j+1 - l_j+1 * inv(D_j) * u_j
  DO i = nxbeg(j), nxend(j)
    l = node(i, j)
    IF(i .GE. nxbeg(jm1) .AND. i .LE. nxend(jm1)) THEN
      ln=node(i,jm1)
      Del2G(:, i) = Diag2g(:, l, k) - MULTI_MATOP(RadOffDiag2g(:, 3, l, k0), AinvD2g(:, i), RadOffDiag2g(:, 1, ln, k0))
    ELSE
      Del2G(:, i) = Diag2g(:, l, k)
    ENDIF

    !Lower Diagonal Part
    Al2G(:, l, k) = RadOffDiag2g(:, 2, l, k0)
    IF(i .NE. nxbeg(j) .AND. (i-1) .GE. nxbeg(jm1) .AND. (i-1) .LE. nxend(jm1)) THEN
      !IF(lnm1 .NE. 0)
      lnm1 = node(i-1, jm1)
      Al2G(:, l, k) = Al2G(:, l, k) - MULTI_MATOP(RadOffDiag2g(:, 3, l, k0), Ainvl2g(:, i), RadOffDiag2g(:, 1, lnm1, k0))
    ENDIF

    !Upper Diagoanl Part
    Au2g(:, i) = RadOffDiag2g(:, 4, l, k0)
    IF(i .NE. nxend(j) .AND. (i+1) .GE. nxbeg(jm1) .AND. (i+1).le.nxend(jm1)) then
      !IF(lnp1 .NE. 0)
      lnp1 = node(i+1, jm1)
      Au2G(:, i) = Au2G(:, i) - MULTI_MATOP(RadOffDiag2g(:, 3, l, k0), Ainvu2G(:, i), RadOffDiag2g(:, 1, lnp1, k0))
    ENDIF
  ENDDO
ENDDO
CALL FactorizeILU1D2g(ny, k)
END SUBROUTINE

SUBROUTINE FactorizeILU1D(irow, iz)
INTEGER :: irow, iz
INTEGER :: i, j, k, l
INTEGER :: lm1, im1

REAL :: ald1
k = iz
i = nxbeg(irow); l = node(i, irow)

DelInv(l, k) = 1._8 / Del(i)
!Calculate Inv(del) * u for later use in the back sub
Deliau(l, k) = DelInv(l, k) * Au(i)

im1 = i
DO i = nxbeg(irow) + 1, nxend(irow)
  lm1 = l                     !previous row index
  l = node(i, irow)           !current row index
  ald1 = al(l, k) * DelInv(lm1, k)
  Del(i) = Del(i) - ald1 * Au(im1)
  DelInv(l, k) = 1 / Del(i)
  !Calculate inv(del)*u for later use in the back sub.
  Deliau(l, k) = DelInv(l, k) * Au(i)
  im1 = i
ENDDO

ENDSUBROUTINE

SUBROUTINE FactorizeILU1D2G(irow, iz)
USE MAT2x2OP
INTEGER :: irow, iz
INTEGER :: i, j, k, l
INTEGER :: lm1, im1

REAL :: ald1(4)
k = iz
i = nxbeg(irow); l = node(i, irow)

DelInv2G(:, l, k) = MATINV(Del2g(:, i))
!Calculate Inv(del) * u for later use in the back sub
Deliau2g(:, l, k) = SUBMATOP(DelInv2g(:, l, k), Au2g(:, i))

im1 = i
DO i = nxbeg(irow) + 1, nxend(irow)
  lm1 = l                     !previous row index
  l = node(i, irow)           !current row index
  ald1(:) = SUBMATOP(al2g(:, l, k), DelInv2g(:, lm1, k))
  Del2g(:, i) = Del2g(:, i) - SUBMATOP(ald1, Au2g(:, im1))
  DelInv2g(:, l, k) = MATINV(Del2g(:, i))
  !Calculate inv(del)*u for later use in the back sub.
  Deliau2g(:, l, k) = SUBMATOP(DelInv2g(:, l, k), Au2g(:, i))
  im1 = i
ENDDO
!im1 = i
!DO i = nxbeg(irow) + 1, nxend(irow)
!  lm1 = l                     !previous row index
!  l = node(i, irow)           !current row index
!  ald1 = al(l, k) * DelInv(lm1, k)
!  Del(i) = Del(i) - ald1 * Au(im1)
!  DelInv(l, k) = 1 / Del(i)
!  !Calculate inv(del)*u for later use in the back sub.
!  Deliau(l, k) = DelInv(l, k) * Au(i)
!  im1 = i
!ENDDO

ENDSUBROUTINE


SUBROUTINE ApproxBlockInv1d(irow, iz)
INTEGER :: irow, iz
INTEGER :: i, j, k, l
INTEGER :: lp1, mp1, ix
REAL :: al1, au1

k = iz

!Last Row
ix = nxend(irow); l=node(ix, irow)
AInvd(ix) = DelInv(l, k)

DO i = nxend(irow) - 1, nxbeg(irow), -1
  lp1 = l;  mp1 = ix
  l = node(i, irow); ix = ix - 1
  !Lower Index
  al1 = AinvD(mp1) * al(lp1, k)
  AInvl(mp1) = - al1 * DelInv(l, k)
  !Upper Index
  au1 = DelInv(l, k) * Au(ix)
  AinvU(ix) = - au1 * AinvD(mp1)
  !Diagonal Part
  AinvD(ix) = DelInv(l, k) - au1 * Ainvl(mp1)
ENDDO

END SUBROUTINE

SUBROUTINE ApproxBlockInv1d2G(irow, iz)
USE MAT2x2OP
INTEGER :: irow, iz
INTEGER :: i, j, k, l
INTEGER :: lp1, mp1, ix
REAL :: al1(4), au1(4)

k = iz

!Last Row
ix = nxend(irow); l=node(ix, irow)
AInvd2g(:, ix) = DelInv2g(:, l, k)

DO i = nxend(irow) -1, nxbeg(irow), -1
  lp1 = l;  mp1 = ix
  l = node(i, irow); ix = ix - 1
  !Lower Index
  Al1(:) = SUBMATOP(AinvD2g(:, mp1), Al2g(:, lp1, k))
  Ainvl2g(:, mp1) = - SUBMATOP(Al1(:), DelInv2g(:, l, k))
  !Upper Indxe
  Au1(:) = SUBMATOP(DelInv2g(:, l, k), Au2g(:, ix))
  AinvU2g(:, ix) = -SUBMATOP(Au1(:), AinvD2g(:, mp1))
  !Diagonal Part
  AinvD2g(:, ix) = DelInv2g(:, l, k) - SUBMATOP(Au1(:), Ainvl2g(:, mp1))
ENDDO

END SUBROUTINE

SUBROUTINE SetBiLUGeom(Core, ng_, PE)
TYPE(CoreInfo_TYPE) :: Core
TYPE(PE_TYPE) :: PE
INTEGER :: ng_
TYPE(Pin_TYPE), POINTER :: PIN(:)
INTEGER :: ix, iy, ixy

PIN => Core%PIN
nx = Core%nx; ny = Core%ny; nxy = Core%nxy
!ng = ng_
CALL Dmalloc0(node, 0, nx+1, 0, ny+1)
CALL Dmalloc0(nxbeg, 0, ny+1); CALL Dmalloc0(nxend, 0, ny+1)
DO ixy = 1, Core%nxy
  ix = Pin(ixy)%ix; iy = Pin(ixy)%iy
  node(ix, iy) = ixy
ENDDO

DO iy = 1, ny
  DO ix = 1, nx
    IF(NODE(ix, iy) .GT. 0) EXIT
  ENDDO
  nxbeg(iy) = ix;
  DO ix = nx, 1, -1
    IF(NODE(ix, iy) .GT. 0) EXIT
  ENDDO
  nxend(iy) = ix
ENDDO

nz = PE%nz; nzfm = PE%nzfm
myzb = PE%myzb; myze = PE%myzef
myzbf = PE%myzbf; myzef = PE%myzef

lSetGeom = .TRUE.
NULLIFY(PIN)
END SUBROUTINE

SUBROUTINE AllocBiLUVar()
IMPLICIT NONE
!REAL, POINTER, PRIVATE :: Au(:), del(:)                  !Size(nx)
!REAL, POINTER, PRIVATE :: ainvl(:), ainvu(:), ainvd(:)            !Size(nx)

CALL Dmalloc(Au, nx);      CALL Dmalloc(Del, nx)
CALL Dmalloc(Ainvl, nx);   CALL Dmalloc(Ainvu, nx)
CALL Dmalloc(AinvD, nx);   CALL Dmalloc(Buf, nxy)
CALL Dmalloc(S, nxy);   CALL Dmalloc(B0, nxy)
CALL Dmalloc(X0, nxy)
lPrivateVarAlloc = .TRUE.
END SUBROUTINE

SUBROUTINE AllocBiLUVar2G()
IMPLICIT NONE
CALL Dmalloc(Au2G, 4, nx);      CALL Dmalloc(Del2G, 4, nx)
CALL Dmalloc(Ainvl2G, 4, nx);   CALL Dmalloc(Ainvu2G, 4, nx)
CALL Dmalloc(AinvD2G, 4, nx);
!CALL Dmalloc(Buf2G, 2, nxy)
CALL Dmalloc(S2g, 2, nxy);   CALL Dmalloc(B02g, 2, nxy)
!CALL Dmalloc(X0, nxy)
lPrivateVar2GAlloc = .TRUE.
END SUBROUTINE

SUBROUTINE AllocBiLUMat(BiLU)
TYPE(BiLU_TYPE) :: BiLU
CALL Dmalloc0(BiLU%DelInv, 1, nxy, myzbf, myzef)
CALL Dmalloc0(BiLU%Deliau, 1, nxy, myzbf, myzef)
CALL Dmalloc0(BiLU%Al, 1, nxy, myzbf, myzef)
BiLU%lAlloc = .TRUE.
END SUBROUTINE

SUBROUTINE AllocBiLUMat2G(BiLU)
TYPE(BiLU_TYPE) :: BiLU
CALL Dmalloc0(BiLU%DelInv2G, 1, 4, 1, nxy, myzbf, myzef)
CALL Dmalloc0(BiLU%Deliau2G, 1, 4, 1, nxy, myzbf, myzef)
CALL Dmalloc0(BiLU%Al2G, 1, 4, 1, nxy, myzbf, myzef)
BiLU%lAlloc = .TRUE.
END SUBROUTINE

SUBROUTINE DATASYNC_OMP(x, idom)
USE LsRadDcmp_MOD,   ONLY : RadDcmp
IMPLICIT NONE
REAL :: x(nxy)
INTEGER :: idom

INTEGER, POINTER :: List(:)
INTEGER :: nxylocal
INTEGER :: i, ixy
nxylocal = RadDcmp%nxylocal(idom)
List => RadDcmp%PinIdx(:, idom)

DO i = 1, nxylocal
  ixy=List(i)
  Buf(ixy) = x(ixy)
ENDDO
!$OMP BARRIER
DO i = 1, nxy
  x(ixy)=Buf(ixy)
ENDDO
!$OMP BARRIER
END SUBROUTINE

SUBROUTINE DATACOPY_OMP(dat1, dat2, idom)
USE LsRadDcmp_MOD,   ONLY : RadDcmp
IMPLICIT NONE
REAL :: dat1(nxy), dat2(nxy)
INTEGER :: idom

INTEGER, POINTER :: List(:)
INTEGER :: nxylocal
INTEGER :: i, ixy
nxylocal = RadDcmp%nxylocal(idom)
List => RadDcmp%PinIdx(:, idom)

DO i = 1, nxylocal
  ixy=List(i)
  dat2(ixy) = dat1(ixy)
ENDDO
!$OMP BARRIER

END SUBROUTINE

END MODULE
