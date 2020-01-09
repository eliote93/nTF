MODULE LsRadDcmp_MOD
USE PARAM
USE TYPEDEF,           ONLY : RadDcmp_type
IMPLICIT NONE
TYPE(RadDcmp_type), SAVE :: RadDcmp
!REAL, POINTER, PRIVATE :: RadCoeffDat(:, :)
!REAL, POINTER :: RadDB2(:, :), Shat(:)
INTEGER, PRIVATE :: nxy
!LOGICAL, PRIVATE :: lAlloc = .FALSE.
!LOGICAL, PRIVATE :: lAllocDb2 = .FALSE.
LOGICAL, PRIVATE :: lInit = .FALSE.
CONTAINS

SUBROUTINE SetRadDcmp(Core, PE)
USE TYPEDEF,        ONLY : CoreInfo_Type,        Pin_Type,        Asy_Type,         PE_TYPE
USE BasicOperation, ONLY : CP_CA
USE UtilFunction,   ONLY : GetRangeDecomp,       GetRadDecomp_SA
USE ALLOCS
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(PE_TYPE) :: PE

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Asy_Type), POINTER :: Asy(:)

INTEGER :: ixya
INTEGER :: i, j, k, k0, ix, iy, ixy, ixy0
INTEGER :: ndivx, ndivy, nxy0
INTEGER :: nxylocal(100), XDomRange(2, 100), YDomRange(2, 100)
INTEGER, POINTER :: node(:, :), LineInfo(:)
LOGICAL :: lZigZag

IF(lInit) RETURN
nxy = Core%nxy
Pin => Core%Pin; Asy => Core%Asy

IF(PE%nCmfdThread .GT. nxy) THEN
  PE%nCmfdThread = 1
ENDIF

IF(PE%nCmfdThread  .GT. Core%ny) THEN
  PE%nCmfdThread = 1
ENDIF


ALLOCATE(node(Core%nx, Core%ny))
ALLOCATE(LineInfo(Core%ny))
CALL CP_CA(node, 0, Core%nx, Core%ny)
CALL CP_CA(LineInfo, 0, Core%ny)
DO ixy = 1, Core%nxy
  ix = Pin(ixy)%ix; iy = Pin(ixy)%iy
  node(ix, iy) = ixy
  LineInfo(iy) = LineInfo(iy) + 1
ENDDO
lZigZag = .FALSE.
DO iy = 2, Core%ny
  IF(LineInfo(iy) .NE. LineInfo(1)) THEN
    lZigZag = .TRUE.
    EXIT
  ENDIF
ENDDO

ndivx = 1; ndivy = PE%nCmfdThread
DO i = 1, ndivx
  CALL GetRangeDecomp(1, Core%nx, ndivx, i-1, XdomRange(1, i), XdomRange(2, i))
ENDDO
DO i = 1, ndivy
  CALL GetRangeDecomp(1, Core%ny, ndivy, i-1, YdomRange(1, i), YdomRange(2, i))
ENDDO
nxylocal = 0

IF(lZigZag) CALL GetRadDecomp_SA(YdomRange, LineInfo,Core%ny, ndivy)
k = 0
DO iy = 1, ndivy
  DO ix = 1, ndivx
    k = k+1
    DO j = YdomRange(1, iy), YdomRange(2, iy)
      DO i = XdomRange(1, ix), XdomRange(2, ix)
        ixy = node(i, j)
        if(ixy .EQ. 0) CYCLE
        nxylocal(k) = nxylocal(k) + 1
      ENDDO
    ENDDO
  ENDDO
ENDDO

!Set Dimension
k0 = k; nxy0 = 0; k = 0
do i = 1, k0
  IF(nxylocal(i) .EQ. 0) CYCLE
  k = k + 1
  nxy0 = max(nxy0, nxylocal(i))
ENDDO
RadDcmp%ndom = k

CALL DMALLOC(RadDcmp%PinIdx, nxy0, k)
CALL DMALLOC0(RadDcmp%PinDom, -1, Core%nxy)
!
RadDcmp%PinDom(-1) = 0
k0 = 0; k = 0
DO iy = 1, ndivy
  DO ix = 1, ndivx
    k0 = k0 + 1;
    IF(nxylocal(k0) .EQ. 0) CYCLE
    k = k + 1; RadDcmp%nxylocal(k) = nxylocal(k0)
    ixy0 = 0
    DO j = YdomRange(1, iy), YdomRange(2, iy)
      DO i = XdomRange(1, ix), XdomRange(2, ix)
        ixy = node(i, j)
        if(ixy .EQ. 0) CYCLE
        ixy0 = ixy0 + 1
        RadDcmp%PinIdx(ixy0, k) = ixy
        RadDcmp%PinDom(ixy) = k
      ENDDO
    ENDDO
  ENDDO
ENDDO

CALL DMALLOC(RadDcmp%nxbeg, Core%ny, RadDcmp%ndom)
CALL DMALLOC(RadDcmp%nxend, Core%ny, RadDcmp%ndom)
CALL DMALLOC(RadDcmp%nybeg, RadDcmp%ndom)
CALL DMALLOC(RadDcmp%nyend, RadDcmp%ndom)
RadDcmp%nxbeg=100000; RadDcmp%nxend=0
RadDcmp%nybeg=100000; RadDcmp%nyend=0
DO iy = 1, Core%ny
  DO ix = 1, Core%nx
    IF(node(ix, iy) .EQ. 0) CYCLE
    ixy = node(ix, iy)
    i = RadDcmp%PinDom(ixy)
    RadDcmp%nxbeg(iy, i) = MIN(RadDcmp%nxbeg(iy, i), ix)
    RadDcmp%nxend(iy, i) = MAX(RadDcmp%nxend(iy, i), ix)
    RadDcmp%nybeg(i) = MIN(RadDcmp%nybeg(i), iy)
    RadDcmp%nyend(i) = MAX(RadDcmp%nyend(i), iy)
  ENDDO
ENDDO
DEALLOCATE(node)
NULLIFY(Asy, Pin)
linit = .TRUE.
END SUBROUTINE

!SUBROUTINE SetRadDcmpLs(Core, A, X, ig, iz, iz0, PE, mod)
!USE PARAM
!USE TYPEDEF,        ONLY : CoreInfo_Type, CMFDLS_Type,  PE_TYPE,         &
!                           Pin_Type
!USE BasicOperation, ONLY : CP_CA
!IMPLICIT NONE
!TYPE(CoreInfo_Type) :: Core
!TYPE(CMFDLS_Type) :: A
!TYPE(PE_TYPE) :: PE
!REAL :: x(nxy)
!INTEGER :: ig, iz, iz0, mod
!
!TYPE(Pin_Type), POINTER :: Pin(:)
!REAL :: Temp
!INTEGER :: ixy, idom, ineigh
!INTEGER :: i, j
!INTEGER :: ndom
!
!IF(.NOT. lAlloc) THEN
!  ALLOCATE(RadCoeffDat(4, nxy))
!  lAlloc = .TRUE.
!ENDIF
!
!ndom = RadDcmp%ndom
!Pin => Core%Pin
!IF(mod .EQ. 1) THEN
!  CALL CP_CA(RadCoeffDat, 0._8, 4, nxy)
!  DO idom = 1, ndom
!    DO i = 1, RadDcmp%nxylocal(idom)
!      ixy = RadDcmp%PinIdx(i, idom)
!      DO j = 1, 4
!       ! IF(idom .EQ. 3 .AND. j==2) CYCLE
!        ineigh = Pin(ixy)%NeighIdx(j)
!        IF(ineigh .LT. 1) CYCLE
!        IF(RadDcmp%PinDom(ineigh) .EQ. idom) CYCLE
!        RadCoeffDat(j, ixy) = A%RadOffDiag(j, ixy, iz0)
!        Temp = A%RadOffDiag(j, ixy, iz0) * X(ineigh) / x(ixy)
!        !temp = 0
!        A%RadOffDiag(j, ixy, iz0) = 0
!        A%Diag(ixy, iz) = A%Diag(ixy, iz) + Temp
!      ENDDO
!    ENDDO
!  ENDDO
!ELSE
!  DO idom = 1, ndom
!    DO i = 1, RadDcmp%nxylocal(idom)
!      ixy = RadDcmp%PinIdx(i, idom)
!      DO j = 1, 4
!        !IF(idom .EQ. 3 .AND. j==2) CYCLE
!        ineigh = Pin(ixy)%NeighIdx(j)
!        IF(ineigh .LT. 1) CYCLE
!        IF(RadDcmp%PinDom(ineigh) .EQ. idom) CYCLE
!        Temp = RadCoeffDat(j, ixy) * X(ineigh) / x(ixy)
!       ! temp = 0
!        A%Diag(ixy, iz) = A%Diag(ixy, iz) - Temp
!        A%RadOffDiag(j, ixy, iz0) = RadCoeffDat(j, ixy)
!      ENDDO
!    ENDDO
!  ENDDO
!ENDIF
!CONTINUE
!END SUBROUTINE
!
!SUBROUTINE SetRadDcmpBiLU2D(Core, RadOffDiag1g, mod)
!USE PARAM
!USE TYPEDEF,        ONLY : CoreInfo_Type, CMFDLS_Type,  PE_TYPE,         &
!                           Pin_Type
!USE BasicOperation, ONLY : CP_CA
!IMPLICIT NONE
!TYPE(CoreInfo_Type) :: Core
!REAL :: RadOffDiag1g(4, nxy)
!INTEGER :: mod
!
!TYPE(Pin_Type), POINTER :: Pin(:)
!REAL :: Temp
!INTEGER :: ixy, idom, ineigh
!INTEGER :: i, j
!INTEGER :: ndom
!
!IF(.NOT. lAlloc) THEN
!  ALLOCATE(RadCoeffDat(4, nxy))
!  lAlloc = .TRUE.
!ENDIF
!
!ndom = RadDcmp%ndom
!Pin => Core%Pin
!IF(mod .EQ. 1) THEN
!  CALL CP_CA(RadCoeffDat, 0._8, 4, nxy)
!  DO idom = 1, ndom
!    DO i = 1, RadDcmp%nxylocal(idom)
!      ixy = RadDcmp%PinIdx(i, idom)
!      DO j = 1, 4
!        ineigh = Pin(ixy)%NeighIdx(j)
!        IF(ineigh .LT. 1) CYCLE
!        IF(RadDcmp%PinDom(ineigh) .EQ. idom) CYCLE
!        RadCoeffDat(j, ixy) = RadOffDiag1g(j, ixy)
!        RadOffDiag1g(j, ixy) = 0
!      ENDDO
!    ENDDO
!  ENDDO
!ELSE
!  DO idom = 1, ndom
!    DO i = 1, RadDcmp%nxylocal(idom)
!      ixy = RadDcmp%PinIdx(i, idom)
!      DO j = 1, 4
!        ineigh = Pin(ixy)%NeighIdx(j)
!        IF(ineigh .LT. 1) CYCLE
!        IF(RadDcmp%PinDom(ineigh) .EQ. idom) CYCLE
!        RadOffDiag1g(j, ixy) = RadCoeffDat(j, ixy)
!      ENDDO
!    ENDDO
!  ENDDO
!ENDIF
!CONTINUE
!END SUBROUTINE
!
!SUBROUTINE UpdtDB2(A, x, b, iz1, iz2, mod)
!USE PARAM
!USE TYPEDEF,   ONLY :  CMFDLS_Type
!USE BasicOperation, ONLY : CP_CA
!USE Allocs
!IMPLICIT NONE
!TYPE(CmfdLs_Type) :: A
!REAL :: x(nxy, iz1:iz2), b(nxy, iz1:iz2)
!INTEGER :: iz1, iz2, mod
!
!REAL :: temp
!INTEGER :: iz, iz0, idom, ixy, ineigh, ibd, nbd
!
!IF(.NOT. lAllocDb2) THEN
!  CALL Dmalloc0(RadDB2, 1, nxy, iz1,iz2)
!  CALL Dmalloc(Shat, nxy)
!  lAllocDb2 = .TRUE.
!ENDIF
!
!CALL CP_CA(RadDB2(1:nxy, iz1:iz2), 0._8, nxy, iz2-iz1+1)
!nbd =4
!DO iz = iz1, iz2
!  iz0 = A%AxialPlaneMap(iz)
!  DO ixy = 1, nxy
!    temp = 0
!    temp = A%Diag(ixy, iz) * x(ixy, iz)
!    DO ibd = 1, nbd
!      ineigh = A%NeighIdx(ibd, ixy)
!
!      IF(ineigh .LE. 0) CYCLE
!      !temp = temp + A%RadOffDiag(ibd, ixy, iz0) * x(ineigh, iz)
!    ENDDO
!    temp =  temp / x(ixy, iz)
!    RadDB2(ixy, iz) = temp
!  ENDDO
!ENDDO
!END SUBROUTINE
!
!SUBROUTINE SetRadDcmpDB2(Core, s, iz)
!USE PARAM
!USE TYPEDEF,        ONLY : CoreInfo_Type, CMFDLS_Type,  PE_TYPE,         &
!                           Pin_Type
!USE BasicOperation, ONLY : CP_CA, CP_VA
!IMPLICIT NONE
!TYPE(CoreInfo_Type) :: Core
!REAL :: s(nxy)
!INTEGER :: iz
!
!TYPE(Pin_Type), POINTER :: Pin(:)
!REAL :: Temp
!INTEGER :: ixy, idom, ineigh
!INTEGER :: i, j
!INTEGER :: ndom
!
!CALL CP_VA(Shat(1:nxy), s(1:nxy), nxy)
!Pin => Core%Pin
!ndom = RadDcmp%ndom
!DO idom = 1, ndom
!  DO i = 1, RadDcmp%nxylocal(idom)
!    ixy = RadDcmp%PinIdx(i, idom)
!    DO j = 1, 4
!      ineigh = Pin(ixy)%NeighIdx(j)
!      IF(ineigh .LT. 1) CYCLE
!      IF(RadDcmp%PinDom(ineigh) .EQ. idom) CYCLE
!      s(ixy) = s(ixy)  - shat(ineigh) / RadDB2(ineigh, iz)
!      !RadCoeffDat(j, ixy) = RadOffDiag1g(j, ixy)
!      !RadOffDiag1g(j, ixy) = 0
!    ENDDO
!  ENDDO
!ENDDO
!END SUBROUTINE


END MODULE

!SUBROUTINE SetRadDcmp()
!USE PARAM
!USE TYPEDEF,        ONLY : RadDcmp_TYPE
!IMPLICIT NONE
!END SUBROUTINE
