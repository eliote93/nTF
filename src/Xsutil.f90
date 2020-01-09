#include <defines.h>
MODULE XsUtil_mod
USE PARAM
USE TYPEDEF,   ONLY : XsMac_Type
USE ALLOCS
IMPLICIT NONE

INTEGER, PARAMETER, PRIVATE :: ndatmax = 200
TYPE(XsMac_Type), POINTER, PRIVATE :: XsMacDat(:)
LOGICAL, PRIVATE :: luse(ndatmax) = .FALSE.
LOGICAL, PRIVATE :: lAlloc_XsMacDat = .FALSE.

INTERFACE LineIntPol
  MODULE PROCEDURE LineIntPolDouble
  MODULE PROCEDURE LineIntPolFloat
END INTERFACE

INTERFACE LineIntPol2
  MODULE PROCEDURE LineIntPol2Double
  MODULE PROCEDURE LineIntPol2Float
END INTERFACE

CONTAINS

SUBROUTINE GetXsMacDat(XsMac, ng, lIsoXsOut)
USE XSLIB_MOD,   ONLY : nelthel
IMPLICIT NONE
TYPE(XsMac_Type), POINTER :: XsMac
LOGICAL :: lIsoXsOut
INTEGER :: ng

INTEGER :: i
LOGICAL :: lget

!$OMP CRITICAL
IF(.NOT. lAlloc_XsMacDat) THEN
  ALLOCATE(XsMacDat(nDatMax))
  lAlloc_XsMacDat = .TRUE.
ENDIF
lget = .FALSE.
DO i = 1, nDatMax
  IF(lUse(i)) CYCLE
  XsMacDat(i)%id = i;
  lUse(i) = .TRUE.
  XsMac => XsMacDat(i)
  IF(XsMac%lalloc .AND. XsMac%ng .NE. ng) THEN
    CALL FreeXsMac(XsMac)
    IF(XsMac%lisoalloc) CALL FreeXsIsoMac(xsMac)
  ENDIF
  IF(.NOT. XsMac%lAlloc) THEN
    XsMac%ng = ng; CALL AllocXsMac(XsMac)
  ENDIF
  IF(lIsoXsOut .AND. .NOT. XsMac%lIsoAlloc) CALL AllocMacIsoXs(XsMac, ng, nelthel)
  EXIT
ENDDO
!$OMP END CRITICAL
END SUBROUTINE

SUBROUTINE ReturnXsMacDat(XsMac)
IMPLICIT NONE
TYPE(XsMac_Type), POINTER :: XsMac
!$OMP CRITICAL
lUse(XsMac%id) = .FALSE.
NULLIFY(XsMac)
!$OMP END CRITICAL
END SUBROUTINE

SUBROUTINE AllocXsMac(XsMac)
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac
INTEGER :: ig, ng
ng = XsMac%ng
CALL Dmalloc(XsMac%xsmact, ng);  CALL Dmalloc(XsMac%xsmactr, ng)
CALL Dmalloc(XsMac%xsmaca, ng);  CALL Dmalloc(XsMac%xsmacf, ng)
CALL Dmalloc(XsMac%xsmacnf, ng); CALL Dmalloc(XsMac%xsmackf, ng)
CALL Dmalloc(XsMac%xsmacsm, ng, ng); 
CALL Dmalloc(XsMac%xsmacp1sm, ng, ng); CALL Dmalloc(XsMac%xsmacp2sm, ng, ng); CALL Dmalloc(XsMac%xsmacp3sm, ng, ng)
XsMac%lAllocSM=TRUE
CALL Dmalloc(XsMac%xsmacs, ng)
CALL Dmalloc(XsMac%xsmacstr, ng); CALL Dmalloc(XsMac%Chi, ng)
XsMac%lalloc = TRUE
END SUBROUTINE

SUBROUTINE FreeXsMac(XsMac)
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac
IF(.NOT. ASSOCIATED(XsMac%xsmact)) RETURN
Deallocate(XsMac%xsmact);  Deallocate(XsMac%xsmactr)
Deallocate(XsMac%xsmaca);  Deallocate(XsMac%xsmacnf)
Deallocate(XsMac%xsmackf);  Deallocate(XsMac%xsmacf) 
Deallocate(XsMac%xsmacsm); Deallocate(XsMac%xsmacs)
Deallocate(XsMac%xsmacstr); Deallocate(XsMac%Chi)
IF(XsMac%lAllocSm)THEN
    DEALLOCATE(XsMac%xsmacp1sm)
    DEALLOCATE(XsMac%xsmacp2sm)
    DEALLOCATE(XsMac%xsmacp3sm)
    XsMac%lAllocSm=.FALSE.
ENDIF

!IF(ALLOCATED(XsMac%xsmacp1sm)) DEALLOCATE(XsMac%xsmacp1sm)
!IF(ALLOCATED(XsMac%xsmacp2sm)) DEALLOCATE(XsMac%xsmacp2sm)
!IF(ALLOCATED(XsMac%xsmacp3sm)) DEALLOCATE(XsMac%xsmacp3sm)
XsMac%lAlloc = .FALSE.
XsMac%ng = 0
END SUBROUTINE

SUBROUTINE AllocMacIsoXs(XsIsoMac, ng, niso)
IMPLICIT NONE
INTEGER :: ng, niso
TYPE(XsMac_Type) :: XsIsoMac
CALL Dmalloc(XsIsoMac%IsoXsMacA, niso, ng)
CALL Dmalloc(XsIsoMac%IsoXsMacnf, niso, ng)
CALL Dmalloc(XsIsoMac%IsoXsMacf, niso, ng)
CALL Dmalloc(XsIsoMac%IsoXsMackf, niso, ng)
CALL Dmalloc(XsIsoMac%IsoXsMacS0, niso, ng)
CALL Dmalloc(XsIsoMac%IsoXsMacS1, niso, ng)
CALL Dmalloc(XsIsoMac%IsoXsMacSS, niso, ng)
CALL Dmalloc(XsIsoMac%IsoXsMacTr, niso, ng)
CALL Dmalloc(XsIsoMac%IsoXsMacT, niso, ng) !BYS edit
CALL Dmalloc(XsIsoMac%IsoXsRadCap, niso, ng)
XsIsoMac%lIsoalloc = .TRUE.
XsIsoMac%niso = niso
END SUBROUTINE

SUBROUTINE FreeXsIsoMac(XsIsoMac)
TYPE(XsMac_Type) :: XsIsoMac
IF(.NOT. XsIsoMac%lIsoAlloc) RETURN
DEALLOCATE(XsIsoMac%IsoXsMacA);  DEALLOCATE(XsIsoMac%IsoXsMacnf)  
DEALLOCATE(XsIsoMac%IsoXsMacf);  DEALLOCATE(XsIsoMac%IsoXsMackf)  
DEALLOCATE(XsIsoMac%IsoXsMacS0); DEALLOCATE(XsIsoMac%IsoXsMacTr)  
DEALLOCATE(XsIsoMac%IsoXsRadCap)
DEALLOCATE(XsIsoMac%IsoXsMacS1); DEALLOCATE(XsIsoMac%IsoXsMacSS); 
DEALLOCATE(XsIsoMac%IsoXsMacT); 
XsIsoMac%lIsoAlloc = .FALSE.
XsIsoMac%niso = 0
END SUBROUTINE

SUBROUTINE XsTempInterpolation(id, isodata, temp, wt1, wt2, it1, it2)
USE PARAM
USE XSLIB_MOD,       ONLY : LIBDATA,itempmap
IMPLICIT NONE
TYPE(LIBDATA), INTENT(IN) :: isodata
INTEGER, INTENT(IN) :: id
REAL, INTENT(IN) :: TEMP
REAL :: wt1, wt2
INTEGER :: it1, it2

it1 = temp; it1 = itempmap(it1, id)
it2 = it1;   wt2 = 1
IF(isodata%ntemp .ne. 1 .and. it1 .lt. isodata%ntemp) THEN
  it2 = it1 + 1
  wt2 = (temp - isodata%temp(it1))/(isodata%temp(it2)-isodata%temp(it1))
ENDIF
wt1 = 1._8 - wt2

END SUBROUTINE

SUBROUTINE P1XsTempInterpolation(id, isodata, temp, wt1, wt2, it1, it2)
USE PARAM
USE XSLIB_MOD,       ONLY : LIBDATA,itempmapp1
IMPLICIT NONE
TYPE(LIBDATA), INTENT(IN) :: isodata
INTEGER, INTENT(IN) :: id
REAL, INTENT(IN) :: TEMP
REAL :: wt1, wt2
INTEGER :: it1, it2

it1 = temp; it1 = itempmapp1(it1, id)
it2 = it1;   wt2 = 1
IF(isodata%np1temp .ne. 1 .and. it1 .lt. isodata%np1temp) THEN
  it2 = it1 + 1
  wt2 = (temp - isodata%p1temp(it1))/(isodata%p1temp(it2)-isodata%p1temp(it1))
ENDIF
wt1 = 1._8 - wt2

END SUBROUTINE

SUBROUTINE P2XsTempInterpolation(id, isodata, temp, wt1, wt2, it1, it2)
USE PARAM
USE XSLIB_MOD,       ONLY : LIBDATA,itempmapp1
IMPLICIT NONE
TYPE(LIBDATA), INTENT(IN) :: isodata
INTEGER, INTENT(IN) :: id
REAL, INTENT(IN) :: TEMP
REAL :: wt1, wt2
INTEGER :: it1, it2

it1 = temp; it1 = itempmapp1(it1, id)
it2 = it1;   wt2 = 1
IF(isodata%np1temp .ne. 1 .and. it1 .lt. isodata%np1temp) THEN
  it2 = it1 + 1
  wt2 = (temp - isodata%p1temp(it1))/(isodata%p1temp(it2)-isodata%p1temp(it1))
ENDIF
wt1 = 1._8 - wt2

END SUBROUTINE

SUBROUTINE P3XsTempInterpolation(id, isodata, temp, wt1, wt2, it1, it2)
USE PARAM
USE XSLIB_MOD,       ONLY : LIBDATA,itempmapp1
IMPLICIT NONE
TYPE(LIBDATA), INTENT(IN) :: isodata
INTEGER, INTENT(IN) :: id
REAL, INTENT(IN) :: TEMP
REAL :: wt1, wt2
INTEGER :: it1, it2

it1 = temp; it1 = itempmapp1(it1, id)
it2 = it1;   wt2 = 1
IF(isodata%np1temp .ne. 1 .and. it1 .lt. isodata%np1temp) THEN
  it2 = it1 + 1
  wt2 = (temp - isodata%p1temp(it1))/(isodata%p1temp(it2)-isodata%p1temp(it1))
ENDIF
wt1 = 1._8 - wt2

END SUBROUTINE

FUNCTION LineIntPolDouble(x, ndat, xdat, ydat)
USE PARAM
IMPLICIT NONE
REAL :: x, LineIntPolDouble
INTEGER :: ndat
REAL :: xdat(ndat), ydat(ndat)
REAL :: wt1, wt2
INTEGER :: i, n1, n2

IF(ndat .EQ. 1) THEN
  LineIntPolDouble = ydat(1); RETURN
ENDIF

DO i = 2, ndat
  if(x .le. xdat(i)) EXIT
ENDDO
IF(i .GT. ndat) i = ndat
n2 = i; n1 = n2 - 1
wt2 = (x - xdat(n1))/(xdat(n2) - xdat(n1))
wt1 = 1 - wt2
LineIntPolDouble = wt2 * ydat(n2) + wt1 * ydat(n1)
END FUNCTION

FUNCTION LineIntPol2Double(x, ndat, xdat, ydat)
USE PARAM
IMPLICIT NONE
REAL :: x, LineIntPol2Double
INTEGER :: ndat
REAL :: xdat(ndat), ydat(ndat)
REAL :: wt1, wt2
INTEGER :: i, n1, n2

IF(ndat .EQ. 1) THEN
  LineIntPol2Double = ydat(1); RETURN
ENDIF

IF (x.le.xdat(1)) THEN
    LineIntPol2Double = ydat(1); RETURN
ELSEIF (x.ge.xdat(ndat)) THEN
    LineIntPol2Double = ydat(ndat); RETURN
ENDIF

DO i = 2, ndat
  if(x .le. xdat(i)) EXIT
ENDDO
IF(i .GT. ndat) i = ndat
n2 = i; n1 = n2 - 1
wt2 = (x - xdat(n1))/(xdat(n2) - xdat(n1))
wt1 = 1 - wt2
LineIntPol2Double = wt2 * ydat(n2) + wt1 * ydat(n1)
END FUNCTION

FUNCTION LineIntPolFloat(x, ndat, xdat, ydat)
USE PARAM
IMPLICIT NONE
REAL :: x, LineIntPolFloat
INTEGER :: ndat
REAL :: xdat(ndat)
REAL(4) :: ydat(ndat)
REAL :: wt1, wt2
INTEGER :: i, n1, n2

IF(ndat .EQ. 1) THEN
  LineIntPolFloat = ydat(1); RETURN
ENDIF

DO i = 2, ndat
  if(x .le. xdat(i)) EXIT
ENDDO
IF(i .GT. ndat) i = ndat
n2 = i; n1 = n2 - 1
wt2 = (x - xdat(n1))/(xdat(n2) - xdat(n1))
wt1 = 1 - wt2
LineIntPolFloat = wt2 * ydat(n2) + wt1 * ydat(n1)
END FUNCTION

FUNCTION LineIntPol2Float(x, ndat, xdat, ydat)
USE PARAM
IMPLICIT NONE
REAL :: x, LineIntPol2Float
INTEGER :: ndat
REAL :: xdat(ndat)
REAL(4) :: ydat(ndat)
REAL :: wt1, wt2
INTEGER :: i, n1, n2

IF(ndat .EQ. 1) THEN
  LineIntPol2Float = ydat(1); RETURN
ENDIF

IF (x.le.xdat(1)) THEN
    LineIntPol2Float = ydat(1); RETURN
ELSEIF (x.ge.xdat(ndat)) THEN
    LineIntPol2Float = ydat(ndat); RETURN
ENDIF

DO i = 2, ndat
  if(x .le. xdat(i)) EXIT
ENDDO
IF(i .GT. ndat) i = ndat
n2 = i; n1 = n2 - 1
wt2 = (x - xdat(n1))/(xdat(n2) - xdat(n1))
wt1 = 1 - wt2
LineIntPol2Float = wt2 * ydat(n2) + wt1 * ydat(n1)
END FUNCTION

SUBROUTINE calcWgt(x,arrX,N,w,idx,flag)
    INTEGER,INTENT(IN) :: N
    CHARACTER*1,INTENT(IN) :: flag
    REAL,INTENT(IN) :: arrX(N)
    REAL,INTENT(IN) :: x
    REAL,INTENT(OUT) :: w(2)
    INTEGER,INTENT(OUT) :: idx(2)
    INTEGER :: i,n1,n2
    IF (N.eq.1) THEN
        w(1)=1._8; w(2)=0._8
        RETURN
    ENDIF
    if (flag.eq.'A') then
        DO i=2,N
            IF (x.le.arrX(i)) EXIT
        ENDDO
        IF (i.gt.N) i=N
        n2=i; n1=n2-1
    elseif (flag.eq.'D') then
        DO i=N-1,1,-1
            IF (x.ge.arrX(i)) EXIT
        ENDDO
        IF (i.lt.1) i=1
        n1=i; n2=n1+1
    endif
    idx(1)=n1
    idx(2)=n2
    w(2) = (x - arrX(n1))/(arrX(n2) - arrX(n1))
    w(1) = 1._8 - w(2)
END SUBROUTINE

SUBROUTINE calcWgt2(x,arrX,N,w,idx,flag)
    INTEGER,INTENT(IN) :: N
    CHARACTER*1,INTENT(IN) :: flag
    REAL,INTENT(IN) :: arrX(N)
    REAL,INTENT(IN) :: x
    REAL,INTENT(OUT) :: w(2)
    INTEGER,INTENT(OUT) :: idx(2)
    INTEGER :: i,n1,n2
    IF (N.eq.1) THEN
        w(1)=1._8; w(2)=0._8
        RETURN
    ENDIF
    if (flag.eq.'A') then
        DO i=2,N
            IF (x.le.arrX(i)) EXIT
        ENDDO
        IF (i.gt.N) i=N
        n2=i; n1=n2-1
    elseif (flag.eq.'D') then
        IF (x.ge.arrX(N)) then
            idx(1)=N-1; idx(2)=N;
            w(1)=0._8; w(2)=1._8;
            RETURN
        ELSEIF (x.lt.arrX(1)) then
            idx(1)=1; idx(2)=2;
            w(1)=1._8; w(2)=0._8;
            RETURN
        ENDIF
        DO i=N-1,1,-1
            IF (x.ge.arrX(i)) EXIT
        ENDDO
        IF (i.lt.1) i=1
        n1=i; n2=n1+1
    endif
    idx(1)=n1
    idx(2)=n2
    w(2) = (x - arrX(n1))/(arrX(n2) - arrX(n1))
    w(1) = 1._8 - w(2)
END SUBROUTINE

SUBROUTINE SetXeDynEnv(idiso, pnum, niso, niso_depl, ntiso)
USE PARAM
USE BasicOperation
IMPLICIT NONE
INTEGER :: IdIso(ntiso)
REAL :: pnum(ntiso)
INTEGER, INTENT(INOUT) :: niso, niso_depl
INTEGER, INTENT(IN) :: ntiso

REAL :: pnum0(1000)
INTEGER :: idiso0(1000)
INTEGER :: niso0, niso_Depl0

LOGICAL :: lXeExist, lIExist
REAL :: pnumXE, pnumI
INTEGER :: i
lXeExist = .FALSE.
DO i = 1, niso
  IF(idiso(i) .EQ. 54635) THEN
    pnumXE = pnum(i)
    lXeExist = .TRUE.
    EXIT
  ENDIF
ENDDO
lIExist = .FALSE.
DO i = 1, niso
  IF(idiso(i) .EQ. 53635) THEN
    pnumI = pnum(i)
    lIExist = .TRUE.
    EXIT
  ENDIF
ENDDO
IF(lXeExist .AND. lIExist) RETURN
niso0 = niso; niso_depl0 = niso_depl
IF(.NOT. lXeExist) THEN
  DO i = niso + 1, niso_depl
     IF(idiso(i) .EQ. 54635) THEN
      pnumXE = pnum(i)
      EXIT
    ENDIF 
  ENDDO
  pnumXE= max(pnumXE, epsm30)
ENDIF

IF(.NOT. lIExist) THEN
  DO i = niso + 1, niso_depl
     IF(idiso(i) .EQ. 53635) THEN
      pnumI = pnum(i)
      EXIT
    ENDIF 
  ENDDO
  pnumI= max(pnumI, epsm30)
ENDIF

CALL CP_VA(idiso0(1:ntiso), idiso(1:ntiso), ntiso)
CALL CP_VA(pnum0(1:ntiso), pnum(1:ntiso), ntiso)
IF(.NOT. lXeExist) THEN
  niso = niso + 1
  IdIso(niso) = 54635; pnum(niso) = pnumXE
ENDIF
IF(.NOT. lIExist) THEN
  niso = niso + 1
  IdIso(niso) = 53635; pnum(niso) = pnumI
ENDIF
niso_depl = niso
DO i = niso0+1, niso_depl0
  IF(.NOT. lXeExist .AND. idiso0(i) .EQ. 54635) CYCLE
  IF(.NOT. lIExist .AND. idiso0(i) .EQ. 53635) CYCLE
  niso_depl = niso_depl + 1
  IdIso(niso_depl) = idiso0(i)
  pnum(niso_depl) = pnum0(i)
ENDDO
END SUBROUTINE

SUBROUTINE CnvrtNID(nid)
  IMPLICIT NONE
  INTEGER :: nid,xxx
  xxx = mod(nid,1000)
  if (xxx.gt.500) nid = nid - 500
END SUBROUTINE

END MODULE