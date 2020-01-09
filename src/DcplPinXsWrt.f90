SUBROUTINE DcplPinXsWrt(iunt, ipln, imod)
USE PARAM
USE TYPEDEF,     ONLY : CoreInfo_Type     ,Pin_Type
USE GEOM,        ONLY : Core              ,ng
IMPLICIT NONE

INTEGER :: iunt, ipln, imod

INTEGER :: nfn
CHARACTER(256) :: fn

TYPE(Pin_Type), POINTER :: Pin(:)

INTEGER :: ixy, itype
INTEGER :: nxy



ipln = 1; nfn = 7
WRITE(fn, '(A)') 'DcplPln'
IF(ipln .GT. 9) THEN
  WRITE(fn(nfn+1:nfn+12), '(I2)') ipln 
  nfn = nfn + 2 
ELSE
  WRITE(fn(nfn + 1:nfn + 1), '(I2)') ipln 
  nfn = nfn + 1  
ENDIF
WRITE(fn(nfn+1:nfn+5), '(A)') '_mode'
nfn= nfn + 5
WRITE(fn(nfn+1:nfn+1), '(A)') imod
nfn = nfn + 1
WRITE(fn(nfn+1:nfn+3), *) '.xs'
nfn = nfn + 3


nxy = Core%nxy
Pin => Core%Pin

OPEN(UNIT = iunt, FILE = fn, STATUS = 'REPLACE')
!GEOM INFO
!nxy, ng 
WRITE(iunt, '(20I10)') nxy, ng, imod
DO ixy = 1, nxy
  IF(Pin(ixy)%lfuel) itype = 1
  IF(Pin(ixy)%lGT) itype = 2
  IF(Pin(ixy)%lRadRef) itype = 3 
ENDDO
!WRITE(iunt, '(A10, 4x, )') 'PIN_MACRO'
!PIN_MACRO
CLOSE(iunt)

END SUBROUTINE

