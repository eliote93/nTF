SUBROUTINE InitSubPlane()
USE PARAM
USE TYPEDEF
USE GEOM,   ONLY : nSubPlane,  nz,      nzfm,    SubPlaneMap,   SubPlaneRange, &
                   hz,         hzInv,   hzfm,    HzFmInv
!                   myzb,      myze, myzbf, myzef
USE CNTL,   ONLY : nTracerCntl                   
USE PE_MOD, ONLY : PE
USE ALLOCS
IMPLICIT NONE
INTEGER :: i
!INTEGER :: myzb, myze, myzbf, myzef
REAL :: rnSubPlane

nzfm = nSubPlane * nz
rnSubPlane = one/nSubPlane

!myzbf = nSubPlane * (myzb - 1) + 1
!myzef = nSubPlane * myze

CALL Dmalloc(SubPlaneMap, nzfm)
CALL Dmalloc(SubPlaneRange, 2, nz)
CALL Dmalloc(hzfm, nzfm)
CALL Dmalloc(hzfmInv, nzfm)

DO i = 1, nz
  SubPlaneRange(1, i) = nSubPlane *(i - 1) + 1
  SubPlaneRange(2, i) = nSubPlane * i
  SubPlaneMap(nSubPlane *(i - 1) + 1 : nSubPlane * i) = i
  SubPlaneMap(nSubPlane *(i - 1) + 1 : nSubPlane * i) = i
  hzfm(nSubPlane *(i - 1) + 1 : nSubPlane * i) = hz(i) * rnSubPlane
  hzfmInv(nSubPlane *(i - 1) + 1 : nSubPlane * i) = HzInv(i)
ENDDO

END SUBROUTINE