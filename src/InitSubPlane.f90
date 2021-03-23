SUBROUTINE InitSubPlane()
USE PARAM
USE TYPEDEF
USE GEOM,   ONLY : nSubPlane,  nz,      nzfm,    SubPlaneMap,   SubPlaneRange, &
                   hz,         hzInv,   hzfm,    HzFmInv, maxhzfm, nSubPlanePtr
USE CNTL,   ONLY : nTracerCntl                   
USE ALLOCS
IMPLICIT NONE
INTEGER :: i
REAL :: inv_hzfm
REAL :: rnSubPlane

IF(maxhzfm .GT. 0.) THEN

  inv_hzfm = one / maxhzfm
  CALL Dmalloc(nSubPlanePtr, nz)
  CALL Dmalloc(SubPlaneRange, 2, nz)

  nzfm = 0
  DO i = 1, nz
    SubPlaneRange(1, i) = nzfm + 1
    nSubPlanePtr(i) = CEILING(hz(i) * inv_hzfm)
    nzfm = nzfm + nSubPlanePtr(i)
    SubPlaneRange(2, i) = nzfm
  END DO

  CALL Dmalloc(SubPlaneMap, nzfm)
  CALL Dmalloc(hzfm, nzfm)
  CALL Dmalloc(hzfmInv, nzfm)
  DO i = 1, nz
    rnSubPlane = one / nSubPlanePtr(i)
    SubPlaneMap(SubPlaneRange(1, i) : SubPlaneRange(2, i)) = i
    hzfm(SubPlaneRange(1, i) : SubPlaneRange(2, i)) = hz(i) * rnSubPlane
    hzfmInv(SubPlaneRange(1, i) : SubPlaneRange(2, i)) = HzInv(i) * nSubPlanePtr(i)
  END DO

ELSE
  nzfm = nSubPlane * nz
  rnSubPlane = one/nSubPlane

  CALL Dmalloc(nSubPlanePtr, nz)
  CALL Dmalloc(SubPlaneMap, nzfm)
  CALL Dmalloc(SubPlaneRange, 2, nz)
  CALL Dmalloc(hzfm, nzfm)
  CALL Dmalloc(hzfmInv, nzfm)

  DO i = 1, nz
    nSubPlanePtr(i) = nSubPlane
    SubPlaneRange(1, i) = nSubPlane *(i - 1) + 1
    SubPlaneRange(2, i) = nSubPlane * i
    SubPlaneMap(nSubPlane *(i - 1) + 1 : nSubPlane * i) = i
    SubPlaneMap(nSubPlane *(i - 1) + 1 : nSubPlane * i) = i
    hzfm(nSubPlane *(i - 1) + 1 : nSubPlane * i) = hz(i) * rnSubPlane
    hzfmInv(nSubPlane *(i - 1) + 1 : nSubPlane * i) = HzInv(i) * nSubPlane
  ENDDO
END IF

END SUBROUTINE
  
SUBROUTINE InitSubPlane_OLD()
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