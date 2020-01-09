#include <defines.h>
SUBROUTINE VolCorChk(AsyGeom, RayInfo, CellInfo, PinUsed, nCellType, itype)
USE PARAM
USE TYPEDEF, ONLY : BasicGeom,          Cell_type,           RayInfo_Type,   &
                    AziAngleInfo_Type,  AsyRayInfo_type,     CellRayInfo_type
USE BasicOperation, ONLY : CP_CA, CP_VA, Multi_CA
USE ALLOCS
IMPLICIT NONE
TYPE(BasicGeom) :: AsyGeom
TYPE(RayInfo_type) :: RayInfo
TYPE(Cell_type), POINTER :: CellInfo(:)
LOGICAL, POINTER :: PinUsed(:, :)
INTEGER :: nCellType
INTEGER :: itype

TYPE(AziAngleInfo_Type), POINTER :: AziAng(:)
TYPE(AsyRayInfo_type), POINTER :: AsyRay(:)
TYPE(CellRayInfo_type), POINTER :: CellRay

INTEGER :: i, j, jbeg, jend, k, icel, iseg
INTEGER :: imray, iaray, iceray, iang, ipin, ireg
INTEGER :: nAsyRay, nAziAng, nxy
INTEGER :: nFsrMax, nFsr, nCellRay, nSeg
REAL, POINTER :: RTVol(:, :, :), Vol(:)
REAL :: len, del, maxv, minv,avg, nsample, noff

EQUIVALENCE(i, iang); EQUIVALENCE(j, iaray)
!EQUIVALENCE(k, iceray)

AziAng => RayInfo%AziAngle
AsyRay => RayInfo%AsyRay

nAziAng = RayInfo%nAziAngle
nAsyRay = RayInfo%nAsyRay
nxy = AsyGeom%nx*AsyGeom%ny
nFsrMax = 500

 ALLOCATE(RTVol(nFsrMax, nCellType, nxy))



DO i = 1, nAziAng  
  !Assembly Ray Sweeping
  maxv = 0._8; minv = 100._8; avg = 0.; nsample = 0.; noff = 0
  jbeg = AziAng(i)%AsyRayIdxSt; jend = jbeg + AziAng(i)%nAsyRay - 1
  del = AziAng(i)%del
  CALL CP_CA(RTVol, 0._8, nFsrMax, nCellType, nxy)
  RTVol = 0.
  DO j  = jbeg, jend
    imray = AsyRay(j)%ModRayIdx              !Corresponding Modular Ray 
    nCellRay = AsyRay(j)%nCellRay            !Number of CellRay
    IF(AsyRay(j)%PartialAsyRayFlag .NE. itype) CYCLE
    DO k = 1, nCellRay
      ipin = AsyRay(j)%PinIdx(k); iceray = AsyRay(j)%PinRayIdx(k)
      DO icel = 1, nCellType
        IF(CellInfo(icel)%BaseCellStr.NE.icel) CYCLE
        IF(.NOT. PinUsed(icel, ipin)) CYCLE
         CellRay => CellInfo(icel)%CellRay(iceray)
        nSeg = CellRay%nSeg
        DO iSeg = 1, nSeg
          iReg = CellRay%LocalFsrIdx(iSeg); len = CellRay%LenSeg(iSeg)
          RTVol(iReg, icel, ipin) = RTVol(iReg, icel, ipin) + len*Del*0.001_8
        ENDDO
      ENDDO !End of Cell Type Sweep
    ENDDO   !End of Cell Ray Sweeping
  ENDDO   !End of Assembly Ray Sweeping 
  !CALL Multi_CA(del, RTVol, nFsrMax, nCellType, nxy)
  
  !Calculate Correction
  DO ipin = 1, nxy   !Pin Sweep
    DO icel = 1, nCellType  !Cell Type Sweep
      IF(CellInfo(icel)%BaseCellStr.NE.icel) CYCLE
      IF(.NOT. PinUsed(icel, ipin)) CYCLE
      nFsr = CellInfo(icel)%nFsr
      Vol => CellInfo(icel)%vol
      DO ireg = 1, nFsr  !Flat Source Region Sweep
        RTVol(iReg, iCel, iPin) = Vol(iReg) / RTVol(iReg, iCel, iPin)
        IF( abs(1._8-RTVol(iReg, iCel, iPin)) .gt. 10E-4_8) noff = noff + 1
        !statistics
        maxv = max(RTVol(iReg, iCel, iPin), maxv)
        minv = min(RTVol(iReg, iCel, iPin), minv)
        nsample = nsample + 1
        avg = avg + RTVol(iReg, iCel, iPin)
        !IF(i .eq. 4) write(89, *)ipin, icel, ireg, RTVol(iReg, iCel, iPin)
      ENDDO  !Flat Source Region 
    ENDDO  !End of Cell Type Sweep
  ENDDO  !End of Pin Sweep
  avg = avg/ nsample
  WRITE(*, '(I5,A5,F10.6)')  i, 'MAX:', MAXv
  WRITE(*, '(5x,A5,F10.6)')    'MIN:', minv
  WRITE(*, '(5x,A5,F10.6)')    'AVG:', avg
  WRITE(*, '(5x,I15,a,I15, F10.6)'), noff, '/', nsample, (1-noff/nsample)*100
  WRITE(1, '(I5,A5,F10.6)')  i, 'MAX:', MAXv
  WRITE(1, '(5x,A5,F10.6)')    'MIN:', minv
  WRITE(1, '(5x,A5,F10.6)')    'AVG:', avg
  WRITE(1, '(5x,I15,a,I15, F10.6)'), noff, '/', nsample, (1-noff/nsample)*100
ENDDO  !End of Azimuthal Angle Sweeping
DEALLOCATE(RTVol)
NULLIFY(AziAng)
NULLIFY(AsyRay)
END SUBROUTINE

SUBROUTINE BaseCellRayVolCor(RayInfo, CellInfo, CellRayBase, nDiv)
!Make sets of Cellrays which 
USE PARAM
USE typedef, ONLY : Cell_type,  CellRayBase_Type,  RayInfo_Type,             &
                    CellRayInfo_type, AziAngleInfo_Type
!USE ALLOCS
USE BasicOperation, ONLY : CP_CA, MUlti_CA
USE Allocs
IMPLICIT NONE
TYPE(RayInfo_type) :: RayInfo
TYPE(Cell_type) :: CellInfo
TYPE(CellRayBase_Type) :: CellRayBase
INTEGER :: nDiv

TYPE(AziAngleInfo_Type), POINTER :: AziAngle(:)
TYPE(CellRayInfo_type), POINTER :: CellRay
INTEGER :: nAziAng, nseg
INTEGER :: iReg, nReg
INTEGER :: i, k, l
INTEGER :: j, jst, jbeg, jend
INTEGER, POINTER :: nlines(:)
REAL :: len, del
REAL :: VOl(10000)

nAziAng = RayInfo%nAziAngle
AziAngle => RayInfo%AziAngle
nReg= CellInfo%nFsr
nlines => CellRayBase%nlines

#define vcdbg  ! BYS edit 16/02/26 --- volcorrection debug on
#define Volcor1
DO i = 1, nAziAng
  del = AziAngle(i)%del
  
#ifdef Volcor1 
  DO jst = 0, nDiv-1
    
    CALL CP_CA(Vol(1:nReg), 0._8, nReg)
    jbeg = CellrayBase%AziAngIdx(i); jend = CellrayBase%AziAngIdx(i) + nlines(i) - 1
    
    DO j = jbeg + jst, jend, nDiv
      CellRay => CellRayBase%CellRay(j)
      nSeg = CellRay%nSeg
      DO k = 1, nSeg
        iReg = CellRay%LocalFsrIdx(k); len = CellRay%LenSeg(k)
        Vol(ireg) = Vol(ireg) + len
      ENDDO
    ENDDO
    
    CALL MULTI_CA(del, Vol(1:nReg), nReg)
    
    !Calculate Volume Correction Factors
    DO ireg = 1, nReg
 !     IF(Vol(iReg) .GT. epsm7) THEN
        Vol(iReg) = CellInfo%Vol(iReg) / Vol(iReg)
!      ELSE
!        Vol(iReg) = one
!      ENDIF
    ENDDO
    !Correct Cell Ray Lengths
    DO j = jbeg + jst, jend, nDiv
      CellRay => CellRayBase%CellRay(j)
      nSeg = CellRay%nSeg
      DO k = 1, nSeg
        iReg = CellRay%LocalFsrIdx(k); len = CellRay%LenSeg(k)
        !vol(iReg) = 1.
        len = len*vol(iReg); CellRay%LenSeg(k) = len * 1000._8
      ENDDO
    ENDDO    

#ifdef vcdbg 
    !Volume correction verification routines
    CALL CP_CA(Vol(1:nReg), 0._8, nReg)
    DO j = jbeg + jst, jend, nDiv
      CellRay => CellRayBase%CellRay(j)
      nSeg = CellRay%nSeg 
      DO k = 1, nSeg
        iReg = CellRay%LocalFsrIdx(k); len = CellRay%LenSeg(k)
        Vol(ireg) = Vol(ireg) + len*0.001_8
      ENDDO
    ENDDO
    
    CALL MULTI_CA(del, Vol(1:nReg), nReg)
    DO ireg = 1, nReg
        Vol(iReg) = CellInfo%Vol(iReg) / Vol(iReg)
        continue
    ENDDO    
    continue
#endif  
  ENDDO  
#else
   del = AziAngle(i)%del/nDIv
  CALL CP_CA(Vol(1:nReg), 0._8, nReg)
  jbeg = CellrayBase%AziAngIdx(i); jend = CellrayBase%AziAngIdx(i) + nlines(i) - 1
  
  DO j = jbeg , jend
    CellRay => CellRayBase%CellRay(j)
    nSeg = CellRay%nSeg
    DO k = 1, nSeg
      iReg = CellRay%LocalFsrIdx(k); len = CellRay%LenSeg(k)
      Vol(ireg) = Vol(ireg) + len
    ENDDO
  ENDDO
  
  CALL MULTI_CA(del, Vol(1:nReg), nReg)
  
  !Calculate Volume Correction Factors
  DO ireg = 1, nReg
!     IF(Vol(iReg) .GT. epsm7) THEN
      Vol(iReg) = CellInfo%Vol(iReg) / Vol(iReg)
      CONTINUE
!      ELSE
!        Vol(iReg) = one
!      ENDIF
  ENDDO
  !Correct Cell Ray Lengths
  DO j = jbeg, jend
    CellRay => CellRayBase%CellRay(j)
    nSeg = CellRay%nSeg
    DO k = 1, nSeg
      iReg = CellRay%LocalFsrIdx(k); len = CellRay%LenSeg(k)
      len = len*vol(iReg); CellRay%LenSeg(k) = len
    ENDDO
  ENDDO    
#ifdef vcdbg 
    !Volume correction verification routines
    CALL CP_CA(Vol(1:nReg), 0._8, nReg)
    DO j = jbeg , jend
      CellRay => CellRayBase%CellRay(j)
      nSeg = CellRay%nSeg 
      DO k = 1, nSeg
        iReg = CellRay%LocalFsrIdx(k); len = CellRay%LenSeg(k)
        Vol(ireg) = Vol(ireg) + len
      ENDDO
    ENDDO
    
    CALL MULTI_CA(del, Vol(1:nReg), nReg)
    DO ireg = 1, nReg
        Vol(iReg) = CellInfo%Vol(iReg) / Vol(iReg)
        continue
    ENDDO    
    continue
#endif  
#endif  
ENDDO

CALL DMALLOC(CellInfo%MaxSeg, nAziAng, nReg)
DO i = 1, nAziAng
  jbeg = CellrayBase%AziAngIdx(i); jend = CellrayBase%AziAngIdx(i) + nlines(i) - 1
  DO j = jbeg, jend
    CellRay => CellRayBase%CellRay(j)
    nSeg = CellRay%nSeg
    DO k = 1, nSeg
      iReg = CellRay%LocalFsrIdx(k); 
      CellInfo%MaxSeg(i, ireg) = max(CellRay%LenSeg(k), CellInfo%MaxSeg(i, ireg))
    ENDDO  
  ENDDO
  DO j = 1, nreg
    CellInfo%MaxSeg(i, j) = CellInfo%MaxSeg(i, j) * 0.001_8
  ENDDO
ENDDO

NULLIFY(AziAngle)
NULLIFY(CellRay)
NULLIFY(nlines)
END SUBROUTINE