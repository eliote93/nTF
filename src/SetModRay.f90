#define modray2
SUBROUTINE SetModularRayAngle(AziAngleInfo, RayInfo, AsyPitch)
USE PARAM
USE TYPEDEF, ONLY : AziAngleInfo_Type,   RayInfo_Type
!USE GEOM,    ONLY : AsyPitch
!USE RAYS,    ONLY : RayInfo,         AziAngle,        PolarAngle,              &
!                    ModRay          
IMPLICIT NONE

TYPE(AziAngleInfo_Type), POINTER :: AziANgleInfo(:)
TYPE(RayInfo_Type) :: RayInfo
REAL :: AsyPitch

INTEGER :: nAziAngle, nModRay
INTEGER :: i, j, k, iazi45
INTEGER :: nx, ny, nx0, ny0
REAL :: DelAng, Del, Del0, AziAngle, AziAngle0
REAL :: tanv, sinv, cosv, ang
REAL :: wt, wtsum
LOGICAL :: lEven

ALLOCATE(AziANgleInfo(RayInfo%nAziAngle))

nAziAngle = RayInfo%nAziAngle/2
Del0 = RayInfo%Del
DelAng = PI/RayInfo%nAziAngle
lEven = FALSE; iazi45 =0
IF(MOD(nAziAngle, 2) .EQ. 1) THEN
  lEven = TRUE; 
  iazi45 = INT(nAziAngle/2) + 1
ENDIF

DO i = 1, nAziAngle
  AziAngle0 = DelAng * (Real(i,8) - 0.5)

  !nx = int(AsyPitch * SIN(AziAngle0) / Del0)
  !ny = in(AsyPitch * COS(AziAngle0) / Del0)
  nx = AsyPitch * SIN(AziAngle0) / Del0
  ny = AsyPitch * COS(AziAngle0) / Del0
  nx = nx/2 + mod(nx,2); nx = 2 * nx
  ny = ny/2 + mod(ny,2); ny = 2* ny
#ifndef modray2    
  IF(mod(nx + ny, 2) .EQ. 0 .and. I .NE. iazi45) THEN
    IF(nx .GT. ny) THEN; nx = nx + 1
    ELSE; ny = ny + 1
    ENDIF
  ENDIF
#else
!  nx0 = nx; ny0 = ny
!  IF(mod(nx0, 2) .EQ. 1 .AND. mod(ny0, 2) .EQ. 1) THEN
!    nx = nx + 10; ny = ny + 10
!  ELSEIF(mod(nx0, 2) .EQ. 1) THEN
!    nx = nx + 9; ny = ny + 8
!  ELSEIF(mod(ny0, 2) .EQ. 1) THEN
!    ny = ny + 9; nx = nx + 8
!  ENDIF
!  IF(mod(nx, 2) .EQ. 1) nx = nx + 1
!  IF(mod(ny, 2) .EQ. 1) ny = ny + 1
!  IF(mod(nx + ny, 2) .EQ. 1 .and. I .NE. iazi45) THEN
!    IF(nx .GT. ny) THEN; nx = nx + 1
!    ELSE; ny = ny + 1
!    ENDIF
!  ENDIF
!CALL FindAngle(nx, ny, AziAngle0, del0, AsyPitch, 0.98_8)
#endif  
  tanv = (1.0_8*nx/ny); ang = atan(tanv)
  del = AsyPitch / sqrt(REAL(nx * nx + ny * ny,8))
#ifndef modray2     
  nModRay = 2 * NINT(REAL(HALF*(nx + ny), 8)) - 1
#else
  nModRay = (nx+ny)
#endif  
  AziAngleInfo(i)%ang = ang; AziAngleInfo(i)%tanv = tanv
  AziAngleInfo(i)%sinv = SIN(ang); AziAngleInfo(i)%cosv = COS(ang); 
  AziAngleInfo(i)%del = del; AziAngleInfo(i)%nModRay = nModRay
ENDDO

DO i = 1, nAziAngle
  j = 2*nAziAngle - i + 1
  AziAngleInfo(j)%ang = PI - AziAngleInfo(i)%ang
  Ang = AziAngleInfo(j)%ang
  AziAngleInfo(j)%Del = AziAngleInfo(i)%Del
  AziAngleInfo(j)%nModRay = AziAngleInfo(i)%nModRay
  AziAngleInfo(j)%sinv = SIN(Ang); AziAngleInfo(j)%cosv = COS(Ang)
  AziAngleInfo(j)%tanv = TAN(Ang)
ENDDO

!WEIGHT SETTING
wtsum = 1.0
AziAngle0 = 0.
DO i = 1, nAziAngle
  j = 2*nAziAngle - i + 1
  AziAngle = Half * (AziAngleInfo(i)%ang + AziAngleInfo(i+1)%ang)
  wt = (AziAngle - AziAngle0) / (2 * PI) * wtsum
  AziAngleInfo(i)%weight = wt; AziAngleInfo(j)%weight = wt;
  AziAngle0 =  AziAngle
ENDDO

wt = 0.
DO i = 1, 2 * nAziAngle
  wt = wt + AziAngleInfo(i)%weight
ENDDO
RayInfo%AziAngle => AziAngleInfo
END SUBROUTINE

SUBROUTINE FindAngle(nx, ny, ang0, del0, l, range)
IMPLICIT NONE
INTEGER :: nx, ny
REAL :: ang0, del0, l, range
REAL :: crit, chk, tanv, ang, errmin, err, del, crit_del
REAL :: del_err
INTEGER :: ix, iy

Crit = l/del0/range
Crit = Crit * Crit
crit_del = del0 * (1-range)
!tanv = (1.0_8*nx/ny); ang = atan(tanv)
ix = 0
errmin = 1.E+8_8
DO 
  ix = ix + 2; iy = 0
  chk = ix * ix
  IF(chk .GT. crit) EXIt
  DO 
    iy = iy + 2
    CHK = IX * IX + IY *IY
    IF(chk .GT. crit) EXIt
    tanv = (1.0_8*ix/iy); ang = atan(tanv)
    del = l / sqrt(REAL(ix * ix + iy * iy,8))
    del_err = abs(del-del0)
    IF(del_err .GT. crit_del) CYCLE
    err = abs(ang-ang0)
    IF(err .LT. errmin) THEN
      nx = ix; ny = iy
      errmin = err
    ENDIF
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE SetPolarRayAngle(PolarAngleInfo, RayInfo)
USE PARAM
USE TYPEDEF, ONLY : PolarAngle_Type,  RayInfo_Type
IMPLICIT NONE
TYPE(PolarAngle_Type), POINTER :: PolarAngleInfo(:)
TYPE(RayInfo_Type) :: RayInfo

INTEGER :: iPol, nPolAngle, nPolAnlgeHemi
INTEGER :: I, J, K
REAL :: DelAng, ang0, ang1
REAL :: wtsum
!#define newweight
#ifndef newweight 
REAL :: WeightSet(4, 3), SinvSet(4, 3)

DATA WeightSet /0.212854_8, 0.787146_8, 0._8,        0._8,        &
                0.046233_8, 0.283619_8, 0.670148_8, 0._8,        &
                0.0116656_8, 0.0855345_8, 0.3079014_8, 0.5948985_8/
                
DATA SinvSet  / 0.363900_8, 0.899900_8, 0._8,        0._8,        &
                0.166648_8, 0.537707_8, 0.932954_8, 0._8,        &
                0.0834435_8, 0.2908650_8, 0.6328485_8, 0.9492286_8/

!DATA WeightSet /0.1565149_8, 0.8434851_8, 0._8,        0._8,        &
!                0.0244934_8, 0.2184170_8, 0.7570896_8, 0._8,        &
!                0.0047107_8, 0.0448068_8, 0.2480350_8, 0.7024475_8/
!                
!DATA SinvSet  / 0.3012200_8, 0.8728932_8, 0._8,        0._8,        &
!                0.1187086_8, 0.4403554_8, 0.9057904_8, 0._8,        &
!                0.0521634_8, 0.1997629_8, 0.5164710_8, 0.9223376_8/
#else
REAL :: WeightSet(5, 4), SinvSet(5, 4)
DATA WeightSet / 0.218863_8, 0.781137_8, 0.000000_8, 0.000000_8, 0.000000_8,  &  
                 0.044355_8, 0.270249_8, 0.685397_8, 0.000000_8, 0.000000_8,  &  
                 0.010165_8, 0.070812_8, 0.289244_8, 0.629779_8, 0.000000_8,  &  
                 0.002782_8, 0.019393_8, 0.089318_8, 0.298508_8, 0.590000_8 /    


DATA SinvSet  /0.373917_8, 0.901601_8, 0.000000_8, 0.000000_8, 0.000000,  &
               0.166542_8, 0.523932_8, 0.928619_8, 0.000000_8, 0.000000,  &
               0.079818_8, 0.266013_8, 0.597389_8, 0.941545_8, 0.000000,  & 
               0.041631_8, 0.140115_8, 0.331811_8, 0.645183_8, 0.949724   /
#endif

#define UPPERHEMISWEEP      
          
#ifdef UPPERHEMISWEEP
  nPolAngle = RayInfo%nPolarAngle
  RayInfo%nPolarAngleHemi = nPolAngle 
#else
  RayInfo%nPolarAngleHemi = RayInfo%nPolarAngle
  RayInfo%nPolarAngle = RayInfo%nPolarAngle * 2
#endif

nPolAnlgeHemi = RayInfo%nPolarAngleHemi
nPolAngle = RayInfo%nPolarAngle

Allocate(PolarAngleInfo(nPolAngle))

ipol = nPolAnlgeHemi - 1
wtsum = 0.
#ifndef newweight
IF(nPolAnlgeHemi .LE. 4) THEN
#else
IF(nPolAnlgeHemi .LE. 5) THEN
#endif
  DO i = 1, nPolAnlgeHemi
    PolarAngleInfo(i)%weight = WeightSet(i, ipol)
    PolarAngleInfo(i)%sinv = SinvSet(i,ipol)
    PolarAngleInfo(i)%ang = ASIN(SinvSet(i,ipol))
    PolarAngleInfo(i)%cosv = COS(PolarAngleInfo(i)%ang)
  ENDDO  
ELSE
  DelAng = HALF * PI / nPolAnlgeHemi
  DO i = 1, nPolAnlgeHemi
    ang1 = i*DelAng; ang0 = ang1 - DelAng
    j = nPolAnlgeHemi + 1 - i
    PolarAngleInfo(j)%weight = SIN(ang1)- SIN(ang0)
    wtsum = wtsum + PolarAngleInfo(j)%weight
    PolarAngleInfo(j)%sinv = (half*(ang1-ang0)+0.25*(sin(2.*ang1)-sin(2.*ang0)))/PolarAngleInfo(j)%weight
    PolarAngleInfo(j)%ang = ASIN(PolarAngleInfo(j)%sinv)
    PolarANgleInfo(j)%cosv = COS(PolarAngleInfo(j)%ang)
  ENDDO
ENDIF

#ifndef UPPERHEMISWEEP
DO i = 1, nPolAnlgeHemi 
  j = nPolAngle - i + 1
  PolarAngleInfo(i)%weight =  Half * PolarAngleInfo(i)%weight
  PolarAngleInfo(j)%weight = PolarAngleInfo(i)%weight
  PolarAngleInfo(j)%ang = PI - PolarAngleInfo(i)%ang
  PolarANgleInfo(j)%sinv = SIN(PolarAngleInfo(i)%ang)
  PolarANgleInfo(j)%cosv = COS(PolarAngleInfo(i)%ang)
ENDDO
#endif

RayInfo%PolarAngle => PolarAngleInfo
END SUBROUTINE

SUBROUTINE SetModRay(ModRay, RayInfo, AsyGeom, AsyPitch, lEdge)
USE PARAM
USE TYPEDEF, ONLY : RayInfo_Type, ModRayInfo_Type, AziAngleInfo_Type,   &
                    BasicGeom
USE GeomTreatment, ONLY : GeomBoundary_RayIntersection
USE SetRay, ONLY : SetModRayLinking
USE UtilFunction                  
USE ioutil, ONLY : terminate                    
IMPLICIT NONE
TYPE(RayInfo_Type) :: RayInfo
TYPE(ModRayInfo_Type), POINTER :: ModRay(:)
TYPE(BasicGeom) :: AsyGeom(0:3)
REAL :: AsyPitch
LOGICAL :: lEdge

TYPE(AziAngleInfo_Type), POINTER :: AziAngle(:)
REAL :: AsyHalfPitch
REAL :: del, del_x, del_y, sinv, cosv, tanv
REAL :: RAY(3)
REAL :: x_intercept0, x_intercept
REAL :: InOutPt(2,2)
INTEGER :: InOutSurf(2)
INTEGER :: npt
INTEGER :: nMRay                     !Total Number of Moduar Ray
INTEGER :: nModRay                   !Number of Modular Ray in a certain angle
INTEGER :: nAziAngle, nAziAngle90
INTEGER :: nHalfModRay
INTEGER :: i, j ,k
INTEGER :: iazi45, iazi135
INTEGER :: imray                     !Modular Ray Index
INTEGER :: idxBeg, idxEnd
LOGICAL :: lEVEN


AziAngle => RayInfo%AziAngle
nAziAngle = RayInfo%nAziAngle
AsyHalfPitch = Half * AsyPitch

nmray = 0
DO i = 1, nAziAngle
  nmray = nmray + AziAngle(i)%nModRay
ENDDO
Allocate(ModRay(nmray))
RayInfo%nModRay = nmray

nAziAngle90 = nAziANgle/2
iazi45 = 0; iazi135 = 0; lEven = TRUE
IF(MOD(nAziAngle90, 2) .EQ. 1) THEN
  lEven = FALSE; 
  iazi45 = INT(nAziAngle90/2) + 1
  iazi135 = nAziAngle90 + iazi45
ENDIF

!Modular Ray In & Out Points
imray = 0
DO i = 1, nAziAngle
  nModRay = AziAngle(i)%nModRay
  nHalfModRay = nModRay/2 + mod(nModRay,2)
  AziAngle(i)%ModRayIdxSt = imray + 1
  del = AziAngle(i)%del
  sinv = AziAngle(i)%sinv; cosv = AziAngle(i)%cosv
  tanv = AziAngle(i)%tanv
  del_x = abs(AziAngle(i)%del / cosv)
#ifndef modray2    
  x_intercept0 = -del_x * nHalfModRay
#else
  x_intercept0 = -del_x * (nHalfModRay) -del_x * half 
#endif  
  DO j = 1, nModRay
#ifndef modray2      
    IF((iazi45 .EQ. i .OR. iazi135 .EQ. i) .AND. nHalfModRay .EQ. j) CYCLE !Skip Conner to Conner Ray
#endif    
    imray = imray + 1
    x_intercept = x_intercept0 + del_x * j
    RAY(1) = tanv; RAY(2) = -1.; RAY(3) = x_intercept
    !IF(i .NE. 6) CYCLE
    CALL GeomBoundary_RayIntersection(AsyGeom(0), InOutPt, InOutSurf, ray, npt)
    ModRay(imray)%InOutPoint = InOutPt;   ModRay(imray)%InOutSurf = InOutSurf
    ModRay(imray)%iAziAngIdx = i
    ModRay(imray)%NextRayIdx = 0
!    if(i .ne. 1) cycle
!    write(88, *) InOutPt(1,1), InOutPt(2,1)
!    write(88, *) InOutPt(1,2), InOutPt(2,2)
  ENDDO
#ifndef modray2   
  IF(iazi45 .EQ. i .OR. iazi135 .EQ. i) AziAngle(i)%nModRay = AziAngle(i)%nModRay - 1
#endif  
ENDDO

DO i = 1, nAziAngle
  nModRay = AziAngle(i)%nModRay
  idxBeg = AziAngle(i)%ModRayIdxSt
  idxEnd = AziAngle(i)%ModRayIdxSt + nModRay - 1
  CALL SetModRayLinking(Modray, idxBeg, idxEnd)
  !Checking 
  DO j = idxBeg, idxEnd
    IF(ModRay(j)%NextRayIdx(1) .EQ. 0) CALL terminate('SetMod Ray.f90: Can Not Make the modular Ray linking Info')
    IF(ModRay(j)%NextRayIdx(2) .EQ. 0) CALL terminate('SetMod Ray.f90: Can Not Make the modular Ray linking Info')
  ENDDO 
ENDDO
RayInfo%ModRay => ModRay
!stop
END SUBROUTINE

SUBROUTINE SetModRayLinking(ModRay, IdxBeg, IdxEnd)
USE PARAM
USE TYPEDEF, ONLY : ModRayInfo_Type
USE GeomTreatment, ONLY : xAxis_symPoint, yAxis_symPoint
USE ioutil, ONLY : terminate
IMPLICIT NONE
TYPE(ModRayInfo_Type), POINTER :: ModRay(:)
INTEGER :: IdxBeg, IdxEnd

INTEGER :: nModRay, imray
INTEGER :: isurf, jsurf, NextSurf(4) ,mp(2)           !SURFACE INDEX
INTEGER :: i, j, k, l
REAL :: x(2), xNext(2)

DATA NextSurf /3, 4, 1, 2/
DATA mp /2, 1/ 
!Modray Linking Information
!#define CompleteSearch
DO imray = idxBeg, idxEnd
  !Forward Direction : y Increase
#ifdef CompleteSearch    
  DO j = 1, 2
    isurf = ModRay(imray)%InOutSurf(j); jsurf = NextSurf(isurf)
    x = ModRay(imray)%InOutPoint(:, j)
    IF(isurf .EQ. WEST .OR. isurf .EQ. EAST) xNext = xAxis_symPoint(x)
    IF(isurf .EQ. NORTH .OR. isurf .EQ. SOUTH) xNext = yAxis_symPoint(x)
    !Linking POINT Search
    DO k = idxBeg, idxEnd
      l = k - idxBeg + 1
      IF(imray .EQ. k) CYCLE
      DO i = 1, 2
        IF(jsurf .NE. ModRay(k)%InOutSurf(i)) CYCLE
        IF(abs(xNext(1)-ModRay(k)%InOutPoint(1, i)) .GT. epsm4) CYCLE
        IF(abs(xNext(2)-ModRay(k)%InOutPoint(2, i)) .GT. epsm4) CYCLE
        goto 100
      ENDDO
    ENDDO

    IF(k .GT. idxEnd .or. i .GT. 2) CALL TERMINATE('SetMod Ray.f90: Can Not Make the modular Ray linking Info')
100 continue    
    ModRay(imray)%NextRayIdx(j) = k
  ENDDO
#else
  j = 2
  isurf = ModRay(imray)%InOutSurf(j); jsurf = NextSurf(isurf)
  x = ModRay(imray)%InOutPoint(:, j)
  IF(isurf .EQ. WEST .OR. isurf .EQ. EAST) xNext = xAxis_symPoint(x)
  IF(isurf .EQ. NORTH .OR. isurf .EQ. SOUTH) xNext = yAxis_symPoint(x)
  !Linking POINT Search
  DO k = idxBeg, idxEnd
    IF(imray .EQ. k) CYCLE
    i = 1
    IF(jsurf .NE. ModRay(k)%InOutSurf(i)) CYCLE
    IF(abs(xNext(1)-ModRay(k)%InOutPoint(1, i)) .GT. epsm4) CYCLE
    IF(abs(xNext(2)-ModRay(k)%InOutPoint(2, i)) .GT. epsm4) CYCLE
    goto 100
  ENDDO

  IF(k .GT. idxEnd .or. i .GT. 2) CALL TERMINATE('SetMod Ray.f90: Can Not Make the modular Ray linking Info')
100 continue    
  ModRay(imray)%NextRayIdx(j) = k
  ModRay(k)%NextRayIdx(mp(j)) = imray
#endif  
ENDDO
END SUBROUTINE

