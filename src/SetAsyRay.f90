SUBROUTINE SetAsyRay(AsyRay, RayInfo, AsyGeom, lEdge)
USE PARAM
USE TYPEDEF, ONLY : RayInfo_Type,       AsyRayInfo_type,   ModRayInfo_type, & 
                    AziAngleInfo_Type,  BasicGeom
USE GeomTreatment, ONLY : GeomBoundary_RayIntersection
IMPLICIT NONE
TYPE(RayInfo_Type) :: RayInfo
TYPE(AsyRayInfo_type), POINTER :: AsyRay(:)
TYPE(BasicGeom) :: AsyGeom(0:3)
LOGICAL :: lEdge

TYPE(ModRayInfo_type), POINTER :: ModRay(:)
TYPE(AziAngleInfo_Type),POINTER :: AziAngle(:)

!RAY INDEX
INTEGER :: nAsyRay, nModRay, nAziAng
INTEGER :: iasyray, imray
INTEGER :: i, j, k, jbeg, jend

!Hal Asy Ray Variables
REAL :: x(2), y(2), xnew(2), ynew(2), xlim(2), ylim(2)
REAL :: InOutPts(2,2),RAY(3)
REAL :: tanv, len
INTEGER :: itype, npt, InOutSurf(2)

AziAngle => RayInfo%AziAngle
ModRay => RayInfo%ModRay
nModRay = RayInfo%nModRay
nAziAng = RayInfo%nAziAngle
IF(lEdge) THEN
  nAsyRay = RayInfo%nModRay
ELSE
  nAsyRay = RayInfo%nModRay*4
ENDIF

ALLOCATE(ASyRay(nAsyRay))

IF(lEdge) Then
  DO i = 1, nAziAng
    nModRay = AziAngle(i)%nModRay
    jbeg = AziAngle(i)%ModRayIdxSt; jend = jbeg + AziAngle(i)%nModRay - 1
    AziAngle(i)%AsyRayIdxSt = jbeg
    AziAngle(i)%nAsyRay = nModray
    DO j = jbeg, jend
      iasyray = j  
      AsyRay(iasyray)%AziAngIdx = i; AsyRay(iasyray)%ModRayidx = j
      AsyRay(iasyray)%NextRayIdx = ModRay(iasyray)%NextRayIdx
      AsyRay(iasyray)%InOutSurf = ModRay(iasyray)%InOutSurf
      AsyRay(iasyray)%InOutPoint = ModRay(iasyray)%InOutPoint
      AsyRay(iasyray)%PartialAsyRayIdx = 0
    ENDDO
  ENDDO
  RayInfo%AsyRay => AsyRay 
  RayInfo%nAsyRay = nAsyRay
  RETURN
ENDIF

nAsyRay = 0
DO i = 1, nAziAng
  nModRay = AziAngle(i)%nModRay
  jbeg = AziAngle(i)%ModRayIdxSt; jend = jbeg + AziAngle(i)%nModRay - 1
  AziAngle(i)%AsyRayIdxSt = nAsyRay + 1
  DO j = jbeg, jend
    nAsyRay = nAsyRay + 1; iasyray = nAsyRay
    AsyRay(iasyray)%AziAngIdx = i; AsyRay(iasyray)%ModRayidx = j
    !AsyRay(iasyray)%NextRayIdx = ModRay(j)%NextRayIdx
    AsyRay(iasyray)%InOutSurf = ModRay(j)%InOutSurf
    AsyRay(iasyray)%InOutPoint = ModRay(j)%InOutPoint
    AsyRay(iasyray)%PartialAsyRayIdx = 0
    AsyRay(iasyray)%PartialAsyRayIdx(0) = j
    AsyRay(iasyray)%NextRayIdx(1) = (ModRay(j)%NextRayIdx(1)-AziAngle(i)%ModRayIdxSt + 1) + AziAngle(i)%AsyRayIdxSt - 1
    AsyRay(iasyray)%NextRayIdx(2) = (ModRay(j)%NextRayIdx(2)-AziAngle(i)%ModRayIdxSt + 1) + AziAngle(i)%AsyRayIdxSt - 1
  ENDDO
  
  jbeg = AziAngle(i)%AsyRayIdxSt
  jend = jbeg + nModRay - 1
  tanv = AziAngle(i)%tanv
  DO itype = 1,3
    xlim = AsyGeom(itype)%x; ylim = AsyGeom(itype)%y
    DO j = jbeg, jend
      imray = AsyRay(j)%ModRayIdx
      x = ModRay(imray)%InOutPoint(1, :); y = ModRay(imray)%InOutPoint(2, :)
      
!      IF(tanv .GT. ZERO) THEN  !Search the modular rays which are neccesarry for setting partial rays
!        IF( (x(1)-xlim(1)) * (x(1)-xlim(2)) .GT. 0.) CYCLE
!        IF( (y(1)-ylim(1)) * (y(1)-ylim(2)) .GT. 0.) CYCLE
!      ELSE
!        IF( (x(2)-xlim(1)) * (x(2)-xlim(2)) .GT. 0.) CYCLE
!        IF( (y(2)-ylim(1)) * (y(2)-ylim(2)) .GT. 0.) CYCLE    
!      ENDIF
!      
      RAY(1) = tanv; RAY(2) = -1.; RAY(3) = -tanv * x(1) + y(1)
      CALL GeomBoundary_RayIntersection(AsyGeom(itype), InOutPts, InOutSurf, ray, npt)
      IF(npt .LE. 1) cycle      
      len =  (InOutPts(1, 1)-InOutPts(1, 2))**2 + &
                (InOutPts(2, 1)-InOutPts(2, 2))**2
      len = sqrt(len)      
      IF(len .lt. epsm8) CYCLE   !Exclude the zero length ray

      nAsyRay = nAsyRay + 1; iasyray = nAsyRay
      AsyRay(nAsyRay)%PartialAsyRayFlag = itype
      AsyRay(j)%PartialAsyRayIdx(itype) = iasyray
      AsyRay(nAsyRay)%AziAngIdx = i; AsyRay(nAsyRay)%ModRayIdx = imray
      AsyRay(iasyray)%InOutPoint = InOutPts;   AsyRay(iasyray)%InOutSurf = InOutSurf
      AsyRay(iasyray)%PartialAsyRayFlag = itype;
      AsyRay(iasyray)%PartialAsyRayIdx(0) = j; AsyRay(iasyray)%PartialAsyRayIdx(itype) = iasyray
      AsyRay(iasyRay)%NextRayIdx = AsyRay(j)%NextRayIdx
#ifdef AsyRayOut
      IF(i .EQ. 3 .and. itype .eq. 2) THEN
        WRITE(88, '(2e20.6)') AsyRay(iasyray)%InOutPoint(1, 1), AsyRay(iasyray)%InOutPoint(2, 1)
        WRITE(88, '(2e20.6)') AsyRay(iasyray)%InOutPoint(1, 2), AsyRay(iasyray)%InOutPoint(2, 2)
      ENDIF     
#endif       
    ENDDO    
  ENDDO
#ifdef DEBUG_AsyRay  
  itype = 1
  DO j = jbeg, jend
    imray = AsyRay(j)%ModRayIdx
    x = ModRay(imray)%InOutPoint(1, :); y = ModRay(imray)%InOutPoint(2, :)
    IF((y(1)*y(2)) .GT. epsm5) CYCLE
    nAsyRay = nAsyRay + 1; iasyray = nAsyRay
    AsyRay(iasyray)%PartialAsyRayFlag = itype;
    !Partail Ray <=> normal asy ray mapping
    AsyRay(iasyray)%PartialAsyRayIdx(0) = j; AsyRay(iasyray)%PartialAsyRayIdx(itype) = iasyray
    AsyRay(j)%PartialAsyRayIdx(itype) = iasyray
    !Finding New in-out point
    AsyRay(iasyray)%InOutSurf = AsyRay(j)%InOutSurf
    AsyRay(iasyRay)%InOutPoint = AsyRay(j)%InOutPoint
    AsyRay(iasyRay)%NextRayIdx = AsyRay(j)%NextRayIdx
    
    AsyRay(iasyray)%InOutSurf(2) = North
    AsyRay(iasyray)%InOutPoint(1, 2) = x(1) - y(1)/tanv
    AsyRay(iasyray)%InOutPoint(2, 2) = 0.
    
    len =  (AsyRay(iasyray)%InOutPoint(1, 1)-AsyRay(iasyray)%InOutPoint(1, 2))**2 + &
            (AsyRay(iasyray)%InOutPoint(2, 1)-AsyRay(iasyray)%InOutPoint(2, 2))**2
    len = sqrt(len)
    IF(len .lt. epsm5) THEN
      nAsyRay = nAsyRay - 1
      AsyRay(j)%PartialAsyRayIdx(itype) = 0
      CYCLE
    ENDIF
!    IF(i .EQ. 9) THEN
!      WRITE(88, '(2e20.6)') AsyRay(iasyray)%InOutPoint(1, 1), AsyRay(iasyray)%InOutPoint(2, 1)
!      WRITE(88, '(2e20.6)') AsyRay(iasyray)%InOutPoint(1, 2), AsyRay(iasyray)%InOutPoint(2, 2)
!    ENDIF     
  ENDDO
  
  itype = 2
  DO j = jbeg, jend
    imray = AsyRay(j)%ModRayIdx
    x = ModRay(imray)%InOutPoint(1, :); y = ModRay(imray)%InOutPoint(2, :)
    IF((x(1)*x(2)) .GT. epsm5) CYCLE !IF modular ray does not pass x =0 -> cycle
    nAsyRay = nAsyRay + 1; iasyray = nAsyRay
    AsyRay(iasyray)%PartialAsyRayFlag = itype;
    !Partail Ray <=> normal asy ray mapping
    AsyRay(iasyray)%PartialAsyRayIdx(0) = j; AsyRay(iasyray)%PartialAsyRayIdx(itype) = iasyray
    AsyRay(j)%PartialAsyRayIdx(itype) = iasyray
    !Finding New in-out point
    AsyRay(iasyray)%InOutSurf = AsyRay(j)%InOutSurf
    AsyRay(iasyRay)%InOutPoint = AsyRay(j)%InOutPoint
    AsyRay(iasyRay)%NextRayIdx = AsyRay(j)%NextRayIdx
    IF(tanv .GT. 0) THEN !Fiding New In-Point
      AsyRay(iasyray)%InOutSurf(1) = WEST
      AsyRay(iasyray)%InOutPoint(1, 1) = 0.
      AsyRay(iasyray)%InOutPoint(2, 1) = -tanv*x(1) + y(1)
    ELSE !Fiding New Out-Point
      AsyRay(iasyray)%InOutSurf(2) = EAST
      AsyRay(iasyray)%InOutPoint(1, 2) = 0.
      AsyRay(iasyray)%InOutPoint(2, 2) = -tanv*x(1) + y(1)    
    ENDIF
    !DEBUG OUTPUT
    len =  (AsyRay(iasyray)%InOutPoint(1, 1)-AsyRay(iasyray)%InOutPoint(1, 2))**2 + &
            (AsyRay(iasyray)%InOutPoint(2, 1)-AsyRay(iasyray)%InOutPoint(2, 2))**2
    len = sqrt(len)
    IF(len .lt. epsm5) THEN
      nAsyRay = nAsyRay - 1
      AsyRay(j)%PartialAsyRayIdx(itype) = 0
    ENDIF
  ENDDO
  
  itype =3
  xlim = AsyGeom(3)%x;  ylim = AsyGeom(3)%y
  DO j = jbeg, jend
    imray = AsyRay(j)%ModRayIdx
    x = ModRay(imray)%InOutPoint(1, :); y = ModRay(imray)%InOutPoint(2, :)
    
    IF(AsyRay(j)%PartialAsyRayIdx(1) .EQ. 0 .AND. AsyRay(j)%PartialAsyRayIdx(1) .EQ. 0 )  CYCLE
    
    xnew(1) =  0;  ynew(1) = -tanv*x(1) + y(1) 
    xnew(2) =  x(1) - y(1)/tanv; ynew(2) = 0.
    
    IF((xnew(2) - xlim(1)) * (xnew(2) - xlim(2)) .GT. epsm5 &
         .AND. (ynew(1) - ylim(1)) * (ynew(1) - ylim(2)) .GT. epsm5) CYCLE

    IF(tanv .GT. 0 .and. ABS(ynew(1)) .LT. epsm5) CYCLE  !Exclude 

    nAsyRay = nAsyRay + 1; iasyray = nAsyRay
    AsyRay(iasyray)%PartialAsyRayFlag = itype;
    !Partail Ray <=> normal asy ray mapping
    AsyRay(iasyray)%PartialAsyRayIdx(0) = j; AsyRay(iasyray)%PartialAsyRayIdx(itype) = iasyray
    AsyRay(j)%PartialAsyRayIdx(itype) = iasyray
    !Finding New in-out point
    AsyRay(iasyray)%InOutSurf = AsyRay(j)%InOutSurf
    AsyRay(iasyRay)%InOutPoint = AsyRay(j)%InOutPoint
    AsyRay(iasyRay)%NextRayIdx = AsyRay(j)%NextRayIdx

    IF( (ynew(1) - ylim(1)) * (ynew(1) - ylim(2)) .LE. epsm5) THEN
      IF(tanv .GT. 0) THEN !Fiding New In-Point
        AsyRay(iasyray)%InOutSurf(1) = WEST
        AsyRay(iasyray)%InOutPoint(1, 1) = 0.
        AsyRay(iasyray)%InOutPoint(2, 1) = -tanv*x(1) + y(1)
      ELSE !Fiding New Out-Point
        AsyRay(iasyray)%InOutSurf(2) = EAST
        AsyRay(iasyray)%InOutPoint(1, 2) = 0.
        AsyRay(iasyray)%InOutPoint(2, 2) = -tanv*x(1) + y(1)    
      ENDIF     
    ENDIF
    
    IF((xnew(2) - xlim(1)) * (xnew(2) - xlim(2)) .LE. epsm5) THEN
      AsyRay(iasyray)%InOutSurf(2) = North
      AsyRay(iasyray)%InOutPoint(1, 2) = x(1) - y(1)/tanv
      AsyRay(iasyray)%InOutPoint(2, 2) = 0.    
    ENDIF
    
    !Exclude Zero length partial modular ray
    len =  (AsyRay(iasyray)%InOutPoint(1, 1)-AsyRay(iasyray)%InOutPoint(1, 2))**2 + &
            (AsyRay(iasyray)%InOutPoint(2, 1)-AsyRay(iasyray)%InOutPoint(2, 2))**2
    len = sqrt(len)
    IF(len .lt. epsm5) THEN
      nAsyRay = nAsyRay - 1
      AsyRay(j)%PartialAsyRayIdx(itype) = 0
      CYCLE
    ENDIF     
  ENDDO  
#endif  
  AziAngle(i)%nAsyRay = nAsyRay - AziAngle(i)%AsyRayIdxSt + 1
ENDDO

RayInfo%nAsyRay = nAsyRay 
RayInfo%AsyRay => AsyRay


END SUBROUTINE