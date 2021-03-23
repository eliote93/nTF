SUBROUTINE CoreRayGen(RayInfo, CoreRay, Core, Asy, AsyGeom, lEdge)
USE PARAM
USE TYPEDEF,        ONLY : RayInfo_type,      ModRayInfo_type,    AsyRayInfo_type,   &
                           AziAngleInfo_Type, CoreRayInfo_type,                      &
                           CoreInfo_type,     Asy_Type,           BasicGeom
USE BasicOperation, ONLY : CP_CA, CP_VA
USE ALLOCS
IMPLICIT NONE                        
TYPE(RayInfo_type) :: RayInfo
TYPE(CoreRayInfo_type), POINTER  :: CoreRay(:)       !Core Ray variable which will be 
TYPE(CoreInfo_type) :: Core                          !Croe Information of the problem
TYPE(Asy_Type),POINTER :: Asy(:)                     !Information of Assemblies which resides inside of core
TYPE(BasicGeom)  :: AsyGeom(0:3)                     !Assembly Geometries Information
Logical :: lEdge


TYPE CoreRay_beg_type
  INTEGER :: ixa, iya, itype
  INTEGER :: iasyray, iang, isurf
ENDTYPE

TYPE(AziAngleInfo_Type), POINTER :: AziAng(:)
TYPE(AsyRayInfo_type), POINTER :: AsyRay(:)
TYPE(ModRayInfo_type), POINTER :: ModRay(:)
!TYPE(CoreRayInfo_type), POINTER  :: CoreRay0(:)          !Intermidate Saving variable
TYPE(CoreRay_beg_type), POINTER :: CoreRay_beg(:)
INTEGER, POINTER :: CoreIdx(:,:)
INTEGER :: nAsyRayinSurf(6, 0:3, 100)                    !Number of Assembly rays which have same incomming surface
INTEGER :: AsyRayinSurfFstIdx(6, 0:3, 100)                !The index of the First Assembly Ray which have same incomming surface
INTEGER :: nAziAng, nxa, nya, nbd, IOSurf(2), isurf, osurf
INTEGER :: nCoreRays(100), nCoreRay, nRay                !Number of Core Rays
INTEGER :: AsyRayIdx(5000), AsyIdx(5000)                 !Temp saving variables


INTEGER :: ixaIdx(5000), iyaIdx(5000)
INTEGER :: Asy_Extd(0:200, 0:200), Dummy_Asy(0:3, 0:3)
INTEGER :: iang, itype, ixa, iya
INTEGER :: i, j, jbeg, jend, jnew
INTEGER :: iasyray, icoreray, iray
INTEGER :: NeighIdx(2,4)
REAL :: tanv 

data Dummy_Asy /0,  1,  0,  1,  &
                0, -1,  0, -1,  &
                2, -1, -1, -1,  &
                2, -1, -1, -1/
DATA NeighIdx(:,SOUTH) /0,  1/; DATA NeighIdx(:,NORTH) /0, -1/
DATA NeighIdx(:,WEST) /-1,  0/; DATA NeighIdx(:,EAST)  /1,  0/                
EQUIVALENCE(IOSurf(1), isurf);EQUIVALENCE(IOSurf(2), osurf)

AziAng => RayInfo%AziAngle
AsyRay => RayInfo%AsyRay
ModRay => RayInfo%ModRay
CoreIdx => Core%CoreIdx

nAziAng = RayInfo%nAziAngle
nxa = Core%nxa; nya = Core%nya

!Calculate the Number of Assembly rays which have same incomming surface
DO iang = 1, nAziAng
  nAsyRayinSurf(:, :, iang) = 0
!  IF(lEdge .and. itype .NE. 0) CYCLE
  jbeg = AziAng(iang)%AsyRayIdxSt; jend = jbeg + AziAng(iang)%nAsyRay - 1
  DO j = jbeg, jend
    itype =AsyRay(j)%PartialAsyRayFlag; IOSurf = AsyRay(j)%InOutsurf
    nAsyRayinSurf(ISurf, itype, iang) = nAsyRayinSurf(ISurf, itype, iang) + 1
    continue
  ENDDO
  
  !The first Assembly Ray Index
  DO itype = 0, 3
    IF(itype .NE. 0 .AND. LEdge) CYCLE
    DO isurf = 1, 4
      AsyRayInSurfFstIdx(isurf, itype, iang) = -1
      DO j = jbeg, jend
        IF(AsyRay(j)%PartialAsyRayFlag .NE. ITYPE) CYCLE
        IF(isurf .EQ. AsyRay(j)%InOutSurf(1)) EXIT
      ENDDO
      IF(j .LE. jend) AsyRayInSurfFstIdx(isurf, itype, iang) = j
    ENDDO
  ENDDO
  CONTINUE
ENDDO

!Extended Assembly Map(including dummy ) in the core
Asy_Extd(0:nxa, 0:nya) = 0.
DO iya = 1, nya
  DO ixa = 1, nxa
    IF(CoreIdx(ixa, iya) .NE. 0) THEN
      j =CoreIdx(ixa, iya) 
      Asy_Extd(ixa, iya) = Asy(j)%PartialAsyFlag
    ELSE
      Asy_Extd(ixa, iya) = Dummy_Asy(Asy_Extd(ixa-1, iya),Asy_Extd(ixa, iya-1))
    ENDIF      
  ENDDO
ENDDO
!Calculate the number of Core Ray
nCoreRay = 0
DO iang = 1, nAziAng
  nCoreRays(iang) = 0
  
  DO ixa = 1, nxa
    itype = Asy_Extd(ixa, 1) !Assembly Geometry Type
    nCoreRays(iang) =  nCoreRays(iang) + nAsyRayinSurf(NORTH, itype, iang)
    itype = Asy_Extd(ixa, nya) !Assembly Geometry Type
    nCoreRays(iang) =  nCoreRays(iang) + nAsyRayinSurf(SOUTH, itype, iang)
  ENDDO

  DO iya = 1, nya
    itype = Asy_Extd(1, iya) !Assembly Geometry Type
    nCoreRays(iang) =  nCoreRays(iang) + nAsyRayinSurf(WEST, itype, iang)
    itype = Asy_Extd(nxa, iya) !Assembly Geometry Type
    nCoreRays(iang) =  nCoreRays(iang) + nAsyRayinSurf(EAST, itype, iang)
  ENDDO
  nCoreRay = nCoreRay + nCoreRays(iang)
ENDDO
CONTINUE

ALLOCATE(CoreRay(nCoreRay))
ALLOCATE(CoreRay_beg(nCoreRay))

RayInfo%nCoreRay = nCoreRay
!Generate the Beginning Assembly
icoreray = 0
DO iang = 1, nAziAng
  jbeg = 1; jend = nCoreRays(iang)
  tanv = AziAng(iang)%tanv
  AziAng(iang)%CoreRayIdxSt = iCoreRay  + 1
  IF(tanv .GT. 0.) THEN
    iya = nya
    DO ixa = nxa, 1, -1
      itype = Asy_Extd(ixa, iya)
      jbeg = AsyRayInSurfFstIdx(SOUTH, itype, iang); jend = jbeg + nAsyRayInSurf(SOUTH, itype, iang) - 1  
      IF( AsyRayInSurfFstIdx(SOUTH, itype, iang) .EQ. -1) CYCLE
      DO j = jbeg, jend
        iCoreRay = iCoreRay + 1; CoreRay_Beg(iCoreRay)%itype = itype
        CoreRay_beg(iCoreRay)%ixa = ixa; CoreRay_beg(iCoreRay)%iya = iya         !Location Information
        CoreRay_Beg(iCoreRay)%iAsyRay = j; CoreRay_Beg(iCoreRay)%iang = iang     !Assembly Ray Index and corresponding Angle
        CoreRay_Beg(iCoreRay)%ISurf = SOUTH
      ENDDO
    ENDDO
    ixa = 1
    DO iya = nya, 1, -1
      itype = Asy_Extd(ixa, iya)
      jbeg = AsyRayInSurfFstIdx(WEST, itype, iang); jend = jbeg + nAsyRayInSurf(WEST, itype, iang) - 1  
      IF (AsyRayInSurfFstIdx(WEST, itype, iang) .EQ. -1) CYCLE
      DO j = jbeg, jend
        iCoreRay = iCoreRay + 1; CoreRay_Beg(iCoreRay)%itype = itype
        CoreRay_beg(iCoreRay)%ixa = ixa; CoreRay_beg(iCoreRay)%iya = iya         !Location Information
        CoreRay_Beg(iCoreRay)%iAsyRay = j; CoreRay_Beg(iCoreRay)%iang = iang     !Assembly Ray Index and corresponding Angle
        CoreRay_Beg(iCoreRay)%ISurf = WEST
      ENDDO
    ENDDO
  ELSE  ! Angle > PI/2
    iya = nya
    DO ixa = 1, nxa 
      itype = Asy_Extd(ixa, iya)
      jbeg = AsyRayInSurfFstIdx(SOUTH, itype, iang); jend = jbeg + nAsyRayInSurf(SOUTH, itype, iang) - 1  
      IF (AsyRayInSurfFstIdx(SOUTH, itype, iang) .EQ. -1) CYCLE
      DO j = jbeg, jend
        iCoreRay = iCoreRay + 1; CoreRay_Beg(iCoreRay)%itype = itype
        CoreRay_beg(iCoreRay)%ixa = ixa; CoreRay_beg(iCoreRay)%iya = iya         !Location Information
        CoreRay_Beg(iCoreRay)%iAsyRay = j; CoreRay_Beg(iCoreRay)%iang = iang     !Assembly Ray Index and corresponding Angle
        CoreRay_Beg(iCoreRay)%ISurf = SOUTH
      ENDDO
    ENDDO
    ixa = nxa
    DO iya = nya, 1, -1
      itype = Asy_Extd(ixa, iya)
      jbeg = AsyRayInSurfFstIdx(EAST, itype, iang); jend = jbeg + nAsyRayInSurf(EAST, itype, iang) - 1  
      IF (AsyRayInSurfFstIdx(EAST, itype, iang) .EQ. -1) CYCLE
      DO j = jbeg, jend
        iCoreRay = iCoreRay + 1; CoreRay_Beg(iCoreRay)%itype = itype
        CoreRay_beg(iCoreRay)%ixa = ixa; CoreRay_beg(iCoreRay)%iya = iya         !Location Information
        CoreRay_Beg(iCoreRay)%iAsyRay = j; CoreRay_Beg(iCoreRay)%iang = iang     !Assembly Ray Index and corresponding Angle
        CoreRay_Beg(iCoreRay)%ISurf = EAST
      ENDDO
    ENDDO  
  ENDIF
  !Number of Core Ray in a given azimuthal ANgle
  AziAng(iang)%nCoreRay = iCoreRay - AziAng(iang)%CoreRayIdxSt + 1
  CONTINUE
ENDDO

!Generate Core Ray
DO i = 1, nCoreRay
  !Begining Points Data
  icoreray = i; iang = CoreRay_Beg(iCoreRay)%iang
  ixa = CoreRay_Beg(i)%ixa; iya = CoreRay_Beg(i)%iya
  j = CoreRay_Beg(iCoreRay)%iAsyRay; iray =1;
  ISurf = CoreRay_Beg(iCoreRay)%ISurf
  iType = Asy_Extd(ixa, iya)
  AsyRayIdx(iray) = j
  AsyIdx(iray) = CoreIdx(ixa, iya)
  ixaIdx(iray) = ixa; iyaIdx(iray) = iya
  DO WHILE(TRUE)
    !Find Next Core Ray
    OSURF = AsyRay(j)%InOutSurf(2);
    ixa = ixa + NeighIdx(1, OSurf); iya = iya + NeighIdx(2, OSurf)
    IF( (ixa - nxa) * (ixa - 1) .GT. 0) EXIT; IF( (iya - nya) * (iya - 1) .GT. 0) EXIT
    !Finding Next Assembly Ray Idx
    itype = Asy_Extd(ixa, iya)
    jnew = AsyRay(j)%NextRayIdx(2); 
    IF(itype .NE. 0) jnew = AsyRay(jnew)%PartialAsyRayIdx(itype)
    IF(jnew .eq. 0) exit
    iRay = iRay + 1; 
    
    AsyRayIdx(iray) = jnew; j = jnew
    AsyIdx(iray) = CoreIdx(ixa, iya)
    
    ixaIdx(iray) = ixa; iyaIdx(iray) = iya
  ENDDO
  !Saving Data
  CoreRay(iCoreRay)%nRay = iRay;
  CALL Dmalloc(CoreRay(iCoreRay)%AsyRayIdx, iRay)
  CALL Dmalloc(CoreRay(iCoreRay)%AsyIdx, iRay)

  CoreRay(iCoreRay)%AsyRayIdx(1:iRay) = AsyRayIdx(1:iRay)
  CoreRay(iCoreRay)%AsyIdx(1:iRay) = AsyIdx(1:iRay)
  
  CALL Dmalloc(CoreRay(iCoreRay)%ixa, iRay)
  CALL Dmalloc(CoreRay(iCoreRay)%iya, iRay)
  CoreRay(icoreray)%ixa(1:iray) = IxaIdx(1:iRay)
  CoreRay(icoreray)%iya(1:iray) = IyaIdx(1:iray)

  CoreRay(iCoreRay)%InOutSurf = (/ISURF, OSURF/)
  CoreRay(icoreRay)%iang = iang
  CONTINUE
  !ALLOCATE(CoreRay(iCoreRay)%AsyRayIdx())
ENDDO

RayInfo%CoreRay => CoreRay
DEALLOCATE(CoreRay_Beg)
#ifdef CrayOut
jbeg = AziAng(1)%CoreRayIdxSt
jend = AziAng(1)%CoreRayIdxSt + AziAng(1)%nCoreRay - 1
DO i = jbeg, jend
  CALL WRITE_CoreRay(i, i)
ENDDO
jbeg = AziAng(8)%CoreRayIdxSt
jend = AziAng(8)%CoreRayIdxSt + AziAng(8)%nCoreRay - 1
DO i = jbeg, jend
  CALL WRITE_CoreRay(i, i)
ENDDO
#endif
END SUBROUTINE

SUBROUTINE WRITE_CORERAY(i, idx, dir, ldetail)
USE PARAM
USE TYPEDEF
USE GEOM
USE RAYS
IMPLICIT NONE
INTEGER :: idx
LOGICAL :: lFirst
data lFirst /TRUE/
REAL, SAVE :: CX(100), CY(100),HX(100),HY(100)
REAL ::x(2),y(2)
INTEGER :: ixa, iya, iray, nRay, dir
INTEGER :: i, j, k
logical :: lDetail
IF(lFirst) Then
  IF(Asy(1)%PartialAsyFlag .EQ. 3) THEN
    CX(1) = 0; CY(1) = 0
    HX(1) = AsyGeom(3)%lx; HY(1) = AsyGeom(3)%ly
  ELSE
    CX(1) = AsyGeom(0)%lx/2; CY(1) = -AsyGeom(0)%ly/2
    HX(1) = AsyGeom(0)%lx; HY(1) = AsyGeom(0)%ly
  ENDIF
  
  DO ixa = 2, Core%nxa
    CX(ixa) = CX(ixa-1) + AsyGeom(0)%lx
    HX(ixa) = Asygeom(0)%lx
  ENDDO

  DO iya = 2, Core%nya
    CY(iya) = CY(iya-1) - AsyGeom(0)%ly
    HY(iya) = Asygeom(0)%ly
  ENDDO
  WRITE(92,'(100e20.5)') HX(1:Core%nxa)
  WRITE(92,'(100e20.5)') HY(1:Core%nya)
  lFirst = .FALSE.
ENDIF

nRay = CoreRay(i)%nRay
IF(ldetail) THEN
  DO j = 1, nRay
    iray = CoreRay(i)%AsyRayIdx(j);
    ixa = CoreRay(i)%ixa(j); iya = CoreRay(i)%iya(j)
    x(1) = Asyray(iray)%InOutPoint(1, 1); x(2) = Asyray(iray)%InOutPoint(1, 2)
    y(1) = Asyray(iray)%InOutPoint(2, 1); y(2) = Asyray(iray)%InOutPoint(2, 2)
    x = x + CX(ixa); y = y + CY(iya)
    WRITE(91, '(2e20.7, 2I9)') x(1), y(1), idx, dir
    WRITE(91, '(2e20.7, 2I9)') x(2), y(2), idx, dir
  ENDDO
else
  j = 1
  iray = CoreRay(i)%AsyRayIdx(j);
  ixa = CoreRay(i)%ixa(j); iya = CoreRay(i)%iya(j)
  x(1) = Asyray(iray)%InOutPoint(1, 1); y(1) = Asyray(iray)%InOutPoint(2, 1)
  x(1) = x(1) + CX(ixa); y(1) = y(1) + CY(iya)
  j = nray
  iray = CoreRay(i)%AsyRayIdx(j);
  ixa = CoreRay(i)%ixa(j); iya = CoreRay(i)%iya(j)
  x(2) = Asyray(iray)%InOutPoint(1, 2); y(2) = Asyray(iray)%InOutPoint(2, 2)
  x(2) = x(2) + CX(ixa); y(2) = y(2) + CY(iya)
  IF(dir .eq. 1) THEN
    WRITE(91, '(2e20.7, 2I9)') x(1), y(1), idx, dir
    WRITE(91, '(2e20.7, 2I9)') x(2), y(2), idx, dir
  ELSE
    WRITE(91, '(2e20.7, 2I9)') x(2), y(2), idx, dir
    WRITE(91, '(2e20.7, 2I9)') x(1), y(1), idx, dir
  ENDIF
endif
END SUBROUTINE