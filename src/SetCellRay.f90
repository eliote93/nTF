SUBROUTINE CellRayGen(RayInfo, CellRayBase, AsyGeom, lEdge)
!#define DBG_Cellray
USE PARAM
USE typedef,        ONLY : Cell_type,        CellRayBase_Type,  RayInfo_Type,             &
                           CellRayInfo_type, AziAngleInfo_Type, AsyRayInfo_type,          &
                           ModRayInfo_type,  BasicGeom,         CellRayBase_SearchTree        
USE GEOM,           ONLY : nAsy,             nAsyType,          nCellType,                &
                           nPinType,         nz,                Core,                     &
                           AsyInfo,          Asy,               AsyInfo,                  &
                           PinInfo,          CellInfo
USE GeomTreatment,  ONLY : Cell_RayIntersection, RectAsyPinIdx, Find_Surface,             &
                          FindCellRayBase
USE SetRay,         ONLY : Copy_CellRayInfo,     VolCorChk             
USE IOUTIL,         ONLY : Terminate
USE BasicOperation, ONLY : CP_CA
USE ALLOCS
IMPLICIT NONE
TYPE(RayInfo_Type) :: RayInfo
TYPE(CellRayBase_Type),POINTER :: CellRayBase(:)
Type(BasicGeom) :: AsyGeom(0:3)
LOGICAL :: lEdge

TYPE(AsyRayInfo_Type), POINTER :: AsyRay(:)
TYPE(AziAngleInfo_Type), POINTER :: AziAng(:)
TYPE(CellRayBase_SearchTree), POINTER :: SearchTree

INTEGER :: nAsyRay, nxya
INTEGER :: itype                          !Partial Assembly Type
INTEGER :: i, j, jbeg, jend, k, l, m
INTEGER :: iang, imray, iaray, ireg, ireg0, iasy, ipin, icel, ix, iy, IOsurf(2)
INTEGER :: npt, nSeg, nAziAng, nx, ny, nxy, AsyType
LOGICAL, POINTER :: PinUsed(:, :, :)        
REAL :: x(2), y(2), tanv                  !Assembly In & Out Points, tanv       
REAL :: RAY(3)                            !Assembly Ray Equation
REAL :: pts(2, 100)                       !Intersection Points
REAL :: CellRayPts(2, 100) 
REAL :: hx(0:500), hy(0:500), cx(1000, 0:3), cy(1000, 0:3), cx0, cy0
INTEGER, POINTER :: bcrmap(:,:,:), sympin(:) ! azi, global pin, rayidx
INTEGER :: jang, jreg, jj
INTEGER :: ibCel

#ifdef DBG_Cellray
REAL :: Err
#endif
!#define DBG_VolCor
#ifdef DBG_VolCor
LOGICAL :: lmod
INTEGER :: ndiv
!INTEGER,POINTER :: idx(:), IDXMAP(:,:,:) !ipin,idx(300), IDXMAP(300,1000,3) !ipin,

INTEGER :: kk,jjj,kkk
INTEGER :: ncellray, iceray
#endif

AziAng => RayInfo%AziAngle
AsyRay => RayInfo%AsyRay

nAziAng = RayInfo%nAziAngle
nAsyRay = RayInfo%nAsyRay

!PinCell Center Location
nxy = AsyGeom(0)%nx * AsyGeom(0)%ny
!ALLOCATE(CellRayStorage(nxy, 0:3))
ALLOCATE(PinUsed(nCellType, nxy, 0:3))
ALLOCATE(BCRMAP(nAziAng, nxy, 0:1000),sympin(nxy))
bcrmap=0;
!ALLOCATE(CELLRayStorage0(nCellType))
DO itype =0, 3
  IF(LEdge .AND. itype .NE. 0) CYCLE
  nx = AsyGeom(itype)%nx; ny = AsyGeom(itype)%ny
  hx(0) = AsyGeom(itype)%x(1); hx(nx) = AsyGeom(itype)%x(2)
  hx(1:nx-1) = -AsyGeom(itype)%Line(3,1:nx-1)
  hy(0) = AsyGeom(itype)%y(1); hy(ny) = AsyGeom(itype)%y(2)
  hy(1:ny-1) = -AsyGeom(itype)%Line(3,nx:nx+ny-2)
  i= 0 
  DO iy = 1, ny
    DO ix = 1, nx
      i = i + 1
      cx(i, itype) = (hx(ix-1) + hx(ix))*Half
      cy(i, itype) = (hy(ny - iy + 1) + hy(ny - iy))*Half
      PinUsed(:, i, itype) = FALSE
    ENDDO
  ENDDO
ENDDO
i=0;
DO iy=1, ny
    DO ix = 1, nx
        i=i+1
        sympin(i)=(ix-1)*ny+iy
    ENDDO
ENDDO
!DO i = 1,AsyGeom(0)%nx-1
!  WRITE(81, *)  -AsyGeom(0)%Line(3,i), -20.0
!  WRITE(81, *)  -AsyGeom(0)%Line(3,i),  20.0
!ENDDO
!DO i = AsyGeom(0)%nx, AsyGeom(0)%nx+AsyGeom(0)%ny-2
!   WRITE(81, *)   -20.0, -AsyGeom(0)%Line(3,i)
!   WRITE(81, *)    20.0, -AsyGeom(0)%Line(3,i)
!ENDDO

nxya = Core%nxya

DO iasy = 1, nxya
  iType = Asy(iasy)%PartialAsyFlag
  AsyType = Asy(iasy)%AsyType
  nxy = AsyInfo(AsyType)%nxy
  DO i = 1, nxy
    ipin = AsyInfo(AsyType)%Pin(i)
    !PinUsed(:, i, iType) = FALSE
    DO k = 1, nz
      icel = PinInfo(ipin)%cell(k)
      PinUsed(icel, i, iType) = TRUE
    ENDDO
  ENDDO
ENDDO


!Find Intersect Point between Ray and Pin Cell Boundary
DO i = 1, nAsyRay                        !Cell Ray Sweep
  iang = AsyRay(i)%AziAngIdx                                     !Azimuthal Angle Index
  imray =  AsyRay(i)%ModRayIdx                                   !Modular Ray Angle Index
  itype = AsyRay(i)%PartialAsyRayFlag
  !IF(itype .NE. 0) CYCLE
  x = AsyRay(i)%InOutPoint(1, :); y = AsyRay(i)%InOutPoint(2, :) !Modular Ray In & Out Point
  !
  tanv = AziAng(iang)%tanv
  RAY(1) = tanv; RAY(2) = -1.; RAY(3) = -tanv * x(1) + y(1)
  CALL Cell_RayIntersection(AsyGeom(itype), RAY, pts, npt, epsm6)
  nSeg = npt + 1
  
  CellRayPts(:, 1) = AsyRay(i)%InOutPoint(:, 1)
  CellRayPts(:, 2:nseg) = pts(:, 1:npt)
  CellRayPts(:, nseg+1) = AsyRay(i)%InOutPoint(:, 2)  
  AsyRay(i)%nCellRay = nSeg
  CALL Dmalloc(AsyRay(i)%PinRaySurf,2,nSeg)
  CALL Dmalloc(AsyRay(i)%PinIdx, nSeg)
  CALL Dmalloc(AsyRay(i)%PinRayIdx, nSeg)
  !REGION SEARCH
  jj=0
  DO j = 1, nSeg
    !Mid - Point of Cell Ray
    x(1) = (CellRayPts(1, j) + CellRayPts(1, j+1)) * HALF
    x(2) = (CellRayPts(2, j) + CellRayPts(2, j+1)) * HALF
    iReg = RectAsyPinIdx(AsyGeom(itype), x)                   !PinCell Index
    !IF(iType .GT. 0) iReg0 = RectAsyPinIdx(AsyGeom(0), x)                   
    DO icel = 1, nCellType
      IF(PinUsed(icel, ireg, itype) ) EXIT
    ENDDO
    !IF(icel .GT. nCellType) Call Terminate('SetCellRay.f90 : Can not find the Cell type')
    IF(icel .GT. nCellType) CYCLE
    !Coordinate Center(Origin) Change : Assembly Center -> Cell Center    
    cx0  = cx(ireg, itype)+(2.*CellInfo(icel)%Geom%cx-CellInfo(icel)%Geom%x(1)-CellInfo(icel)%Geom%x(2))*HALF
    cy0  = cy(ireg, itype)+(2.*CellInfo(icel)%Geom%cy-CellInfo(icel)%Geom%y(1)-CellInfo(icel)%Geom%y(2))*HALF    
    !Incomming and outgoing points coordinates change
    x = CellRayPts(1, j:j+1) - cx0
    y = CellRayPts(2, j:j+1) - cy0    
    !In-Out Surface Searching
    IOSurf(1) = Find_Surface(CellInfo(icel)%Geom,(/x(1), y(1)/))
    IOSurf(2) = Find_Surface(CellInfo(icel)%Geom,(/x(2), y(2)/))

    !BaseCellRay Searching
!    ibcel=CellRayBase(icel)%CellType
    ibcel = CellInfo(icel)%basecellstr
    SearchTree => CellRayBase(ibcel)%SearchTree(iang)
    l = FindCellRayBase(SearchTree, (/x(1), y(1)/), IOSurf(1))
    l = l + CellRayBase(ibcel)%AziAngIdx(iang)-1        
    AsyRay(i)%PinRaySurf(:, j) = IOSurf
    AsyRay(i)%PinIdx(j) = ireg; AsyRay(i)%PinRayIdx(j) = l
        !l = FindCellRayBase(SearchTree, (/x(1), y(1)/), IOSurf(:),iang,nAziAng)
        !if( l .NE. 0 )THEN
        !l = l + CellRayBase(ibcel)%AziAngIdx(iang)-1        
        !endif
        !IF( bcrmap(iang, ireg, 1) .EQ. 0 )THEN
        !    IF(   iang.GT. nAziAng*1/4 .AND. iang.LE. nAziAng*2/4)THEN
        !        jang=nAziAng/2+1-iang
        !        jreg=sympin(ireg)
        !        k=bcrmap(jang,jreg,0)
        !        l=bcrmap(jang,jreg,1)-CellRayBase(ibcel)%AziAngIdx(jang)+1
        !        l = l + CellRayBase(ibcel)%AziAngIdx(iang)-1   
        !    ELSEIF(iang .GT. nAziAng*3/4)THEN
        !        jang=nAziAng*3/2+1-iang
        !        jreg=sympin(ireg)
        !        k=bcrmap(jang,jreg,0)
        !        l=bcrmap(jang,jreg,k)-CellRayBase(ibcel)%AziAngIdx(jang)+1
        !        l=cellraybase(ibcel)%nlines(iang)-l+1
        !        l = l + CellRayBase(ibcel)%AziAngIdx(iang)-1   
        !    ENDIF
        !    bcrmap(iang, ireg, 1)= l
        !    bcrmap(iang, ireg, 0)= 1                
        !ELSE
        !    k=bcrmap(iang, ireg, 0)
        !    bcrmap(iang, ireg, k+1)=bcrmap(iang, ireg, k)+searchtree%ndiv
        !    l=bcrmap(iang, ireg, k+1)
        !    bcrmap(iang, ireg, 0)=bcrmap(iang, ireg, 0)+1
        !ENDIF
        !IF( l .GT. CellRayBase(ibcel)%AziAngIdx(iang)+CellRayBase(ibcel)%nlines(iang)-1 .OR. l .LT. CellRayBase(ibcel)%AziAngIdx(iang))THEN
        !    AsyRay(i)%nCellRay=AsyRay(i)%nCellRay-1
        !    nseg=nseg-1
        !ELSE
        !    jj=jj+1
        !    AsyRay(i)%PinRaySurf(:, jj) = IOSurf
        !    AsyRay(i)%PinIdx(jj) = ireg; AsyRay(i)%PinRayIdx(jj) = l
        !ENDIF
#ifdef DBG_Cellray
    m = CellRayBase(icel)%CellRay(l)%nseg
    err = 0
    err = (x(1)-CellRayBase(icel)%CellRay(l)%pts(1, 1))**2 + (y(1)-CellRayBase(icel)%CellRay(l)%pts(2, 1))**2
    err = err + (x(2)-CellRayBase(icel)%CellRay(l)%pts(1, m+1))**2 + (y(2)-CellRayBase(icel)%CellRay(l)%pts(2, m+1))**2
    err = SQRT(err)
    IF(err .gt. 0.01_8) THEN
      CALL Terminate('SetCellRay.f90 : Can Not fine the CellRay')
    ENDIF  
#endif    
  ENDDO
ENDDO

!----temp
#ifdef DBG_VolCor
!IF( lMod )THEN
write(*,*) 'DBG_VolCor ON'
ALLOCATE(idx(nxy),idxmap(nxy,1000,3))
ndiv=50
DO i = 1, nAziAng  
 idx=0
  jbeg = AziAng(i)%AsyRayIdxSt; jend = jbeg + AziAng(i)%nAsyRay - 1
  DO j  = jbeg, jend ! assembly ray sweep
    IF(AsyRay(j)%PartialAsyRayFlag .NE. itype) CYCLE
    nCellRay = AsyRay(j)%nCellRay            !Number of CellRay
    DO k = 1, nCellRay ! fractional cellray in AsyRay
      ipin = AsyRay(j)%PinIdx(k)      !pin    index of cellray
      iceray = AsyRay(j)%PinRayIdx(k) !pinray index of cellray
      idx(ipin)=idx(ipin)+1
      idxmap(ipin,idx(ipin),1)=j
      idxmap(ipin,idx(ipin),2)=k
      idxmap(ipin,idx(ipin),3)=iceray
    ENDDO
  ENDDO
  DO ipin=1,nxy
    DO icel = 1, nCellType  !Cell Type Sweep
      IF(.NOT. PinUsed(icel, ipin, 0)) CYCLE
        DO j=2,idx(ipin)
            jj=idxmap(ipin,j,1);jjj=idxmap(ipin,j-1,1)
            kk=idxmap(ipin,j,2);kkk=idxmap(ipin,j-1,2)
            AsyRay(jj)%PinRayIdx(kk)=AsyRay(jjj)%PinRayIdx(kkk)+ndiv
            !IF(AsyRay(jj)%PinRayIdx(kk).GT. cellraybase(icel)%nline  )THEN
                !AsyRay(jj)%PinRayIdx(kk)=0
            !ENDIF        
        ENDDO        
    ENDDO
  ENDDO
ENDDO
DEALLOCATE(idx,idxmap)
!ENDIF
#endif

!#define vc_chk   ! BYS edit 16/02/26 --- volcorrection debug on
#ifdef vc_chk
CALL VolCorCHK(AsyGeom(0), RayInfo, CellInfo, PinUsed(:, :, 0), nCellType, 0)
IF(.NOT. ledge) THEN
  CALL VolCorCHK(AsyGeom(1), RayInfo, CellInfo, PinUsed(:, :, 1), nCellType, 1)
  CALL VolCorCHK(AsyGeom(3), RayInfo, CellInfo, PinUsed(:, :, 3), nCellType, 2)
  CALL VolCorCHK(AsyGeom(3), RayInfo, CellInfo, PinUsed(:, :, 3), nCellType, 3)
ENDIF
#endif

DEALLOCATE(PinUsed)

END SUBROUTINE

SUBROUTINE Copy_CellRayInfo(CellRay0, CellRay)
!Copy CellRay structures
!CellRay -> CellRay
USE typedef, only : CellRayInfo_type
USE BasicOperation, only : CP_CA, CP_VA
USE ALLOCS, ONLY : Dmalloc
IMPLICIT NONE
TYPE(CellRayInfo_type), intent(IN) :: CellRay0
TYPE(CellRayInfo_type), intent(OUT) :: CellRay
INTEGER :: nseg

nseg = CellRay0%nseg
CellRay%nSeg = nSeg
CALL Dmalloc(CellRay%LocalFsrIdx, nseg); CALL Dmalloc(CellRay%LenSeg, nseg)

CALL CP_VA(CellRay%LocalFsrIdx, CellRay0%LocalFsrIdx, nseg)
CALL CP_VA(CellRay%LenSeg, CellRay0%LenSeg, nseg)

!CALL Dmalloc(CellRay%LenSegInt, nseg)

END SUBROUTINE