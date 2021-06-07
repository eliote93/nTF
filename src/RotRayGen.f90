#define RotRayDBG
#include <defines.h>
SUBROUTINE RotationRayGen(RayInfo, RotRay, Core, Asy, lEdge)
USE PARAM
USE TYPEDEF,        ONLY : RayInfo_type,      ModRayInfo_type,    AsyRayInfo_type,   &
                           AziAngleInfo_Type, CoreRayInfo_type,   RotRayInfo_Type,   &
                           CoreInfo_type,     Asy_Type,           BasicGeom
USE IOUTIL,         ONLY : Terminate
USE BasicOperation, ONLY : CP_VA
USE ALLOCS
USE cntl,           ONLY : nTracerCntl
IMPLICIT NONE

TYPE(RayInfo_type) :: RayInfo
TYPE(RotRayInfo_Type), POINTER  :: RotRay(:)       !Core Ray variable which will be
TYPE(CoreInfo_type) :: Core                        !Croe Information of the problem
TYPE(Asy_Type),POINTER :: Asy(:)                   !Information of Assemblies which resides inside of core
LOGICAL :: lEdge

TYPE CoreRaySearchTree_Type
  REAL :: X(2), Y(2)
  INTEGER :: ixa(2), iya(2)
  INTEGER :: iang
  INTEGER :: IOSURF(2)
  LOGICAL :: Lused = .FALSE.
END TYPE

TYPE(RotRayInfo_Type), POINTER :: RotRay0(:)       !Rotational Ray
TYPE(AziAngleInfo_Type), POINTER :: AziAng(:)
TYPE(AsyRayInfo_type), POINTER :: AsyRay(:)
TYPE(CoreRayInfo_type), POINTER :: CoreRay(:)
TYPE(CoreRaySearchTree_Type), POINTER :: SearchTree(:)

INTEGER :: nCoreRay, nRotRay, nRay, nAziAng, nxa, nya
INTEGER, POINTER :: RayIdx(:), DirIdx(:)
INTEGER, POINTER :: PhiangBegIdx(:, :), PhiangEndIdx(:,:)               !incomming and Outgoing save adress
INTEGER :: RayRefAng(1000), RayRotAng(1000), mp(2)
INTEGER :: ISurf, OSurf, RefSurf, IOSurf(2)                             !Surface
INTEGER :: dir, iray, iCoreRay, ixa, iya, iang, irefang
INTEGER :: i ,j ,k, m, jbeg, jend, kbeg, kend, jnew
INTEGER :: iFlag, ibegidx, iEndIdx, nSvIdx
REAL :: xref, yref, x, y, xr(2), yr(2), Err
LOGICAL :: ForwardBC, BackwardBC, AllRef
EQUIVALENCE(ISurf, IOSurf(1)); EQUIVALENCE(OSurf, IOSurf(2))
DATA MP /2, 1/
AziAng => RayInfo%AziAngle
AsyRay => RayInfo%AsyRay
CoreRay => RayInfo%CoreRay
nCoreRay = RayInfo%nCoreRay
nAziAng = RayInfo%nAziAngle
AllRef = TRUE
nxa = Core%nxa; nya = Core%nya

DO i = 1, 4
  AllRef = AllRef .and. Core%RadSym(i)
ENDDO

!Making Search Tree

ALLOCATE(SearchTree(nCoreRay))
xr = 0; yr = 0
DO i =1, nCoreRay
  nRay = CoreRay(i)%nRay
  !iang = CoreRay(i)%iang
  SearchTree(i)%IOSURF = CoreRay(i)%InOutSurf
  jbeg = CoreRay(i)%AsyRayIdx(1); jend = CoreRay(i)%AsyRayIdx(nray)

  SearchTree(i)%ixa(1) = CoreRay(i)%ixa(1); SearchTree(i)%ixa(2) = CoreRay(i)%ixa(nRay)
  SearchTree(i)%iya(1) = CoreRay(i)%iya(1); SearchTree(i)%iya(2) = CoreRay(i)%iya(nRay)

  ixa =  SearchTree(i)%ixa(1); iya = SearchTree(i)%iya(1)
  SearchTree(i)%X(1) = AsyRay(jbeg)%InOutPoint(1,1) + Core%AsyCentX(ixa)
  SearchTree(i)%Y(1) = AsyRay(jbeg)%InOutPoint(2,1) + Core%AsyCentY(iya)
  xr(1) = min(xr(1), SearchTree(i)%X(1))
  ixa =  SearchTree(i)%ixa(2); iya = SearchTree(i)%iya(2)
  SearchTree(i)%X(2) = AsyRay(jend)%InOutPoint(1,2) + Core%AsyCentX(ixa)
  SearchTree(i)%Y(2) = AsyRay(jend)%InOutPoint(2,2) + Core%AsyCentY(iya)

  xr(1) = min(xr(1), SearchTree(i)%X(1)); xr(1) = min(xr(1), SearchTree(i)%X(2))
  xr(2) = max(xr(1), SearchTree(i)%X(1)); xr(2) = max(xr(2), SearchTree(i)%X(2))
  yr(1) = max(yr(1), SearchTree(i)%Y(1)); yr(1) = max(yr(1), SearchTree(i)%Y(2))
  yr(2) = min(yr(2), SearchTree(i)%Y(1)); yr(2) = min(yr(2), SearchTree(i)%Y(2))
ENDDO

DO i = 1, nAziAng
  jbeg = AziAng(i)%CoreRayIdxSt; jend = jbeg + AziAng(i)%nCoreRay - 1
  DO j = jbeg, jend
    SearchTree(j)%iang = i
  ENDDO
ENDDO

!Angle Reflection Mapping
DO i = 1, nAziAng/2
  j = nAziAng - i + 1
  RayRefAng(i) = j; RayRefAng(j) = i
ENDDO

!Angle Rotation Mapping
DO i = 1, nAziAng/2
  j = nAziAng / 2 + i
  RayRotAng(i) = j; RayRotAng(j) = i
ENDDO
!Generate Modular Ray
ALLOCATE(RayIdx(nCoreRay))
ALLOCATE(DirIdx(nCoreRay))
ALLOCATE(RotRay0(nCoreRay))
ALLOCATE(PhiAngBegIdx(nCoreRay, 2))
ALLOCATE(PhiAngEndIdx(nCoreRay, 2))
nRotRay = 0;
ibegidx = 1; iendidx = 1;
nSvIdx = 1;
DO i = 1, nCoreRay
  IF(SearchTree(i)%lUsed) CYCLE
  iRay = 1
  j = i
  IOSURF = SearchTree(i)%IOSurf
  ForwardBC = Core%RadSym(OSurf); BackwardBC = Core%RadSym(ISurf)
  !Initial Direction Determination(1 : Forward, 2: Backward)
  IF(ForwardBC .and. BackwardBC) THEN
    DIR =  1
    IF(.NOT. AllRef) CYCLE
    nSvIdx = nSvidx + 2
    iBegIdx = nSvIdx - 1
    iEndIdx = nSvIdx
  ELSEIF (.NOT. ForwardBC .AND. .NOT. BackwardBC) THEN
    DIR =  1
    nSvIdx = nSvidx + 1;
    iBegIdx = 1; iEndIdx = nSvIdx
  ELSEIF (ForwardBC) THEN
    DIR =  1
    nSvIdx = nSvidx + 1;
    iBegIdx = 1; iEndIdx = nSvIdx
  ELSE
    DIR = 2
    nSvIdx = nSvidx + 1;
    iBegIdx = 1; iEndIdx = nSvIdx
  ENDIF
  nRotRay = nRotRay + 1
#ifndef  RotRayDBG
  SearchTree(i)%lUsed = TRUE
#endif
  RayIdx(iray) = i; DirIdx(iray) = Dir
  OSurf = SearchTree(j)%IOSURF(mp(Dir))
  iang = SearchTree(j)%iang
  PhiAngBegIdx(nRotRay, 1) =  iBegIdx  !Forward Sweep Incoming Angular flux saved adress
  PhiAngEndIdx(nRotRay, 2) = iEndIdx    !Backwad Sweep Outgoing Angular flux saving adress
  DO WHILE(TRUE)
    !Termination Condition
    IF(.NOT. Core%RadSym(OSurf)) THEN
      nSvIdx = nSvIdx + 1;
      iBegIdx = 1; iEndIdx = nSvIdx
      iEndIdx  = nSvIdx        !Outgoing angular flux saving adress(Will not be used as incomming)

      EXIT
    ENDIF
    ixa = SearchTree(j)%ixa(mp(DIR)); iya = SearchTree(j)%iya(mp(Dir))       !Assembly Location Index at the end of Core Ray
    xref = SearchTree(j)%x(mp(DIR)); yref = SearchTree(j)%y(mp(DIR))         !Reflection
    iang = SearchTree(j)%iang
    !Find Next CoreRay
    IF(Core%lROT .AND. (OSURF .EQ. WEST .OR. OSURF .EQ. NORTH)) THEN
      irefang = RayRotAng(iang)
      CALL SetNextRaySt_ROT(OSurf, RefSurf, xref, yref, ixa, iya, Dir, SearchTree(j))
    ELSEIF(Core%lCbd) THEN
      irefang = iang
      CALL SetNextRaySt_CBD(OSurf, RefSurf, xref, yref, ixa, iya, xr, yr, nxa, nya)
    ELSE
      irefang = RayRefAng(iang)  !Reflection Angle
      RefSurf = OSurf
    ENDIF
    kbeg = AziAng(irefang)%CoreRayIdxSt; kend = kbeg + AziAng(irefang)%nCoreRay - 1
    !Core Ray Searching
    DO k = kbeg, kend
      IF(SearchTree(k)%lUsed) CYCLE
      DO m =1, 2
        iFlag = (SearchTree(k)%IOSURF(m)-RefSurf)**2
        iFlag = iFlag + (SearchTree(k)%ixa(m)-ixa)**2 + (SearchTree(k)%iya(m)-iya)**2
        !iFlag = 0
        IF(iFlag .NE. 0) CYCLE
        !Points Comparision
        Err = (xref - SearchTree(k)%X(m))**2 + (yref - SearchTree(k)%Y(m))**2
        Err = sqrt(err)
        IF(Err .GT. epsm3) CYCLE
        ! IOSURF(2)
        GOTO 100
      ENDDO
    ENDDO
    IF(k .GT. kend) THEN
      CALL Terminate('RotRayGen.f90 : Can Not Find Relfective Ray')
    ENDIF
100 CONTINUE
    DIR = m; jnew = k;
    IF(jnew .EQ. RayIdx(1)) THEN
      IBegIdx = PhiAngEndIdx(nRotRay, 2)  !Backward Sweep Outgoing == Backward Sweep Incomming
      iEndIdx = PhiAngBegIdx(nRotRay, 1)  !Forward Sweep Outgoing == Forward Sweep Incomming
      EXIT
    ENDIF
    !Saving New Ray
    iray = iray + 1
    RayIdx(iray) = jnew; DirIdx(iray) = Dir
    OSurf = SearchTree(jnew)%IOSURF(mp(Dir))
    j = jnew
#ifndef  RotRayDBG
    SearchTree(j)%lUsed = TRUE
#endif
  ENDDO
    !COPY To the temporal Ray Structure
  ALLOCATE(RotRay0(nRotRay)%RayIdx(iRay))
  ALLOCATE(RotRay0(nRotRay)%OmpRayIdx(iRay))
  ALLOCATE(RotRay0(nRotRay)%Dir(iRay))
  RotRay0(nRotRay)%nray = iray
  CALL CP_VA(RotRay0(nRotRay)%RayIdx, RayIdx(1:iray), iray)
  !CALL CP_VA(RotRay0(nRotRay)%OmpRayIdx, 0, iray)
  CALL CP_VA(RotRay0(nRotRay)%Dir, DirIdx(1:iray), iray)

  PhiAngBegIdx(nRotRay, 2) =  iBegIdx   !Backward Sweep Incoming Angular flux saved adress
  PhiAngEndIdx(nRotRay, 1) = iEndIdx    !Forward Sweep Outgoing Angular flux saving adress

  DO m = 1, iray
    jnew = RayIdx(m)
    SearchTree(jnew)%lUSed = TRUE
  ENDDO
  CONTINUE

ENDDO

DO i = 1, nCoreRay
  IF(.NOT. SearchTree(i)%lUsed) THEN
    CONTINUE
    CALL Terminate('RotRay.f90 : Incomplete set of Rotation Ray')
  ENDIF
ENDDO

ALLOCATE(RotRay(nRotRay))
RayInfo%nRotRay = nRotRay

DO i = 1, nRotRay
  RotRay(i)%nRay = RotRay0(i)%nRay
  CALL DMALLOC(RotRay(i)%RayIdx, RotRay(i)%nRay)
  CALL DMALLOC(RotRay(i)%OmpRayIdx, RotRay(i)%nRay)
  CALL DMALLOC(RotRay(i)%Dir, RotRay(i)%nRay)

  CALL CP_VA(RotRay(i)%RayIdx, RotRay0(i)%RayIdx, RotRay(i)%nRay)
  CALL CP_VA(RotRay(i)%Dir, RotRay0(i)%Dir, RotRay(i)%nRay)
  DEALLOCATE(RotRay0(i)%Dir)
  DEALLOCATE(RotRay0(i)%RayIdx)
  CONTINUE
ENDDO
RayInfo%RotRay => RotRay
RayInfo%nPhiAngSv = nSvIdx
Call Dmalloc(RayInfo%PhiangInSvIdx, nROtRay, 2)
Call Dmalloc(RayInfo%PhiangOutSvIdx, nROtRay, 2)
CALL CP_VA(RayInfo%PhiangInSvIdx, PhiAngBegIdx(1:nRotRay, 1:2), nROtRay, 2)
CALL CP_VA(RayInfo%PhiangOutSvIdx, PhiAngEndIdx(1:nRotRay, 1:2), nROtRay, 2)
!#define RotRayDBG
!#define RotRayOut
#ifdef RotRayOut
!WRITING CORERAY
i=1
DO i = 1, nRotRay
  DO iray = 1, RotRay(i)%nRAY
    j = RotRay(i)%RayIdx(iRay)
    IF(SearchTree(j)%iang .NE. 1 .AND. SearchTree(j)%iang .NE. 8) CYCLE
    CALL Write_CoreRay(j, i, RotRay(i)%Dir(iRay), FALSE)
  ENDDO
ENDDO
#endif

!--- CNJ Edit : Node Majors, GPU Acceleration
CALL Dmalloc0(RayInfo%RotRayAziList, 0, nRotRay, 1, nAziAng / 2)
DO iAng = 1, nAziAng / 2
  DO i = 1, nRotRay
    DO iRay = 1, RotRay(i)%nRay
      iCoreRay = RotRay(i)%RayIdx(iRay)
      IF (CoreRay(iCoreRay)%iAng .EQ. iAng .OR. CoreRay(iCoreRay)%iAng .EQ. nAziAng - iAng + 1) THEN
        RayInfo%RotRayAziList(0, iAng) = RayInfo%RotRayAziList(0, iAng) + 1
        RayInfo%RotRayAziList(RayInfo%RotRayAziList(0, iAng), iAng) = i
        EXIT
      ENDIF
    ENDDO
  ENDDO
ENDDO

!--- CNJ Edit : Angular Multigrid Ray Tracing
IF (.NOT. nTracerCntl%lMultigrid) THEN
  nTracerCntl%MultigridLV = 1
  nTracerCntl%gridNum = 0
ENDIF
CALL MultigridRayGen(RayInfo)

!RMEDIT_2015_08_17
DO i = 1, nRotRay 
  IF(RotRay(i)%nRay .EQ. 1) THEN
    RotRay(i)%Ang1 = CoreRay(RotRay(i)%RayIdx(1))%iang
    RotRay(i)%Ang2 = 0
  ELSE
    RotRay(i)%Ang1 = CoreRay(RotRay(i)%RayIdx(1))%iang
    RotRay(i)%Ang2 = CoreRay(RotRay(i)%RayIdx(2))%iang
  ENDIF
ENDDO
nTracerCntl%lScatBd = AllRef
!RMEDIT_2015_08_17 END

DEALLOCATE(SearchTree)
DEALLOCATE(RayIdx, DirIdx, PhiAngBegIdx, PhiAngEndIdx)

CONTAINS
SUBROUTINE SetNextRaySt_CBD(OutSurf0, RefSurf0, xref0, yref0, ixa0, iya0, xr0, yr0, nxa0, nya0)
IMPLICIT NONE
INTEGER :: OutSurf0, RefSurf0
REAL :: xref0, yref0, xr0(2), yr0(2)
INTEGER :: ixa0, iya0, nxa0, nya0
INTEGER :: CBDSurf(4)
INTEGER :: f_xref(4), f_yref(4), f_xr_1(4), f_xr_2(4), f_yr_1(4), f_yr_2(4)
CBDSurf = (/NORTH, EAST, SOUTH, WEST/)
f_xref =  (/1, 0, 1, 0/)
f_yref =  (/0, 1, 0, 1/)
f_xr_1  = (/0, 1, 0, 0/)
f_xr_2  = (/0, 0, 0, 1/)
f_yr_1  = (/0, 0, 1, 0/)
f_yr_2  = (/1, 0, 0, 0/)
RefSurf0 = CBDSurf(OutSurf0)
xref0 = xref0 * f_xref(RefSurf0) + xr0(1) * f_xr_1(RefSurf0) + xr0(2) * f_xr_2(RefSurf0)
yref0 = yref0 * f_yref(RefSurf0) + yr0(1) * f_yr_1(RefSurf0) + yr0(2) * f_yr_2(RefSurf0)
ixa0 =  ixa0 * f_xref(RefSurf0) + f_xr_1(RefSurf0) + nxa0 * f_xr_2(RefSurf0)
iya0 =  iya0 * f_yref(RefSurf0) + f_yr_1(RefSurf0) + nya0 * f_yr_2(RefSurf0)
END SUBROUTINE

SUBROUTINE SetNextRaySt_ROT(OutSurf0, RefSurf0, xref0, yref0, ixa0, iya0, Dir0, SearchTree0)
IMPLICIT NONE
INTEGER :: OutSurf0, RefSurf0
REAL :: xref0, yref0
INTEGER :: ixa0, iya0, dir0
TYPE(CoreRaySearchTree_Type) :: SearchTree0
IF(OutSurf0 .EQ. WEST) RefSurf0 = NORTH
IF(OutSurf0 .EQ. NORTH) RefSurf0 = WEST
ixa0 = SearchTree0%iya(mp(Dir0)); iya0 = SearchTree0%ixa(mp(DIR0))
yref0 = -SearchTree0%x(mp(Dir0)); xref0 = -SearchTree0%y(mp(DIR0))
END SUBROUTINE
END SUBROUTINE

!--- CNJ Edit
SUBROUTINE RayInfoMaxSize(Core, RayInfo, myzb, myze)
USE PARAM
USE TYPEDEF,        ONLY : RayInfo_type,      ModRayInfo_type,    AsyRayInfo_type,   &
                           AziAngleInfo_Type, CoreRayInfo_type,   RotRayInfo_Type,   &
                           CoreInfo_type,     Asy_Type,           Cell_Type,         &
                           Pin_Type
USE CNTL,           ONLY : nTracerCntl
USE MOC_MOD,        ONLY : nMaxRaySeg,        nMaxCellRay,        nMaxAsyRay,        &
                           nMaxCoreRay
USE GEOM,           ONLY : BaseCellInfo
IMPLICIT NONE

TYPE(CoreInfo_type) :: Core
TYPE(RayInfo_type) :: RayInfo
INTEGER :: myzb, myze

TYPE(AziAngleInfo_Type), POINTER :: AziAng(:)
TYPE(RotRayInfo_Type), POINTER  :: RotRay(:)
TYPE(CoreRayInfo_type), POINTER :: CoreRay(:)
TYPE(AsyRayInfo_type), POINTER :: AsyRay(:)

TYPE(Asy_Type),POINTER :: Asy(:)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)

INTEGER :: nRotRay, nCoreRay, nAsyRay, nAziAng, nAsySeg, nAsyCell
INTEGER :: irotray, icoray, iasyray, iceray
INTEGER :: nseg, nseg0, ncell, nasy
INTEGER :: iasy, ipin, icel, iz, itype
INTEGER :: i, j, k, l, m
INTEGER , Pointer :: nRaySegAsyRay
INTEGER :: ibcel

Asy => Core%Asy
Pin => Core%Pin
CellInfo => Core%CellInfo

AziAng => RayInfo%AziAngle
RotRay => RayInfo%RotRay
CoreRay => RayInfo%CoreRay
AsyRay => RayInfo%AsyRay

nRotRay = RayInfo%nRotRay
nAziAng = RayInfo%nAziAngle

!nAsyRay = RayInfo%nAsyRay
nMaxRaySeg = 0; nMaxCellRay = 0
nMaxAsyRay = 0; nMaxCoreRay = 0
!$OMP PARALLEL PRIVATE(iRotRay, icoray, iasyray, iasy, ipin, itype, iceray, icel, ibcel, ncell, nseg, nasy, nAsySeg, nseg0)
!$OMP DO SCHEDULE(GUIDED) REDUCTION(MAX : nMaxRaySeg, nMaxCellRay, nMaxAsyRay, nMaxCoreRay)
DO i = 1, nRotRay
  iRotRay = i
  RotRay(i)%nSeg = 0
  DO j = 1, RotRay(i)%nRay
    ncell = 0; nseg = 0; nasy = 0
    icoray = RotRay(i)%RayIdx(j)
    DO k = 1, CoreRay(icoray)%nRay
      nAsySeg = 0
      iasyray = CoreRay(icoray)%AsyRayIdx(k)
      iasy = CoreRay(icoray)%AsyIdx(k)
      nAsyCell = AsyRay(iasyray)%nCellRay
      if(iasy .eq. 0) cycle
      DO l = 1, AsyRay(iasyray)%nCellRay
        ipin = AsyRay(iasyray)%PinIdx(l)
        itype = AsyRay(iasyray)%PartialAsyRayFlag
        iceray = AsyRay(iasyray)%PinRayIdx(l)
        ipin = Asy(iasy)%GlobalPinIdx(ipin)
        nseg0 = 0
        DO iz = myzb, myze  !Find Maximum number of ray segments in the given Pin
          icel = Pin(ipin)%cell(iz)
          ibcel=CellInfo(iCel)%basecellstr
          IF (nTracerCntl%lnTIGRst) THEN
            nseg0 = max(nseg0, BaseCellInfo(ibcel)%Cellray(iceray)%nSeg)
          ELSE
            nseg0 = max(nseg0, CellInfo(ibcel)%Cellray(iceray)%nSeg)
          ENDIF
        ENDDO
        nseg = nseg + nseg0
        nAsySeg = nAsySeg + nseg0
      ENDDO !End of Cell Ray Sweep
      ncell = ncell + AsyRay(iasyray)%nCellRay  !Add Number of Cell Ray
    ENDDO !End of Asy Ray
    nasy = nasy + CoreRay(icoray)%nRay
    nMaxRaySeg = Max(nSeg, nMaxRaySeg)
    nMaxCellRay = Max(nCell, nMaxCellRay)
    nMaxAsyRay = Max(nAsy, nMaxAsyRay)
    RotRay(i)%nseg = RotRay(i)%nseg + nseg
  ENDDO !End of Coreray

  nMaxCoreRay = Max(RotRay(i)%nRay, nMaxCoreRay)
ENDDO !End of Rotray
!$OMP END DO
!$OMP END PARALLEL

CONTINUE
END SUBROUTINE

!--- CNJ Edit : Domain Decomposition
SUBROUTINE RayInfo4CmfdGen(RayInfo, RayInfo4Cmfd, Core)
USE PARAM
USE TYPEDEF,        ONLY : RayInfo_type,      RayInfo4Cmfd_Type,  CoreInfo_type,        &
                           ModRayInfo_type,   AsyRayInfo_type,    DcmpAsyRayInfo_Type,  &
                           AziAngleInfo_Type, CoreRayInfo_type,   RotRayInfo_Type,      &
                           Asy_Type,           Cell_Type,         Pin_Type
USE CNTL,           ONLY : nTracerCntl
USE ALLOCS
IMPLICIT NONE

TYPE(CoreInfo_type) :: Core
TYPE(RayInfo_type) :: RayInfo
TYPE(RayInfo4Cmfd_type), POINTER :: RayInfo4Cmfd

TYPE(AziAngleInfo_Type), POINTER :: AziAng(:)
TYPE(RotRayInfo_Type), POINTER  :: RotRay(:)
TYPE(CoreRayInfo_type), POINTER :: CoreRay(:)
TYPE(AsyRayInfo_type), POINTER :: AsyRay(:)
TYPE(DcmpAsyRayInfo_Type), POINTER :: DcmpAsyRay(:, :)

TYPE(Asy_Type),POINTER :: Asy(:)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)

INTEGER :: nAsy
INTEGER :: nRotRay, nCoreRay, nModRay, nPinRay, nAziAng
INTEGER :: irotray, icoray, iasyray, iceray
INTEGER :: iAsy, ipin, itype, idir, isv
INTEGER :: i, j, k, l, m
INTEGER :: iAsyRayBeg, iAsyRayEnd, iAsyBeg, iAsyEnd

Asy => Core%Asy
Pin => Core%Pin
CellInfo => Core%CellInfo

AziAng => RayInfo%AziAngle
RotRay => RayInfo%RotRay
CoreRay => RayInfo%CoreRay
AsyRay => RayInfo%AsyRay
DcmpAsyRay => RayInfo%DcmpAsyRay

nAsy = Core%nxya

nRotRay = RayInfo%nRotRay
nModRay = RayInfo%nModRay
nAziAng = RayInfo%nAziAngle

nRotRay = RayInfo%nRotRay
nAziAng = RayInfo%nAziAngle

ALLOCATE(RayInfo4CMFD)

RayInfo4CMFD%nRotRay = nRotRay
RayInfo4CMFD%nCoreRay = RayInfo%nCoreRay
RayInfo4CMFD%nPolAngle = RayInfo%nPolarAngle

CALL DMALLOC(RayInfo4CMFD%RotRayInOutCell, nRotRay, 2)

DO i = 1, nRotRay
  !1, RotRay(i)%nRay
  nCoreRay = RotRay(i)%nRay
  !Begining Cell Info
  isv = RayInfo%PhiangInSvIdx(i, 1); j = 1
  icoray = RotRay(i)%RayIdx(j); idir = RotRay(i)%DIR(j)
  k = 1; IF(idir .eq. 2) k =  CoreRay(icoray)%nRay
  iasyray = CoreRay(icoray)%AsyRayIdx(k); iasy = CoreRay(icoray)%AsyIdx(k)
  IF(iasy .NE. 0) THEN
    nPinRay = AsyRay(iAsyRay)%nCellRay; itype = Asy(iasy)%PartialAsyFlag
    l = 1; IF(idir .eq. 2) l = nPinRay
    ipin = AsyRay(iAsyRay)%PinIdx(l); ipin = Asy(iAsy)%GlobalPinIdx(ipin)
    IF(isv .NE. 1) RayInfo4CMfd%RotRayInOutCell(i, 1) = ipin
  ENDIF
  !End Cell Info
  isv = RayInfo%PhiangInSvIdx(i, 2); j = nCoreRay
  icoray = RotRay(i)%RayIdx(j); idir = RotRay(i)%DIR(j)
  k = CoreRay(icoray)%nRay; IF(idir .eq. 2) k =  1
  iasyray = CoreRay(icoray)%AsyRayIdx(k); iasy = CoreRay(icoray)%AsyIdx(k)
  IF(iasy .NE.0) THEN
    nPinRay = AsyRay(iAsyRay)%nCellRay; itype = Asy(iasy)%PartialAsyFlag
    l = nPinRay; IF(idir .eq. 2) l = 1
    ipin = AsyRay(iAsyRay)%PinIdx(l); ipin = Asy(iAsy)%GlobalPinIdx(ipin)
    IF(isv .NE. 1) RayInfo4CMfd%RotRayInOutCell(i, 2) = ipin
  ENDIF
ENDDO

IF (nTracerCntl%lDomainDcmp) THEN
  CALL Dmalloc(RayInfo4CMFD%DcmpAsyRayInCell, 2, nModRay, nAsy)
  CALL Dmalloc(RayInfo4CMFD%DcmpAsyRayInSurf, 2, nModRay, nAsy)
  RayInfo4CMFD%DcmpAsyRayCount => RayInfo%DcmpAsyRayCount
  DO i = 1, nAsy
    DO j = 1, RayInfo%DcmpAsyRayCount(i)
      iAsyRayBeg = DcmpAsyRay(j, i)%AsyRayList(1)
      iAsyRayEnd = DcmpAsyRay(j, i)%AsyRayList(DcmpAsyRay(j, i)%nAsyRay)
      IF (DcmpAsyRay(j, i)%DirList(1) .EQ. 1) THEN
        RayInfo4CMFD%DcmpAsyRayInSurf(1, j, i) = AsyRay(iAsyRayBeg)%InOutSurf(1)
        RayInfo4CMFD%DcmpAsyRayInCell(1, j, i) = AsyRay(iAsyRayBeg)%PinIdx(1)
      ELSE
        RayInfo4CMFD%DcmpAsyRayInSurf(1, j, i) = AsyRay(iAsyRayBeg)%InOutSurf(2)
        RayInfo4CMFD%DcmpAsyRayInCell(1, j, i) = AsyRay(iAsyRayBeg)%PinIdx(AsyRay(iAsyRayBeg)%nCellRay)
      ENDIF
      IF (DcmpAsyRay(j, i)%DirList(DcmpAsyRay(j, i)%nAsyRay) .EQ. 1) THEN
        RayInfo4CMFD%DcmpAsyRayInSurf(2, j, i) = AsyRay(iAsyRayEnd)%InOutSurf(2)
        RayInfo4CMFD%DcmpAsyRayInCell(2, j, i) = AsyRay(iAsyRayEnd)%PinIdx(AsyRay(iAsyRayEnd)%nCellRay)
      ELSE
        RayInfo4CMFD%DcmpAsyRayInSurf(2, j, i) = AsyRay(iAsyRayEnd)%InOutSurf(1)
        RayInfo4CMFD%DcmpAsyRayInCell(2, j, i) = AsyRay(iAsyRayEnd)%PinIdx(1)
      ENDIF
      RayInfo4CMFD%DcmpAsyRayInCell(1, j, i) = Asy(i)%GlobalPinIdx(RayInfo4CMFD%DcmpAsyRayInCell(1, j, i))
      RayInfo4CMFD%DcmpAsyRayInCell(2, j, i) = Asy(i)%GlobalPinIdx(RayInfo4CMFD%DcmpAsyRayInCell(2, j, i))
    ENDDO
  ENDDO
ENDIF

!RayInfo%RayInfo4Cmfd => RayInfo4CMfd
RayInfo4CMFD%PhiAngInSvIdx => RayInfo%PhiAngInSvIdx
RayInfo4CMFD%PolarAngle => RayInfo%PolarAngle

END SUBROUTINE

!--- CNJ Edit : Angular Multigrid Ray Tracing
SUBROUTINE MultigridRayGen(RayInfo)
USE TYPEDEF,        ONLY : RayInfo_Type,    AziAngleInfo_Type,  MultigridInfo_Type
USE CNTL,           ONLY : nTracerCntl
USE MOC_MOD,        ONLY : AziMap
USE IOUTIL,         ONLY : terminate
IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo

TYPE(AziAngleInfo_Type), POINTER :: AziAngle(:)
TYPE(MultigridInfo_Type), POINTER :: MultigridInfo
REAL, POINTER :: wtang(:, :), wtsurf(:, :, :), mwt(:, :, :), mwt2(:, :, :), Comp(:, :, :)
REAL :: WeightSet(4, 4), SinvSet(4, 4)
REAL :: wt(2), wttemp, wtsin2, wtcos, wtpolar
REAL :: wtAzi(400), wtPol(4)
REAL :: psinv(4), pcosv(4)
INTEGER :: lvr(10, 0 : 10, 2)
INTEGER :: i, j, k, ilv, jlv, ipol, iazi, AziIdx
INTEGER :: nlv, nAzi, nPolar
REAL :: Inclv

DATA WeightSet / 1.0000000, 0.0000000, 0.0000000, 0.0000000, &
                 0.2128540, 0.7871460, 0.0000000, 0.0000000, &
                 0.0462330, 0.2836190, 0.6701480, 0.0000000, &
                 0.0116656, 0.0855345, 0.3079014, 0.5948985  /                
DATA SinvSet   / 0.7071068, 0.0000000, 0.0000000, 0.0000000, &
                 0.3639000, 0.8999000, 0.0000000, 0.0000000, &
                 0.1666480, 0.5377070, 0.9329540, 0.0000000, &
                 0.0834435, 0.2908650, 0.6328485, 0.9492286  /
               
ALLOCATE(RayInfo%MultigridInfo(nTracerCntl%MultigridLV))

AziAngle => RayInfo%AziAngle
nlv = nTracerCntl%MultigridLV
nAzi = RayInfo%nAziAngle
nPolar = RayInfo%nPolarAngle

IF (mod(nAzi / 2, nlv) .NE. 0) CALL terminate("Multigrid Setup Fail : # of Azimuthal Angles Cannot be Divided by Grid Level")

DO iazi = 1, nAzi / 2
  AziMap(iazi, 1) = 1
  AziMap(iazi, 2) = 2
  AziMap(nAzi - iazi + 1, 1) = 3
  AziMap(nAzi - iazi + 1, 2) = 4
ENDDO
  
IF (nlv .EQ. 1) THEN
  lvr(1, 1, :) = 0
ELSE
  DO ilv = 1, nlv
    Inclv = dble(nlv) / dble(ilv)
    lvr(ilv, 0, 2) = -1
    DO jlv = 1, ilv
      lvr(ilv, jlv, 1) = lvr(ilv, jlv - 1, 2) + 1
      lvr(ilv, jlv, 2) = nint(Inclv * jlv) - 1
    ENDDO
  ENDDO
ENDIF

DO ilv = 1, nlv

  MultigridInfo => RayInfo%MultigridInfo(ilv)
  MultigridInfo%nAzi = nAzi * ilv / nlv
  MultigridInfo%nPolar = nPolar ! max(2, nPolar * ilv / nlv)
  ALLOCATE(MultigridInfo%AziList(MultigridInfo%nAzi))
  wtAzi = 0.0; iazi = 0
  DO i = 1, nAzi / 2, nlv
    DO jlv = 1, ilv
      iazi = iazi + 1
      wt = 0.0
      DO k = i + lvr(ilv, jlv, 1), i + lvr(ilv, jlv, 2)
        wt(1) = wt(1) + AziAngle(k)%weight
        wt(2) = wt(2) + AziAngle(nAzi - k + 1)%weight
      ENDDO
      AziIdx = i + lvr(ilv, jlv, 1)
      wtAzi(AziIdx) = wt(1)
      wtAzi(nAzi - AziIdx + 1) = wt(2)
      MultigridInfo%AziList(iazi) = AziIdx
      MultigridInfo%AziList(MultigridInfo%nAzi - iazi + 1) = nAzi - AziIdx + 1
    ENDDO
  ENDDO
  
  wtPol = WeightSet(:, MultigridInfo%nPolar)
  psinv = SinvSet(:, MultigridInfo%nPolar)
  pcosv = sqrt(1.0 - psinv ** 2)
  
  ALLOCATE(MultigridInfo%EXPA(MultigridInfo%nPolar, -40000 : 0))
  ALLOCATE(MultigridInfo%EXPB(MultigridInfo%nPolar, -40000 : 0))
  CALL ExpTabulation(MultigridInfo%nPolar, psinv, MultigridInfo%EXPA, MultigridInfo%EXPB)
  
  ALLOCATE(MultigridInfo%wtang(nPolar, nAzi)); wtang => MultigridInfo%wtang
  ALLOCATE(MultigridInfo%wtsurf(nPolar, nAzi, 4)); wtsurf => MultigridInfo%wtsurf
  ALLOCATE(MultigridInfo%Comp(9, nPolar, nAzi)); Comp => MultigridInfo%Comp
  ALLOCATE(MultigridInfo%mwt(9, nPolar, nAzi)); mwt => MultigridInfo%mwt
  ALLOCATE(MultigridInfo%mwt2(9, nPolar, nAzi)); mwt2 => MultigridInfo%mwt2
  DO ipol = 1, MultigridInfo%nPolar
    wttemp = wtPol(ipol) * psinv(ipol)
    DO iazi = 1, MultigridInfo%nAzi
      AziIdx = MultigridInfo%AziList(iazi)
      wtang(ipol, AziIdx) = wttemp * wtAzi(AziIdx) * AziAngle(AziIdx)%del
      wtsurf(ipol, AziIdx, 1) = wtPol(ipol) * wtAzi(AziIdx) * AziAngle(AziIdx)%del / abs(AziAngle(AziIdx)%sinv)
      wtsurf(ipol, AziIdx, 3) = wtPol(ipol) * wtAzi(AziIdx) * AziAngle(AziIdx)%del / abs(AziAngle(AziIdx)%sinv)
      wtsurf(ipol, AziIdx, 2) = wtPol(ipol) * wtAzi(AziIdx) * AziAngle(AziIdx)%del / abs(AziAngle(AziIdx)%cosv)
      wtsurf(ipol, AziIdx, 4) = wtPol(ipol) * wtAzi(AziIdx) * AziAngle(AziIdx)%del / abs(AziAngle(AziIdx)%cosv)
    ENDDO
  ENDDO
  IF (nTracerCntl%ScatOd .GE. 1) THEN
    DO ipol = 1, MultigridInfo%nPolar
      wttemp = psinv(ipol)
      DO iazi = 1, MultigridInfo%nAzi
        AziIdx = MultigridInfo%AziList(iazi)
        Comp(1, ipol, AziIdx) = wttemp * AziAngle(AziIdx)%cosv
        Comp(2, ipol, AziIdx) = wttemp * AziAngle(AziIdx)%sinv
        mwt(1:2, ipol, AziIdx) = Comp(1:2, ipol, AziIdx) * wtang(ipol, AziIdx)
        
        mwt2(1:2, ipol, AziIdx) = -mwt(1:2, ipol, AziIdx)
      ENDDO
    ENDDO
  ENDIF
  IF (nTracerCntl%ScatOd .GE. 2) THEN
    DO ipol = 1, MultigridInfo%nPolar
      wttemp = psinv(ipol)
      wtsin2 = psinv(ipol) * psinv(ipol)
      wtcos = pcosv(ipol)
      wtpolar = 1.5 * pcosv(ipol) * pcosv(ipol) - 0.5
      DO iazi = 1, MultigridInfo%nAzi
        AziIdx = MultigridInfo%AziList(iazi)
        Comp(3, ipol, AziIdx) = wtpolar
        Comp(4, ipol, AziIdx) = wtsin2 * (1.0 - 2.0 * AziAngle(AziIdx)%sinv * AziAngle(AziIdx)%sinv)
        Comp(5, ipol, AziIdx) = wtsin2 * (2.0 * AziAngle(AziIdx)%sinv * AziAngle(AziIdx)%cosv)
        mwt(3, ipol, AziIdx) = Comp(3, ipol, AziIdx) * wtang(ipol, AziIdx)
        mwt(4:5, ipol, AziIdx) = 0.75 * Comp(4:5, ipol, AziIdx) * wtang(ipol, AziIdx)
        
        mwt2(3:5, ipol, AziIdx) = mwt(3:5, ipol, AziIdx)
      ENDDO
    ENDDO
  ENDIF
  IF (nTracerCntl%ScatOd .EQ. 3) THEN
    DO ipol = 1, MultigridInfo%nPolar
      wttemp = psinv(ipol)
      DO iazi = 1, MultigridInfo%nAzi
        AziIdx = MultigridInfo%AziList(iazi)
        Comp(6, ipol, AziIdx) = (5.0 * pcosv(ipol) * pcosv(ipol) - 1.0) * wttemp * AziAngle(AziIdx)%cosv
        Comp(7, ipol, AziIdx) = (5.0 * pcosv(ipol) * pcosv(ipol) - 1.0) * wttemp * AziAngle(AziIdx)%sinv
        Comp(8, ipol, AziIdx) = (wttemp ** 3.0) * (4.0 * (AziAngle(AziIdx)%cosv ** 3.0) - 3.0 * AziAngle(AziIdx)%cosv)
        Comp(9, ipol, AziIdx) = (wttemp ** 3.0) * (- 4.0 * (AziAngle(AziIdx)%sinv ** 3.0) + 3.0 * AziAngle(AziIdx)%sinv)
        mwt(6:7, ipol, AziIdx) = 0.375 * Comp(6:7, ipol, AziIdx) * wtang(ipol, AziIdx)
        mwt(8:9, ipol, AziIdx) = 0.625 * Comp(8:9, ipol, AziIdx) * wtang(ipol, AziIdx)
        
        mwt2(6:9, ipol, AziIdx) = -mwt(6:9, ipol, AziIdx)
      ENDDO
    ENDDO
  ENDIF
  
ENDDO

CONTAINS

SUBROUTINE ExpTabulation(npr, sinvpol, expa, expb)
USE PARAM
IMPLICIT NONE

INTEGER :: npr
INTEGER :: ipr,i
REAL :: ttt, xval1,xval2, dx
REAL :: expa(npr,-40000:0), expb(npr,-40000:0)
REAL :: yval1(npr), yval2(npr), sinvpol(4), rsinvpol(npr)

DO ipr = 1, npr
  rsinvpol(ipr) = one/sinvpol(ipr)
ENDDO

dx=epsm3
xval1=-40.0_8
do ipr=1,npr
  expa(ipr,-40000)=0
  expb(ipr,-40000)=1
  ttt=rsinvpol(ipr)*xval1
  if(ttt.lt.-200.0_8) then
    yval1(ipr)=1
  else
    yval1(ipr)=1-exp(ttt)
  endif  
enddo
do i=-39999,0
  xval2=xval1+dx
  do ipr=1,npr
    ttt=rsinvpol(ipr)*xval2
    if(ttt.lt.-100) then
      yval2(ipr)=1
    else
      yval2(ipr)=1-exp(ttt)
    endif
    expa(ipr,i)=1000._8*(yval2(ipr)-yval1(ipr))
    expb(ipr,i)=yval1(ipr)-expa(ipr,i)*xval1
    yval1(ipr)=yval2(ipr)
  enddo
  xval1=xval2
enddo

DO ipr = 1, npr
  DO i = -39999,0
    expa(ipr,i) = dx * expa(ipr,i)
  ENDDO
ENDDO 

END SUBROUTINE

END SUBROUTINE