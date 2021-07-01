#include <defines.h>
SUBROUTINE DcmpRayGen(Core, RayInfo, DcmpAsyRay)

USE ALLOCS
USE PARAM,   ONLY : TRUE, BACKWARD, FORWARD, RED, BLACK
USE TYPEDEF, ONLY : RayInfo_type, RotRayInfo_type, CoreRayInfo_type, AsyRayInfo_Type, CoreInfo_type, DcmpAsyRayInfo_Type, Asy_Type, Pin_Type, Cell_Type
USE PE_Mod,  ONLY : PE
USE MOC_MOD, ONLY : DcmpColorAsy

IMPLICIT NONE

TYPE (CoreInfo_type) :: Core
TYPE (RayInfo_type)  :: RayInfo

TYPE (DcmpAsyRayInfo_type), POINTER, DIMENSION(:,:) :: DcmpAsyRay

TYPE (RotRayInfo_Type),  POINTER, DIMENSION(:) :: RotRay
TYPE (CoreRayInfo_Type), POINTER, DIMENSION(:) :: CoreRay
TYPE (AsyRayInfo_Type),  POINTER, DIMENSION(:) :: AsyRay
TYPE (Asy_Type),         POINTER, DIMENSION(:) :: Asy
TYPE (Pin_Type),         POINTER, DIMENSION(:) :: Pin
TYPE (Cell_Type),        POINTER, DIMENSION(:) :: Cell

INTEGER, POINTER, DIMENSION(:)       :: DcmpAsyRayCount, AsyRayList, DirList, AziList
INTEGER, POINTER, DIMENSION(:,:,:)   :: DcmpAsyAziList
INTEGER, POINTER, DIMENSION(:,:,:,:) :: DcmpAsyLinkInfo

INTEGER :: nRotRay, nCoreRay, nAsyRay, nModRay, nDummyRay, nAsy, nMaxAziModRay, nMaxCellRay, nMaxRaySeg, nPinRay, nRaySeg, nRaySeg0, nAziAngle
INTEGER :: iRotRay, iCoreRay, jCoreRay, iAsyRay, jAsyRay, icelray, jcelray, iAzi, iDir, iz, iasy, icel, ibcel, ipin, iCnt, icolor
INTEGER :: AsyRayBeg, AsyRayEnd, AsyRayInc, myzb, myze, prvAsy, prvCnt, nRef
! ----------------------------------------------------

nAsy  = Core%nxya
Asy  => Core%Asy
Pin  => Core%Pin
Cell => Core%CellInfo

nModRay   = RayInfo%nModRay
nRotRay   = RayInfo%nRotRay
RotRay   => RayInfo%RotRay
CoreRay  => RayInfo%CoreRay
AsyRay   => RayInfo%AsyRay
nAziAngle = RayInfo%nAziAngle

myzb = PE%myzb
myze = PE%myze

ALLOCATE (DcmpAsyRay(nModRay, nAsy))

CALL dmalloc (RayInfo%DcmpAsyRayCount, nAsy)
CALL dmalloc (RayInfo%DcmpAsyLinkInfo, 2, 2, nModRay, nAsy)
CALL dmalloc0(RayInfo%DcmpAsyAziList,  0, nModRay, 1, nAziAngle/2, 1, nAsy)

DcmpAsyRayCount => RayInfo%DcmpAsyRayCount
DcmpAsyAziList  => RayInfo%DcmpAsyAziList
DcmpAsyLinkInfo => RayInfo%DcmpAsyLinkInfo
! ----------------------------------------------------
DO iRotRay = 1, nRotRay
  nCoreRay = RotRay(iRotRay)%nRay
  prvAsy   = 0
  prvCnt   = 0
  
  DO iCoreRay = 1, nCoreRay
    jCoreRay  = RotRay(iRotRay)%RayIdx(iCoreRay)
    nAsyRay   = CoreRay(jCoreRay)%nRay
    nDummyRay = 0
    iDir      = RotRay(iRotRay)%Dir(iCoreRay)
    
    IF (iDir .EQ. BACKWARD) THEN
      AsyRayBeg = nAsyRay; AsyRayEnd = 1; AsyRayInc = -1
    ELSE
      AsyRayBeg = 1; AsyRayEnd = nAsyRay; AsyRayInc = 1
    END IF
    
    DO iAsyRay = AsyRayBeg, AsyRayEnd, AsyRayInc
      jAsyRay = CoreRay(jCoreRay)%AsyRayIdx(iAsyRay)
      iAsy    = CoreRay(jCoreRay)%AsyIdx(iAsyRay)
      
      IF (iAsy .EQ. 0) THEN
        nDummyRay = nDummyRay + AsyRayInc
        
        CYCLE
      END IF
      
      iCnt = DcmpAsyRayCount(iAsy)
      
      IF (iAsy .EQ. prvAsy) THEN
        nRef    = nRef + 1
        nRaySeg = 0
        
        DO icelray = 1, AsyRay(jAsyRay)%nCellRay
          ipin     = AsyRay(jAsyRay)%PinIdx(icelray)
          jcelray  = AsyRay(jAsyRay)%PinRayIdx(icelray)
          ipin     = Asy(iAsy)%GlobalPinIdx(ipin)
          
          nRaySeg0 = 0
          DO iz = myzb, myze
            icel     = Pin(ipin)%cell(iz)
            ibcel    = Cell(iCel)%basecellstr
            nRaySeg0 = max(nRaySeg0, Cell(ibcel)%Cellray(jcelray)%nSeg)
          END DO
          
          nRaySeg = nRaySeg + nRaySeg0
        END DO
        
        nMaxRaySeg  = max(nMaxRaySeg,  nRaySeg) ! # of Seg. in Asy. Ray
        nMaxCellRay = max(nMaxCellRay, AsyRay(jAsyRay)%nCellRay)
        
        AsyRayList => DcmpAsyRay(iCnt, iAsy)%AsyRayList
        DirList    => DcmpAsyRay(iCnt, iAsy)%DirList
        AziList    => DcmpAsyRay(iCnt, iAsy)%AziList
        
        DcmpAsyRay(iCnt, iAsy)%AsyRayList => NULL()
        DcmpAsyRay(iCnt, iAsy)%DirList    => NULL()
        DcmpAsyRay(iCnt, iAsy)%AziList    => NULL()
        
        CALL dmalloc(DcmpAsyRay(iCnt, iAsy)%AsyRayList, nRef)
        CALL dmalloc(DcmpAsyRay(iCnt, iAsy)%DirList,    nRef)
        CALL dmalloc(DcmpAsyRay(iCnt, iAsy)%AziList,    nRef)
        
        DcmpAsyRay(iCnt, iAsy)%AsyRayList(1:nRef-1) = AsyRayList(:)
        DcmpAsyRay(iCnt, iAsy)%DirList   (1:nRef-1) = DirList(:)
        DcmpAsyRay(iCnt, iAsy)%AziList   (1:nRef-1) = AziList(:)
        
        DEALLOCATE (AsyRayList, DirList, AziList)
        NULLIFY    (AsyRayList, DirList, AziList)
      ELSE
        nRef    = 1
        nRaySeg = 0
        
        DcmpAsyRayCount(iAsy) = DcmpAsyRayCount(iAsy) + 1
        
        iCnt = iCnt + 1
        
        DO icelray = 1, AsyRay(jAsyRay)%nCellRay
          ipin     = AsyRay(jAsyRay)%PinIdx(icelray)
          jcelray  = AsyRay(jAsyRay)%PinRayIdx(icelray)
          ipin     = Asy(iAsy)%GlobalPinIdx(ipin)
          
          nRaySeg0 = 0
          DO iz = myzb, myze
            icel     = Pin(ipin)%cell(iz)
            ibcel    = Cell(iCel)%basecellstr
            nRaySeg0 = max(nRaySeg0, Cell(ibcel)%Cellray(jcelray)%nSeg)
          END DO
          
          nRaySeg = nRaySeg + nRaySeg0
        END DO
        
        nMaxRaySeg  = nRaySeg
        nMaxCellRay = AsyRay(jAsyRay)%nCellRay
        
        CALL dmalloc(DcmpAsyRay(iCnt, iAsy)%AsyRayList, nRef)
        CALL dmalloc(DcmpAsyRay(iCnt, iAsy)%DirList,    nRef)
        CALL dmalloc(DcmpAsyRay(iCnt, iAsy)%AziList,    nRef)
        
        IF (prvAsy.NE.0 .AND. prvCnt.NE.0) THEN
          DcmpAsyLinkInfo(1, FORWARD, prvCnt, prvAsy) = iAsy
          DcmpAsyLinkInfo(2, FORWARD, prvCnt, prvAsy) = iCnt
          
          DcmpAsyLinkInfo(1, BACKWARD, iCnt, iAsy) = prvAsy
          DcmpAsyLinkInfo(2, BACKWARD, iCnt, iAsy) = prvCnt
        END IF
      END IF
            
      !DcmpAsyRay(iCnt, iAsy)%nMaxRaySeg  = nMaxRaySeg, Useless
      DcmpAsyRay(iCnt, iAsy)%nMaxCellRay = nMaxCellRay
      
      DcmpAsyRay(iCnt, iAsy)%nAsyRay = nRef
      DcmpAsyRay(iCnt, iAsy)%iRotRay = iRotRay
      DcmpAsyRay(iCnt, iAsy)%iAsy    = iAsy
      DcmpAsyRay(iCnt, iAsy)%iRay    = iCnt
      
      DcmpAsyRay(iCnt, iAsy)%AsyRayList(nRef) = jAsyRay
      DcmpAsyRay(iCnt, iAsy)%AziList   (nRef) = CoreRay(jCoreRay)%iAng
      DcmpAsyRay(iCnt, iAsy)%DirList   (nRef) = iDir
      
      IF (iCoreRay.EQ.1 .AND. iAsyRay.EQ.(AsyRayBeg + nDummyRay)) THEN
        DcmpAsyRay(iCnt, iAsy)%lRotRayBeg(FORWARD) = TRUE
        
        DcmpAsyLinkInfo(1, BACKWARD, iCnt, iAsy) = 0
        DcmpAsyLinkInfo(2, BACKWARD, iCnt, iAsy) = RayInfo%PhiAngOutSvIdx(iRotRay, BACKWARD)
      END IF
      
      prvAsy = iAsy
      prvCnt = iCnt
    END DO
  END DO
  
  IF (iAsy.NE.0 .AND. iCnt.NE.0) THEN
    DcmpAsyRay(iCnt, iAsy)%lRotRayBeg(BACKWARD) = TRUE
    
    DcmpAsyLinkInfo(1, FORWARD, iCnt, iAsy) = 0
    DcmpAsyLinkInfo(2, FORWARD, iCnt, iAsy) = RayInfo%PhiAngOutSvIdx(iRotRay, FORWARD)
  END IF
END DO
! ----------------------------------------------------
CALL dmalloc0(DcmpColorAsy, 0, nAsy, 1, 2)

DO iAsy = 1, nAsy
  IF (mod(Asy(iAsy)%ixa, 2) .EQ. mod(Asy(iAsy)%iya, 2)) THEN
    icolor = RED
  ELSE
    icolor = BLACK
  END IF
  
  Asy(iAsy)%color = icolor
  
  DcmpColorAsy(0, icolor) = DcmpColorAsy(0, icolor) + 1
  
  DcmpColorAsy(DcmpColorAsy(0, icolor), icolor) = iAsy
    
  DO iCnt = 1, DcmpAsyRayCount(iAsy)
    DO iAzi = 1, nAziAngle / 2
      IF (ANY(DcmpAsyRay(iCnt, iAsy)%AziList.NE.iAzi .OR. DcmpAsyRay(iCnt, iAsy)%AziList.NE.nAziAngle - iAzi + 1)) CYCLE
      
      DcmpAsyAziList(0, iAzi, iAsy) = DcmpAsyAziList(0, iAzi, iAsy) + 1
      
      DcmpAsyAziList(DcmpAsyAziList(0, iAzi, iAsy), iAzi, iAsy) = iCnt
      
      EXIT
    END DO
  END DO
END DO

RayInfo%DcmpAsyRay => DcmpAsyRay
! ----------------------------------------------------

END SUBROUTINE DcmpRayGen