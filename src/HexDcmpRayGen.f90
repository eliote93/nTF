#include <defines.h>
SUBROUTINE HexDcmpRayGen(Core, RayInfo, DcmpAsyRay)

USE ALLOCS
USE PARAM,   ONLY : TRUE, BACKWARD, FORWARD, RED, BLACK, GREEN
USE TYPEDEF, ONLY : RayInfo_type, CoreInfo_type, DcmpAsyRayInfo_Type, Pin_Type
USE PE_Mod,  ONLY : PE
USE MOC_MOD, ONLY : DcmpAsyClr, DcmpAziRay
USE HexData, ONLY : hAsy, hLgc
USE HexType, ONLY : Type_HexAsyRay, Type_HexCelRay, Type_HexCoreRay, Type_HexRotRay
USE HexData, ONLY : haRay, hcRay, hRotRay, hAsyTypInfo, hLgc

IMPLICIT NONE

TYPE (CoreInfo_type) :: Core
TYPE (RayInfo_type)  :: RayInfo

TYPE (DcmpAsyRayInfo_type), POINTER, DIMENSION(:,:) :: DcmpAsyRay
! ----------------------------------------------------
INTEGER, POINTER, DIMENSION(:)       :: DcmpAsyRayCount, AsyRayList, DirList, AziList
INTEGER, POINTER, DIMENSION(:,:,:,:) :: DcmpAsyLinkInfo

INTEGER :: nRotRay, nCoreRay, nAsyRay, nModRay, nDummyRay, nAsy, nMaxCellRay
INTEGER :: iRotRay, iCoreRay, jCoreRay, iAsyRay, jAsyRay, icelray, iAzi, iz, iasy, ipin, iCnt, iAsyTyp, iGeoTyp, icBss, iDir, itmp, iClr
INTEGER :: AsyRayBeg, AsyRayEnd, AsyRayInc, myzb, myze, prvAsy, prvCnt, nRef

INTEGER :: PinSt, PinEd, FsrSt, FsrEd, maxNumPin, maxNumFsr ! DEBUG

TYPE (Pin_Type), POINTER, DIMENSION(:) :: Pin

TYPE (Type_HexAsyRay), POINTER :: haRay_Loc
TYPE (Type_HexCelRay), POINTER :: CelRay_Loc
TYPE (Type_HexRotRay), POINTER :: hRotRay_Loc
! ----------------------------------------------------

nAsy = Core%nxya
Pin => Core%Pin

nModRay = RayInfo%nModRay
nRotRay = RayInfo%nRotRay

myzb = PE%myzb
myze = PE%myze

ALLOCATE (DcmpAsyRay(nModRay, nAsy))

CALL dmalloc (RayInfo%DcmpAsyRayCount, nAsy) ! # of mRay (Considering reflections as one)
CALL dmalloc (RayInfo%DcmpAsyLinkInfo, 2, 2, nModRay, nAsy) ! (iAsy/iCnt, iDir, ~)

DcmpAsyRayCount => RayInfo%DcmpAsyRayCount
DcmpAsyLinkInfo => RayInfo%DcmpAsyLinkInfo
! ----------------------------------------------------
DO iRotRay = 1, nRotRay
  hRotRay_Loc => hRotRay(iRotRay)
  nCoreRay     = hRotRay_Loc%ncRay
  
  prvAsy = 0
  prvCnt = 0
  
  DO iCoreRay = 1, nCoreRay
    jCoreRay  = hRotRay_Loc%cRayIdx(iCoreRay)
    nAsyRay   = hcRay(abs(jCoreRay))%nmRay
    nDummyRay = 0
    
    IF (jCoreRay .LT. 0) THEN
      AsyRayBeg = nAsyRay; AsyRayEnd = 1; AsyRayInc = -1; iDir = BACKWARD
    ELSE
      AsyRayBeg = 1; AsyRayEnd = nAsyRay; AsyRayInc = 1;  iDir = FORWARD
    END IF
    
    DO iAsyRay = AsyRayBeg, AsyRayEnd, AsyRayInc
      jAsyRay = hcRay(abs(jCoreRay))%mRayIdx(iAsyRay)
      iAsy    = hcRay(abs(jCoreRay))%AsyIdx (iAsyRay)
      
      IF (iAsy .EQ. 0) THEN
        nDummyRay = nDummyRay + AsyRayInc
        
        CYCLE
      END IF
      
      iAsyTyp = hAsy(iAsy)%AsyTyp
      iGeoTyp = hAsy(iAsy)%GeoTyp
      icBss   = hAsyTypInfo(iAsyTyp)%iBss
      iCnt    = DcmpAsyRayCount(iAsy)
      
      haRay_Loc => haRay(iGeoTyp, icBss, jAsyRay)
      ! --------------------------------------------------
      IF (iAsy .EQ. prvAsy) THEN ! Case : Same Asy.
        nMaxCellRay = max(nMaxCellRay, haRay_Loc%nhpRay)
        
        nRef = nRef + 1
        
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
      ELSE ! Case : New Asy.
        nMaxCellRay = haRay_Loc%nhpRay
        
        nRef = 1
        
        DcmpAsyRayCount(iAsy) = DcmpAsyRayCount(iAsy) + 1
        iCnt = iCnt + 1
        
        CALL dmalloc(DcmpAsyRay(iCnt, iAsy)%AsyRayList, 1)
        CALL dmalloc(DcmpAsyRay(iCnt, iAsy)%DirList,    1)
        CALL dmalloc(DcmpAsyRay(iCnt, iAsy)%AziList,    1)
        
        IF (prvAsy.NE.0 .AND. prvCnt.NE.0) THEN
          DcmpAsyLinkInfo(1, FORWARD, prvCnt, prvAsy) = iAsy
          DcmpAsyLinkInfo(2, FORWARD, prvCnt, prvAsy) = iCnt
          
          DcmpAsyLinkInfo(1, BACKWARD, iCnt, iAsy) = prvAsy
          DcmpAsyLinkInfo(2, BACKWARD, iCnt, iAsy) = prvCnt
        END IF
        
        DcmpAsyRay(iCnt, iAsy)%iRotRay = iRotRay
        DcmpAsyRay(iCnt, iAsy)%iAsy    = iAsy
        DcmpAsyRay(iCnt, iAsy)%iRay    = iCnt
      END IF
      ! --------------------------------------------------
      DcmpAsyRay(iCnt, iAsy)%nMaxCellRay = nMaxCellRay
      DcmpAsyRay(iCnt, iAsy)%nAsyRay     = nRef
      
      DcmpAsyRay(iCnt, iAsy)%AsyRayList(nRef) = jAsyRay
      DcmpAsyRay(iCnt, iAsy)%AziList   (nRef) = hcRay(abs(jCoreRay))%AzmIdx
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
SELECT CASE (hLgc%idcmpclr)
CASE (0) ! NC
  CALL dmalloc0(DcmpAsyClr, 0, nAsy, 1, 1)
  
  DcmpAsyClr(0, 1) = nAsy
  
  DO iAsy = 1, nAsy
    DcmpAsyClr(iAsy, 1) = iAsy
    Core%Asy(iAsy)%color  = 1
  END DO
CASE (1) ! RB 1
  CALL dmalloc0(DcmpAsyClr, 0, nAsy, 1, 2)
  
  DO iAsy = 1, nAsy
    itmp = mod(iAsy, 2)
    
    SELECT CASE (itmp)
    CASE (1); iClr = RED
    CASE (0); iClr = BLACK
    END SELECT
    
    Core%Asy(iAsy)%color = iClr
    
    DcmpAsyClr(0, iClr) = DcmpAsyClr(0, iClr) + 1
    
    DcmpAsyClr(DcmpAsyClr(0, iClr), iClr) = iAsy
  END DO
CASE (2) ! RB 2
  CALL dmalloc0(DcmpAsyClr, 0, nAsy, 1, 2)
  
  DO iAsy = 1, nAsy
    itmp = mod(hAsy(iAsy)%iaX + hAsy(iAsy)%iaY, 2)
    
    SELECT CASE (itmp)
    CASE (1); iClr = RED
    CASE (0); iClr = BLACK
    END SELECT
    
    Core%Asy(iAsy)%color = iClr
    
    DcmpAsyClr(0, iClr) = DcmpAsyClr(0, iClr) + 1
    
    DcmpAsyClr(DcmpAsyClr(0, iClr), iClr) = iAsy
  END DO
CASE (3, 4) ! RGB
  CALL dmalloc0(DcmpAsyClr, 0, nAsy, 1, 3)
  
  DO iAsy = 1, nAsy
    itmp = mod(hAsy(iAsy)%iaX + hAsy(iAsy)%iaY, 3)
    
    SELECT CASE (itmp)
    CASE (1); iClr = RED
    CASE (2); iClr = BLACK
    CASE (0); iClr = GREEN
    END SELECT
    
    Core%Asy(iAsy)%color = iClr
    
    DcmpAsyClr(0, iClr) = DcmpAsyClr(0, iClr) + 1
    
    DcmpAsyClr(DcmpAsyClr(0, iClr), iClr) = iAsy
  END DO
END SELECT
! ----------------------------------------------------
! DEBUG
!nRef = 0
!
!DO iAsy = 1, nAsy
!  DO iCnt = 1, nModRay
!    nRef = nRef + DcmpAsyRay(iCnt, iAsy)%nAsyRay
!  END DO
!END DO
!
!maxNumPin = 0
!maxNumFsr = 0
!
!DO iAsy = 1, nAsy
!  PinSt = hAsy(iAsy)%PinIdxSt
!  PinEd = hAsy(iAsy)%PinIdxSt + hAsy(iAsy)%nTotPin - 1
!  
!  maxNumPin = max(maxNumPin, PinEd - PinSt + 1)
!  
!  DO iz = 1, Core%nz
!    FsrSt = Core%Pin(PinSt)%FsrIdxSt
!    FsrEd = Core%Pin(PinEd)%FsrIdxSt + Core%Cellinfo(Core%Pin(PinEd)%Cell(iz))%nFsr - 1
!    
!    maxNumFsr = max(maxNumFsr, FsrEd - FsrSt + 1)
!  END DO
!END DO
!
!PRINT *, nModRay, Core%nCoreFsr, Core%nxy, nAsy, nRef, maxNumPin, maxNumFsr
!STOP
! ----------------------------------------------------
IF (hLgc%ldcmpad) THEN
  CALL dmalloc0(DcmpAziRay, 0, nModRay, 1, RayInfo%nAziAngle, 1, nAsy)
  
  DO iAsy = 1, nAsy
    DO iAzi = 1, RayInfo%nAziAngle
      DO iCnt = 1, DcmpAsyRayCount(iAsy)
        IF (DcmpAsyRay(iCnt, iAsy)%AziList(1) .NE. iAzi) CYCLE
        
        DcmpAziRay(0, iAzi, iAsy) = DcmpAziRay(0, iAzi, iAsy) + 1
        
        DcmpAziRay(DcmpAziRay(0, iAzi, iAsy), iAzi, iAsy) = iCnt
      END DO
    END DO
  END DO
END IF

RayInfo%DcmpAsyRay => DcmpAsyRay
! ----------------------------------------------------
NULLIFY (DcmpAsyRayCount)
NULLIFY (AsyRayList)
NULLIFY (DirList)
NULLIFY (AziList)
NULLIFY (DcmpAsyLinkInfo)
NULLIFY (Pin)
NULLIFY (haRay_Loc)
NULLIFY (CelRay_Loc)
NULLIFY (hRotRay_Loc)
! ----------------------------------------------------

END SUBROUTINE HexDcmpRayGen