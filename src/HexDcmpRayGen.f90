#include <defines.h>
SUBROUTINE HexDcmpRayGen(Core, RayInfo, DcmpAsyRay)

USE ALLOCS
USE PARAM,   ONLY : TRUE, BACKWARD, FORWARD, RED, BLACK, GREEN
USE TYPEDEF, ONLY : RayInfo_type, CoreInfo_type, DcmpAsyRayInfo_Type, Pin_Type
USE PE_Mod,  ONLY : PE
USE MOC_MOD, ONLY : nMaxDcmpRaySeg, nMaxDcmpCellRay, nMaxDcmpAsyRay
USE HexData, ONLY : hAsy
USE HexType, ONLY : Type_HexAsyRay, Type_HexCelRay, Type_HexCoreRay, Type_HexRotRay
USE HexData, ONLY : haRay, hcRay, hRotRay, hAsyTypInfo

IMPLICIT NONE

TYPE (CoreInfo_type) :: Core
TYPE (RayInfo_type)  :: RayInfo

TYPE (DcmpAsyRayInfo_type), POINTER, DIMENSION(:,:) :: DcmpAsyRay

INTEGER, POINTER, DIMENSION(:)       :: DcmpAsyRayCount, AsyRayList, DirList, AziList
INTEGER, POINTER, DIMENSION(:,:,:)   :: DcmpAsyAziList
INTEGER, POINTER, DIMENSION(:,:,:,:) :: DcmpAsyLinkInfo

INTEGER :: nRotRay, nCoreRay, nAsyRay, nModRay, nDummyRay, nAsy, nMaxAziModRay, nMaxCellRay, nMaxRaySeg, nRaySeg, nRaySeg0, nAziAngle
INTEGER :: iRotRay, iCoreRay, jCoreRay, iAsyRay, jAsyRay, icelray, iAzi, iz, iasy, ipin, iCnt, iAsyTyp, iGeoTyp, icBss, jcBss, iDir, itmp
INTEGER :: AsyRayBeg, AsyRayEnd, AsyRayInc, myzb, myze, prvAsy, prvCnt, nRef

TYPE (Pin_Type), POINTER, DIMENSION(:) :: Pin

TYPE (Type_HexAsyRay), POINTER :: haRay_Loc
TYPE (Type_HexCelRay), POINTER :: CelRay_Loc
TYPE (Type_HexRotRay), POINTER :: hRotRay_Loc
! ----------------------------------------------------

nAsy = Core%nxya
Pin => Core%Pin

nModRay   = RayInfo%nModRay
nRotRay   = RayInfo%nRotRay
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

nMaxDcmpRaySeg  = 0
nMaxDcmpCellRay = 0
nMaxDcmpAsyRay  = 0
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
      
      IF (iAsy .EQ. prvAsy) THEN
        nRef    = nRef + 1
        nRaySeg = 0
        
        DO icelray = 1, haRay_Loc%nhpRay
          iPin = haRay_Loc%CelRay(icelray)%hPinIdx
          iPin = hAsy(iAsy)%PinIdxSt + hAsyTypInfo(iAsyTyp)%PinLocIdx(iGeoTyp, iPin) - 1
          
          nRaySeg0 = 0
          DO iz = myzb, myze
            jcBss = Pin(iPin)%hCelGeo(iz)
            
            CelRay_Loc => haRay(iGeoTyp, jcBss, jAsyRay)%CelRay(icelRay)
            
            nRaySeg0 = max(nRaySeg0, CelRay_Loc%nSegRay)
          END DO
          
          nRaySeg = nRaySeg + nRaySeg0
        END DO
        
        nMaxRaySeg  = max(nMaxRaySeg,  nRaySeg) ! # of Seg. in Asy. Ray
        nMaxCellRay = max(nMaxCellRay, haRay_Loc%nhpRay)
        
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
        
        DO icelray = 1, haRay_Loc%nhpRay
          iPin = haRay_Loc%CelRay(icelray)%hPinIdx
          iPin = hAsy(iAsy)%PinIdxSt + hAsyTypInfo(iAsyTyp)%PinLocIdx(iGeoTyp, iPin) - 1
          
          nRaySeg0 = 0
          DO iz = myzb, myze
            jcBss = Pin(iPin)%hCelGeo(iz)
            
            CelRay_Loc => haRay(iGeoTyp, jcBss, jAsyRay)%CelRay(icelRay)
            
            nRaySeg0 = max(nRaySeg0, CelRay_Loc%nSegRay)
          END DO
          
          nRaySeg = nRaySeg + nRaySeg0
        END DO
        
        nMaxRaySeg  = nRaySeg
        nMaxCellRay = haRay_Loc%nhpRay
        
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
      
      nMaxDcmpRaySeg  = max(nMaxDcmpRaySeg,  nMaxRaySeg)
      nMaxDcmpCellRay = max(nMaxDcmpCellRay, nMaxCellRay)
      nMaxDcmpAsyRay  = max(nMaxDcmpAsyRay,  nRef)
      
      DcmpAsyRay(iCnt, iAsy)%nMaxRaySeg  = nMaxRaySeg
      DcmpAsyRay(iCnt, iAsy)%nMaxCellRay = nMaxCellRay
      
      DcmpAsyRay(iCnt, iAsy)%nAsyRay = nRef
      DcmpAsyRay(iCnt, iAsy)%iRotRay = iRotRay
      DcmpAsyRay(iCnt, iAsy)%iAsy    = iAsy
      DcmpAsyRay(iCnt, iAsy)%iRay    = iCnt
      
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
! Hex Color
DO iAsy = 1, nAsy
  itmp = mod(hAsy(iAsy)%iaX + hAsy(iAsy)%iaY, 3)
  
  SELECT CASE (itmp)
    CASE (1); Core%Asy(iAsy)%color = RED
    CASE (2); Core%Asy(iAsy)%color = BLACK
    CASE (0); Core%Asy(iAsy)%color = GREEN
  END SELECT
END DO
! ----------------------------------------------------
DO iAsy = 1, nAsy
  DO iCnt = 1, DcmpAsyRayCount(iAsy)
    DO iAzi = 1, nAziAngle / 2
      IF (ANY(DcmpAsyRay(iCnt, iAsy)%AziList.NE.iAzi .AND. DcmpAsyRay(iCnt, iAsy)%AziList.NE.nAziAngle - iAzi + 1)) CYCLE
      
      DcmpAsyAziList(0, iAzi, iAsy) = DcmpAsyAziList(0, iAzi, iAsy) + 1
      
      DcmpAsyAziList(DcmpAsyAziList(0, iAzi, iAsy), iAzi, iAsy) = iCnt
      
      EXIT
    END DO
  END DO
END DO

RayInfo%DcmpAsyRay => DcmpAsyRay
! ----------------------------------------------------
NULLIFY (DcmpAsyRayCount)
NULLIFY (AsyRayList)
NULLIFY (DirList)
NULLIFY (AziList)
NULLIFY (DcmpAsyAziList)
NULLIFY (DcmpAsyLinkInfo)
NULLIFY (Pin)
NULLIFY (haRay_Loc)
NULLIFY (CelRay_Loc)
NULLIFY (hRotRay_Loc)
! ----------------------------------------------------

END SUBROUTINE HexDcmpRayGen