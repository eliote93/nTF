#include <defines.h>
!--- CNJ Edit : Domain Decomposition
SUBROUTINE DcmpRayGen(Core, RayInfo, DcmpAsyRay)
USE PARAM
USE TYPEDEF,        ONLY : RayInfo_type,      RotRayInfo_type,      CoreRayInfo_type,                       &
                           AsyRayInfo_Type,   CoreInfo_type,        DcmpAsyRayInfo_Type,                    &
                           Asy_Type,          Pin_Type,             Cell_Type
USE ALLOCS
USE PE_Mod,         ONLY : PE
USE MOC_MOD,        ONLY : nMaxDcmpRaySeg,    nMaxDcmpCellRay,      nMaxDcmpAsyRay
IMPLICIT NONE
TYPE(CoreInfo_type) :: Core
TYPE(RayInfo_type) :: RayInfo
TYPE(DcmpAsyRayInfo_type), POINTER :: DcmpAsyRay(:, :)

TYPE(RotRayInfo_Type), POINTER :: RotRay(:)
TYPE(CoreRayInfo_Type), POINTER :: CoreRay(:)
TYPE(AsyRayInfo_Type), POINTER :: AsyRay(:)
TYPE(Asy_Type), POINTER :: Asy(:)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)

INTEGER, POINTER :: DcmpAsyRayIdx(:), DcmpAsyAziList(:, :, :), DcmpAsyLinkInfo(:, :, :, :)
INTEGER, POINTER :: AsyRayList(:), DirList(:), AziList(:)
INTEGER :: i, j, k, l
INTEGER :: AsyRayBeg, AsyRayEnd, AsyRayInc
INTEGER :: nRotRay, nCoreRay, nAsyRay, nModRay, nDummyRay, nAsy
INTEGER :: nMaxAziModRay, nMaxCellRay, nMaxRaySeg, nPinRay, nRaySeg, nRaySeg0
INTEGER :: nAziAngle, nPolarAngle
INTEGER :: iRotRay, iCoreRay, iAsyRay, iRay, iceray, iAzi, iDir, iz
INTEGER :: iasy, icel, ibcel, ipin
INTEGER :: myzb, myze
INTEGER :: prevAsy, prevRay
INTEGER :: Reflection

nAsy = Core%nxya
Asy => Core%Asy
Pin => Core%Pin
Cell => Core%CellInfo

nModRay = RayInfo%nModRay
nRotRay = RayInfo%nRotRay
RotRay => RayInfo%RotRay
CoreRay => RayInfo%CoreRay
AsyRay => RayInfo%AsyRay

nAziAngle = RayInfo%nAziAngle
nPolarAngle = RayInfo%nPolarAngle
myzb = PE%myzb; myze = PE%myze

ALLOCATE(DcmpAsyRay(nModRay, nAsy))
CALL Dmalloc(RayInfo%DcmpAsyRayCount, nAsy)
CALL Dmalloc(RayInfo%DcmpAsyLinkInfo, 2, 2, nModRay, nAsy)
CALL Dmalloc0(RayInfo%DcmpAsyAziList, 0, nModRay, 1, nAziAngle / 2, 1, nAsy)

DcmpAsyRayIdx => RayInfo%DcmpAsyRayCount
DcmpAsyAziList => RayInfo%DcmpAsyAziList
DcmpAsyLinkInfo => RayInfo%DcmpAsyLinkInfo

nMaxDcmpRaySeg = 0; nMaxDcmpCellRay = 0; nMaxDcmpAsyRay = 0

DO i = 1, nRotRay
  iRotRay = i
  nCoreRay = RotRay(i)%nRay
  iRay = 0; nAsyRay = 0
  prevAsy = 0; prevRay = 0
  DO j = 1, nCoreRay
    iCoreRay = RotRay(iRotRay)%RayIdx(j)
    nAsyRay = CoreRay(iCoreRay)%nRay
    nDummyRay = 0
    iDir = RotRay(i)%Dir(j)
    IF (iDir .EQ. backward) THEN
      AsyRayBeg = nAsyRay; AsyRayEnd = 1; AsyRayInc = -1
    ELSE
      AsyRayBeg = 1; AsyRayEnd = nAsyRay; AsyRayInc = 1
    ENDIF
    DO k = AsyRayBeg, AsyRayEnd, AsyRayInc
      iRay = iRay + 1
      iAsyRay = CoreRay(iCoreRay)%AsyRayIdx(k)
      iAsy = CoreRay(iCoreRay)%AsyIdx(k)
      IF (iAsy .EQ. 0) THEN
        nDummyRay = nDummyRay + AsyRayInc; CYCLE
      ENDIF
      IF (iAsy .EQ. prevAsy) THEN
        Reflection = Reflection + 1
        nRaySeg = 0
        DO l = 1, AsyRay(iAsyRay)%nCellRay
          ipin = AsyRay(iAsyRay)%PinIdx(l)
          iceray = AsyRay(iAsyRay)%PinRayIdx(l)
          ipin = Asy(iAsy)%GlobalPinIdx(ipin)
          nRaySeg0 = 0
          DO iz = myzb, myze
            icel = Pin(ipin)%cell(iz)
            ibcel = Cell(iCel)%basecellstr
            nRaySeg0 = max(nRaySeg0, Cell(ibcel)%Cellray(iceray)%nSeg)
          ENDDO
          nRaySeg = nRaySeg + nRaySeg0
        ENDDO
        nMaxRaySeg = max(nRaySeg, nMaxRaySeg)
        nMaxDcmpRaySeg = max(nMaxRaySeg, nMaxDcmpRaySeg)
        nMaxCellRay = max(AsyRay(iAsyRay)%nCellRay, nMaxCellRay)
        nMaxDcmpCellRay = max(nMaxCellRay, nMaxDcmpCellRay)
        nMaxDcmpAsyRay = max(Reflection, nMaxDcmpAsyRay)
        DcmpAsyRay(DcmpAsyRayIdx(iAsy), iAsy)%nMaxRaySeg = nMaxRaySeg
        DcmpAsyRay(DcmpAsyRayIdx(iAsy), iAsy)%nMaxCellRay = nMaxCellRay
        DcmpAsyRay(DcmpAsyRayIdx(iAsy), iAsy)%nAsyRay = Reflection
        AsyRayList => DcmpAsyRay(DcmpAsyRayIdx(iAsy), iAsy)%AsyRayList
        DirList => DcmpAsyRay(DcmpAsyRayIdx(iAsy), iAsy)%DirList
        AziList => DcmpAsyRay(DcmpAsyRayIdx(iAsy), iAsy)%AziList
        DcmpAsyRay(DcmpAsyRayIdx(iAsy), iAsy)%AsyRayList => NULL()
        DcmpAsyRay(DcmpAsyRayIdx(iAsy), iAsy)%DirList => NULL()
        DcmpAsyRay(DcmpAsyRayIdx(iAsy), iAsy)%AziList => NULL()
        ALLOCATE(DcmpAsyRay(DcmpAsyRayIdx(iAsy), iAsy)%AsyRayList(Reflection))
        ALLOCATE(DcmpAsyRay(DcmpAsyRayIdx(iAsy), iAsy)%DirList(Reflection))
        ALLOCATE(DcmpAsyRay(DcmpAsyRayIdx(iAsy), iAsy)%AziList(Reflection))
        DcmpAsyRay(DcmpAsyRayIdx(iAsy), iAsy)%AsyRayList(1 : Reflection - 1) = AsyRayList(:)
        DcmpAsyRay(DcmpAsyRayIdx(iAsy), iAsy)%DirList(1 : Reflection - 1) = DirList(:)
        DcmpAsyRay(DcmpAsyRayIdx(iAsy), iAsy)%AziList(1 : Reflection - 1) = AziList(:)
        DEALLOCATE(AsyRayList, DirList, AziList)
        NULLIFY(AsyRayList, DirList, AziList)
      ELSE
        Reflection = 1
        DcmpAsyRayIdx(iAsy) = DcmpAsyRayIdx(iAsy) + 1
        DcmpAsyRay(DcmpAsyRayIdx(iAsy), iAsy)%nAsyRay = Reflection
        nRaySeg = 0
        DO l = 1, AsyRay(iAsyRay)%nCellRay
          ipin = AsyRay(iAsyRay)%PinIdx(l)
          iceray = AsyRay(iAsyRay)%PinRayIdx(l)
          ipin = Asy(iAsy)%GlobalPinIdx(ipin)
          nRaySeg0 = 0
          DO iz = myzb, myze
            icel = Pin(ipin)%cell(iz)
            ibcel = Cell(iCel)%basecellstr
            nRaySeg0 = max(nRaySeg0, Cell(ibcel)%Cellray(iceray)%nSeg)
          ENDDO
          nRaySeg = nRaySeg + nRaySeg0
        ENDDO
        nMaxRaySeg = nRaySeg
        nMaxDcmpRaySeg = max(nMaxRaySeg, nMaxDcmpRaySeg)
        nMaxCellRay = AsyRay(iAsyRay)%nCellRay
        nMaxDcmpCellRay = max(nMaxCellRay, nMaxDcmpCellRay)
        nMaxDcmpAsyRay = max(Reflection, nMaxDcmpAsyRay)
        DcmpAsyRay(DcmpAsyRayIdx(iAsy), iAsy)%nMaxRaySeg = nMaxRaySeg
        DcmpAsyRay(DcmpAsyRayIdx(iAsy), iAsy)%nMaxCellRay = nMaxCellRay
        ALLOCATE(DcmpAsyRay(DcmpAsyRayIdx(iAsy), iAsy)%AsyRayList(Reflection))
        ALLOCATE(DcmpAsyRay(DcmpAsyRayIdx(iAsy), iAsy)%DirList(Reflection))
        ALLOCATE(DcmpAsyRay(DcmpAsyRayIdx(iAsy), iAsy)%AziList(Reflection))
        IF (prevAsy .NE. 0 .AND. prevRay .NE. 0) THEN       
          DcmpAsyLinkInfo(1, forward, prevRay, prevAsy) = iAsy
          DcmpAsyLinkInfo(2, forward, prevRay, prevAsy) = DcmpAsyRayIdx(iAsy)
          DcmpAsyLinkInfo(1, backward, DcmpAsyRayIdx(iAsy), iAsy) = prevAsy
          DcmpAsyLinkInfo(2, backward, DcmpAsyRayIdx(iAsy), iAsy) = prevRay
        ENDIF
      ENDIF
      DcmpAsyRay(DcmpAsyRayIdx(iAsy), iAsy)%iRotRay = iRotRay
      DcmpAsyRay(DcmpAsyRayIdx(iAsy), iAsy)%iAsy = iAsy
      DcmpAsyRay(DcmpAsyRayIdx(iAsy), iAsy)%iRay = DcmpAsyRayIdx(iAsy)
      DcmpAsyRay(DcmpAsyRayIdx(iAsy), iAsy)%AsyRayList(Reflection) = iAsyRay
      DcmpAsyRay(DcmpAsyRayIdx(iAsy), iAsy)%AziList(Reflection) = CoreRay(iCoreRay)%iAng
      DcmpAsyRay(DcmpAsyRayIdx(iAsy), iAsy)%DirList(Reflection) = iDir
      IF (j .EQ. 1 .AND. k .EQ. (AsyRayBeg + nDummyRay)) THEN
        DcmpAsyRay(DcmpAsyRayIdx(iAsy), iAsy)%lRotRayBeg(forward) = TRUE
        DcmpAsyLinkInfo(1, backward, DcmpAsyRayIdx(iAsy), iAsy) = 0
        DcmpAsyLinkInfo(2, backward, DcmpAsyRayIdx(iAsy), iAsy) = RayInfo%PhiAngOutSvIdx(iRotRay, backward)
      ENDIF
      prevAsy = iAsy; prevRay = DcmpAsyRayIdx(iAsy)
    ENDDO
  ENDDO
  IF (prevAsy .NE. 0 .AND. prevRay .NE. 0) THEN
    DcmpAsyRay(prevRay, prevAsy)%lRotRayBeg(backward) = TRUE
    DcmpAsyLinkInfo(1, forward, prevRay, prevAsy) = 0
    DcmpAsyLinkInfo(2, forward, prevRay, prevAsy) = RayInfo%PhiAngOutSvIdx(iRotRay, forward)
  ENDIF
ENDDO

DO iAsy = 1, nAsy
    
  IF (mod(Asy(iAsy)%ixa, 2) .EQ. mod(Asy(iAsy)%iya, 2)) THEN
    Asy(iAsy)%color = red
  ELSE
    Asy(iAsy)%color = black
  ENDIF
  
  DO iAsyRay = 1, DcmpAsyRayIdx(iAsy)
    DO iAzi = 1, nAziAngle / 2
      IF (ANY(DcmpAsyRay(iAsyRay, iAsy)%AziList .EQ. iAzi .OR. DcmpAsyRay(iAsyRay, iAsy)%AziList .EQ. nAziAngle - iAzi + 1)) THEN
        DcmpAsyAziList(0, iAzi, iAsy) = DcmpAsyAziList(0, iAzi, iAsy) + 1
        DcmpAsyAziList(DcmpAsyAziList(0, iAzi, iAsy), iAzi, iAsy) = iAsyRay
        EXIT
      ENDIF
    ENDDO
  ENDDO
  
ENDDO

RayInfo%DcmpAsyRay => DcmpAsyRay

END SUBROUTINE