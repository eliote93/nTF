#include <defines.h>
MODULE CrCsp_Mod
USE PARAM
USE TYPEDEF,           ONLY : CoreInfo_Type,           FxrInfo_Type,           PE_TYPE
USE CNTL,              ONLY : nTracerCntl_Type
IMPLICIT NONE
TYPE CrCspDat_Type
  LOGICAL :: lUse = .FALSE.
  CHARACTER(256) :: fn
  INTEGER :: nzdat, ncrdat, icrbeg, icrend, ng
  INTEGER, POINTER :: icrpos(:)
  REAL, POINTER :: h(:), hz(:), crpos(:)
  REAL, POINTER :: Phi(:, :, :), Phi4th(:, :, :, :)
  
  TYPE(IntSpec_Type), POINTER :: IntSpec(:)
ENDTYPE

TYPE IntSpec_Type
  INTEGER :: nFxr
  REAL, POINTER :: PhiDist(:, :, :)
END TYPE


TYPE(CrCspDat_Type), POINTER, PRIVATE :: CrCspDat(:)
INTEGER, PRIVATE :: nCrCspDat = 0
LOGICAL, PRIVATE :: lAlloc_CrCspDat = .FALSE.
INTEGER, POINTER :: CrCspMap(:)
LOGICAL, PRIVATE :: lAlloc_CrCspMap = .FALSE.
LOGICAL, PRIVATE :: CspFlag(2) = .FALSE.

TYPE(CrCspDat_Type), POINTER, PRIVATE :: Dat
REAL, POINTER, PRIVATE :: hz0(:)
REAL, POINTER, PRIVATE :: h0(:)
INTEGER, PRIVATE :: ng , nz0

CONTAINS

SUBROUTINE MacXsSmCsp(XsMac, XsMac1, XsMac2, CspFxr, igb, ige, ng)
USE TYPEDEF,           ONLY : CspFxr_Type, XsMac_Type
USE BasicOperation,    ONLY : CP_CA
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac, XsMac1, XsMac2
TYPE(CspFxr_Type) :: CspFxr
INTEGER :: igb, ige, ng

REAL :: wt1, wt2
INTEGER :: ig, ig2

DO ig = igb, ige
  DO ig2 = 1, ng
    wt1 = CspFxr%fcsp(ig2,1); wt2 = CspFxr%fcsp(ig2, 2)
    XsMac%xsmacSm(ig2, ig) = XsMac1%xsmacSm(ig2, ig) * wt1 + XsMac2%xsmacSm(ig2, ig) * wt2
  ENDDO
  wt1 = CspFxr%fcsp(ig,1); wt2 = CspFxr%fcsp(ig, 2)
  XsMac%xsmacs(ig) = XsMac1%xsmacs(ig) * wt1 + XsMac2%xsmacs(ig) * wt2
  XsMac%xsmacstr(ig) = XsMac1%xsmacstr(ig) * wt1 + XsMac2%xsmacstr(ig) * wt2
ENDDO

END SUBROUTINE

SUBROUTINE MacXsP1SmCsp(XsMac, XsMac1, XsMac2, CspFxr, igb, ige, ng)
USE TYPEDEF,           ONLY : CspFxr_Type, XsMac_Type
USE BasicOperation,    ONLY : CP_CA
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac, XsMac1, XsMac2
TYPE(CspFxr_Type) :: CspFxr
INTEGER :: igb, ige, ng

REAL :: wt1, wt2
INTEGER :: ig, ig2

DO ig = igb, ige
  DO ig2 = 1, ng
    wt1 = CspFxr%fcsp(ig2,1); wt2 = CspFxr%fcsp(ig2, 2)
    XsMac%xsmacp1sm(ig2, ig) = XsMac1%xsmacp1sm(ig2, ig) * wt1 + XsMac2%xsmacp1sm(ig2, ig) * wt2
  ENDDO
ENDDO

END SUBROUTINE


SUBROUTINE MacXsP2SmCsp(XsMac, XsMac1, XsMac2, CspFxr, igb, ige, ng)
USE TYPEDEF,           ONLY : CspFxr_Type, XsMac_Type
USE BasicOperation,    ONLY : CP_CA
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac, XsMac1, XsMac2
TYPE(CspFxr_Type) :: CspFxr
INTEGER :: igb, ige, ng

REAL :: wt1, wt2
INTEGER :: ig, ig2

DO ig = igb, ige
  DO ig2 = 1, ng
    wt1 = CspFxr%fcsp(ig2,1); wt2 = CspFxr%fcsp(ig2, 2)
    XsMac%xsmacp2sm(ig2, ig) = XsMac1%xsmacp2sm(ig2, ig) * wt1 + XsMac2%xsmacp2sm(ig2, ig) * wt2
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE MacXsP3SmCsp(XsMac, XsMac1, XsMac2, CspFxr, igb, ige, ng)
USE TYPEDEF,           ONLY : CspFxr_Type, XsMac_Type
USE BasicOperation,    ONLY : CP_CA
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac, XsMac1, XsMac2
TYPE(CspFxr_Type) :: CspFxr
INTEGER :: igb, ige, ng

REAL :: wt1, wt2
INTEGER :: ig, ig2

DO ig = igb, ige
  DO ig2 = 1, ng
    wt1 = CspFxr%fcsp(ig2,1); wt2 = CspFxr%fcsp(ig2, 2)
    XsMac%xsmacp3sm(ig2, ig) = XsMac1%xsmacp3sm(ig2, ig) * wt1 + XsMac2%xsmacp3sm(ig2, ig) * wt2
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE MacXsStrCsp(XsMac, XsMac1, XsMac2, CspFxr, igb, ige)
USE TYPEDEF,           ONLY : CspFxr_Type, XsMac_Type
USE BasicOperation,    ONLY : CP_CA
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac, XsMac1, XsMac2
TYPE(CspFxr_Type) :: CspFxr
INTEGER :: igb, ige

REAL :: wt1, wt2
INTEGER :: ig

DO ig = igb, ige
  wt1 = CspFxr%fcsp(ig, 1); wt2 = CspFxr%fcsp(ig, 2)
  XsMac%XsMacstr(ig) = XsMac1%XsMacstr(ig) * wt1 + XsMac2%XsMacstr(ig) * wt2
ENDDO
END SUBROUTINE

SUBROUTINE MacXskfCsp(XsMac, XsMac1, XsMac2, CspFxr, igb, ige)
USE TYPEDEF,           ONLY : CspFxr_Type, XsMac_Type
USE BasicOperation,    ONLY : CP_CA
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac, XsMac1, XsMac2
TYPE(CspFxr_Type) :: CspFxr
INTEGER :: igb, ige

REAL :: wt1, wt2
INTEGER :: ig

DO ig = igb, ige
  wt1 = CspFxr%fcsp(ig, 1); wt2 = CspFxr%fcsp(ig, 2)
  XsMac%XsMackf(ig) = XsMac1%XsMackf(ig) * wt1 + XsMac2%XsMackf(ig) * wt2
ENDDO
END SUBROUTINE


SUBROUTINE MacXsNfCsp(XsMac, XsMac1, XsMac2, CspFxr, igb, ige)
USE TYPEDEF,           ONLY : CspFxr_Type, XsMac_Type
USE BasicOperation,    ONLY : CP_CA
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac, XsMac1, XsMac2
TYPE(CspFxr_Type) :: CspFxr
INTEGER :: igb, ige

REAL :: wt1, wt2
INTEGER :: ig

DO ig = igb, ige
  wt1 = CspFxr%fcsp(ig, 1); wt2 = CspFxr%fcsp(ig, 2)
  XsMac%XsMacNf(ig) = XsMac1%XsMacNf(ig) * wt1 + XsMac2%XsMacNf(ig) * wt2
ENDDO
END SUBROUTINE

SUBROUTINE MacXsBaseCsp(XsMac, XsMac1, XsMac2, CspFxr, niso0, igb, ige, lIsoXsOut)
USE TYPEDEF,           ONLY : CspFxr_Type, XsMac_Type
USE BasicOperation,    ONLY : CP_CA
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac, XsMac1, XsMac2
TYPE(CspFxr_Type) :: CspFxr

INTEGER :: igb, ige, niso0
LOGICAL :: lIsoXsOut

REAL :: wt1, wt2
INTEGER :: ig, iso, i, iloc

DO ig = igb, ige
  wt1 = CspFxr%fcsp(ig, 1); wt2 = CspFxr%fcsp(ig, 2)
  XsMac%XsMacA(ig) = XsMac1%XsMacA(ig) * wt1 + XsMac2%XsMacA(ig) *wt2
  XsMac%XsMacNf(ig) = XsMac1%XsMacNf(ig) * wt1 + XsMac2%XsMacNf(ig) *wt2
  XsMac%XsMacf(ig) = XsMac1%XsMacf(ig) * wt1 + XsMac2%XsMacf(ig) *wt2
  XsMac%XsMacKf(ig) = XsMac1%XsMacKf(ig) * wt1 + XsMac2%XsMacKf(ig) *wt2
  XsMac%XsMacS(ig) = XsMac1%XsMacS(ig) * wt1 + XsMac2%XsMacS(ig) *wt2
  XsMac%XsMacTr(ig) = XsMac1%XsMacTr(ig) * wt1 + XsMac2%XsMacTr(ig) *wt2
  XsMac%XsMacStr(ig) = XsMac1%XsMacStr(ig) * wt1 + XsMac2%XsMacStr(ig) *wt2
  XsMac%XsMacT(ig) = XsMac1%XsMacT(ig) * wt1 + XsMac2%XsMacT(ig) *wt2
ENDDO
IF(.NOT. lIsoXsOut) RETURN
CALL CP_CA(XsMac%IsoXsMacA(1:niso0, igb:ige), 0._8, niso0, ige-igb+1)
CALL CP_CA(XsMac%IsoXsMacNf(1:niso0, igb:ige), 0._8, niso0, ige-igb+1)
CALL CP_CA(XsMac%IsoXsMacKf(1:niso0, igb:ige), 0._8, niso0, ige-igb+1)
CALL CP_CA(XsMac%IsoXsMacf(1:niso0, igb:ige), 0._8, niso0, ige-igb+1)
CALL CP_CA(XsMac%IsoXsMacS0(1:niso0, igb:ige), 0._8, niso0, ige-igb+1)
CALL CP_CA(XsMac%IsoXsMactr(1:niso0, igb:ige), 0._8, niso0, ige-igb+1)

DO ig = igb, ige
  DO iso = 1, CspFxr%niso(1)
    iloc = CspFxr%isomap(iso, 1)
    XsMac%IsoXsMacA(iloc, ig) = XsMac1%IsoXsMacA(iso, ig) * CspFxr%fcsp(ig, 1)
    XsMac%IsoXsMacNf(iloc, ig) = XsMac1%IsoXsMacNf(iso, ig) * CspFxr%fcsp(ig, 1)
    XsMac%IsoXsMacKf(iloc, ig) = XsMac1%IsoXsMacKf(iso, ig) * CspFxr%fcsp(ig, 1)
    XsMac%IsoXsMacf(iloc, ig) = XsMac1%IsoXsMacf(iso, ig) * CspFxr%fcsp(ig, 1)
    XsMac%IsoXsMacS0(iloc, ig) = XsMac1%IsoXsMacS0(iso, ig) * CspFxr%fcsp(ig, 1)
    XsMac%IsoXsMactr(iloc, ig) = XsMac1%IsoXsMactr(iso, ig) * CspFxr%fcsp(ig, 1)
  ENDDO
  DO iso = 1, CspFxr%niso(2)
    iloc = CspFxr%isomap(iso, 2)
    XsMac%IsoXsMacA(iloc, ig) = XsMac%IsoXsMacA(iloc, ig) + XsMac2%IsoXsMacA(iso, ig) * CspFxr%fcsp(ig, 2)  
    XsMac%IsoXsMacNf(iloc, ig) = XsMac%IsoXsMacNf(iloc, ig) + XsMac2%IsoXsMacNf(iso, ig) * CspFxr%fcsp(ig, 2)
    XsMac%IsoXsMacKf(iloc, ig) = XsMac%IsoXsMacKf(iloc, ig) + XsMac2%IsoXsMacKf(iso, ig) * CspFxr%fcsp(ig, 2)
    XsMac%IsoXsMacf(iloc, ig) = XsMac%IsoXsMacf(iloc, ig) + XsMac2%IsoXsMacf(iso, ig) * CspFxr%fcsp(ig, 2)  
    XsMac%IsoXsMacS0(iloc, ig) = XsMac%IsoXsMacS0(iloc, ig) + XsMac2%IsoXsMacS0(iso, ig) * CspFxr%fcsp(ig, 2)
    XsMac%IsoXsMactr(iloc, ig) = XsMac%IsoXsMactr(iloc, ig) + XsMac2%IsoXsMactr(iso, ig) * CspFxr%fcsp(ig, 2)
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE SetCrCspFtn(Fxr, nFxr, CrPos, CrFrac, iasy, iz, nTRACERCntl, PE)
USE BasicOperation,        ONLY : CP_VA
IMPLICIT NONE
!TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type) :: Fxr(nFxr)
TYPE(PE_Type) :: PE
TYPE(nTracerCntl_Type) :: nTracerCntl
REAL :: CrPos, CrFrac
INTEGER :: nFxr, iasy, iz, ig

INTEGER :: idat, ifxr, ifxr0, icase1, icase2
INTEGER :: PlnIdx(3)
REAL :: htop, hbot, IntW1, IntW2, Frac(3), wt1(ng,2), wt2(ng,2)
REAL :: wt, wtbar
LOGICAL :: lMatch, lOutDat

idat = CrCspMap(iasy);  !Cusping Function Type
Dat => CrCspDat(idat) 

htop = h0(iz); hbot = h0(iz-1)

CALL SetCspFluxDat(CrPos, icase1, icase2, IntW1, lmatch, loutDat)
CALL SetCspAxPln(htop, CrPos, hbot, PlnIdx, Frac)

CALL CALFluxWeight(wt1(:, 1), wt2(:, 1), icase1, plnIdx, Frac)
CALL CALFluxWeight(wt1(:, 2), wt2(:, 2), icase2, plnIdx, Frac)
!DO ifxr = 1, ng
!  WRITE(99, '(f10.5)') wt1(ifxr,1)
!ENDDO

IntW2 = 1._8 - IntW1
Do ig = 1, ng
  wt1(ig, 1) = IntW1 * wt1(ig, 1) + IntW2 * wt1(ig, 2)  
  wt2(ig, 1) = IntW1 * wt2(ig, 1) + IntW2 * wt2(ig, 2)
ENDDO
ifxr0 = 0
DO ifxr = 1, nFxr
  IF(.NOT. Fxr(ifxr)%lCrFxr) CYCLE
  ifxr0 = ifxr0 + 1
  !wt1(1:ng,1) = 0.1; wt2(1:ng,1) = 0.9
  idat = idat + 1
  DO ig = 1, ng
    wt = wt1(ig, 1) * Dat%IntSpec(icase1)%PhiDist(ifxr0, 1, ig)
    wtbar = wt2(ig, 1) * Dat%IntSpec(icase1)%PhiDist(ifxr0, 2, ig)
    wt = wt / (wt+wtbar)
    wtbar = 1._8 -wt
    Fxr(ifxr)%CspFxr%fcsp(ig, 1) = wt     !Cr_In
    Fxr(ifxr)%cspfxr%fcsp(ig, 2) = wtbar  !Cr_out    
  ENDDO
  !CALL CP_VA(Fxr(ifxr)%CspFxr%fcsp(1:ng, 1), wt1(1:ng,1), ng)     !Cr_In
  !CALL CP_VA(Fxr(ifxr)%CspFxr%fcsp(1:ng, 2), wt2(1:ng,1), ng)     !Cr Out
ENDDO
END SUBROUTINE

SUBROUTINE CalFluxWeight(wt1, wt2, icase, PlnIdx, Frac)
USE UtilFunction,          ONLY : Lp4thIntegral
IMPLICIT NONE
REAL :: wt1(ng), wt2(ng), Frac(3)
INTEGER :: icase, PlnIdx(3)

INTEGER :: ig, iz
REAL :: CrPlnInt1, CrPlnInt2, Coeff(0:4), phi_crin, phi_crout
REAL :: x1, x2
DO ig = 1, ng
  
  !Cr Interface Integral
  x1= -1._8; x2 = -1._8 + 2._8* Frac(2);
  iz = PlnIdx(2)
  coeff = Dat%Phi4th(0:4, iz, ig, icase)
  CrPlnInt2 = Lp4thIntegral(Coeff, x2, x1) * 0.5_8 * Dat%hz(iz)   !Unrodded Regoin
  x1 = -1._8; x2 = 1
  CrPlnInt1 = Lp4thIntegral(Coeff, x2, x1) * 0.5_8 * Dat%hz(iz)   !Rodded Region
  CrPlnInt1 = CrPlnInt1 - CrPlnInt2
  Phi_crin = CrPlnInt1; Phi_Crout = CrPlnInt2
  !CR In Region Integral
  DO iz = PlnIdx(2) + 1, PlnIdx(1)
    Phi_Crin = Phi_CrIn + Dat%hz(iz) * Dat%Phi4th(0, iz, ig, icase)
  ENDDO
  iz = PlnIdx(1)
  x1 = -1._8 + 2._8* Frac(1); x2 = 1._8; coeff = Dat%Phi4th(0:4, iz, ig , icase)
  CrPlnInt1 = Lp4thIntegral(Coeff, x2, x1) * 0.5_8 * Dat%hz(iz)
  PHi_CrIn = Phi_CrIn - CrPlnInt1;
  
  DO iz = PlnIdx(3), PlnIdx(2) - 1
    Phi_CrOut = Phi_CrOut + Dat%hz(iz) * Dat%Phi4th(0, iz, ig, icase)   
  ENDDO
  iz = PlnIdx(3)
  x1= -1._8; x2 = -1._8  + 2._8* Frac(3); Coeff = Dat%Phi4th(0:4, iz, ig, icase)
  CrPlnInt2 = Lp4thIntegral(Coeff, x2, x1) * 0.5_8 * Dat%hz(iz)
  Phi_CrOut = Phi_CrOut - CrPlnInt2
  
  wt1(ig) = Phi_CrIn / (Phi_CrOut + Phi_CrIn)
  wt2(ig) = 1 - wt1(ig)
ENDDO

END SUBROUTINE

SUBROUTINE SetCspAxPln(htop, crPos, hbot, plnidx, frac)
IMPLICIT NONE
REAL :: hbot, htop, crpos
INTEGER :: plnIdx(3)
REAL :: frac(3)

INTEGER :: iz, i
REAL :: pos(3), dh(2)

pos = (/htop, crpos, hbot/)
DO i = 1, 3
  DO iz = 1, Dat%nzdat
    dh(1) = pos(i) - Dat%h(iz-1)
    dh(2) = pos(i) - Dat%h(iz)
    IF(dh(1) .LT. 0._8 .OR. dh(2) .GE. 0._8) CYCLE
      PlnIdx(i) = iz; Frac(i) = dh(1) / Dat%hz(iz)
    EXIT
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE SetCspFluxDat(CrPos, icase1, icase2, IntW, lmatch, loutDat)
IMPLICIT NONE
REAL :: CrPos
INTEGER :: idat, icase1, icase2
LOGICAL :: lMatch           !Exist the same control rod position in the data 
LOGICAL :: loutdat          !Crposition is out of the data


INTEGER :: icrbeg, icrend, iz, icase
REAL :: dhtop, dhbot, hbot, htop
REAL :: DelZ, IntW

loutdat = .TRUE.; lmatch = .FALSE.

icrbeg = Dat%icrbeg; icrend = Dat%icrend
icase = 0
DO iz = icrbeg, icrend
  icase = icase + 1
  htop = Dat%h(iz); hbot = Dat%h(iz-1)
  dhbot = crpos - hbot; dhtop = crpos - htop
  IF(dhbot .GT. -0.0001_8 .AND. dhtop .LT. 0.0001_8) THEN
    loutdat = .FALSE.
    IF(abs(dhbot) .LT. 0.0001_8) lMatch = .TRUE.
    DelZ= abs(dhtop - dhbot)
    IntW = dhbot / DelZ
    IF(lMatch) IntW = 1.0_8
    EXIT
  ENDIF
ENDDO
IF(.NOT. loutdat) THEN
  icase1 = icase; icase2 = icase+1
  IF(icase1 .EQ. Dat%ncrdat) icase2= icase1
ENDIF
END SUBROUTINE

SUBROUTINE CalFluxSum

END SUBROUTINE

SUBROUTINE InitCrCsp(Core, GroupInfo, nTracerCntl, PE)
USE TYPEDEF,    ONLY : GroupInfo_Type
USE ioutil,     ONLY : message
USE FILES,      ONLY : io8
USE ALLOCS
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE
INTEGER :: iz

IF(.NOT. CspFlag(1) .or. .NOT. CspFlag(2)) THEN
  WRITE(mesg, '(A)') '  Control Rod Cusping Module is not active... '
  IF(PE%MASTER) CALL message(io8, .TRUE., .TRUE. , MESG)
  RETURN
ENDIF
nTracerCntl%lCrCsp = .TRUE.
nTracerCntl%lCrCspFtn = .TRUE.
IF(PE%MASTER) THEN
  WRITE(mesg, '(A)') 'Initialize Control Rod Cusping Module...'
  IF(PE%MASTER) CALL message(io8, .TRUE., .TRUE. , MESG)
ENDIF

ng = GroupInfo%ng; hz0 => Core%hz
nz0 = Core%nz; 
CALL DMALLOC0(h0, 0, nz0)
h0(0) = 0
DO iz = 1, nz0
  h0(iz) = h0(iz-1) + hz0(iz)
ENDDO


CALL ReadCrCspDat(PE)
END SUBROUTINE

SUBROUTINE GetInp_CSPFILE(oneline)
USE IOUTIL, ONLY : GetFN
IMPLICIT NONE
CHARACTER(256) :: oneline
CHARACTER(20) :: astring
INTEGER :: i

IF(.NOT. lAlloc_CrCspDat) THEN
  ALLOCATE(CrCspDat(nCrCspDat))
  lAlloc_CrCspDat = .TRUE.
ENDIF

READ(oneline, *) astring, i
CALL GetFn(oneline, 3, CrCspDat(i)%Fn)
CrCspDat(i)%luse = .TRUE.
CspFlag(1) = .TRUE.
END SUBROUTINE  

SUBROUTINE GetInp_CSPMAP(nxya, nxa, nya, Indev, outdev, PE)
USE ioutil,         only : terminate,   toupper,       IFnumeric,   nfields,     message
USE inputcards ,    only : oneline,     probe
USE BasicOperation, ONLY : CP_CA,       CP_VA
USE ALLOCS
IMPLICIT NONE
INTEGER, INTENT(IN) :: nxya, nxa, nya, InDev, OutDev
TYPE(PE_TYPE), INTENT(IN) :: PE

INTEGER :: inpdat(nxa*nya)
INTEGER :: ixa, iya, ixya, ndat, nFieldsLine
INTEGER :: i 
IF(.NOT. lAlloc_CrCspMap) THEN
  CALL DMALLOC(CrCspMap, nxya)
  lAlloc_CrCspMap = .TRUE.
ENDIF
ndat = 0; ixa = 0; iya = 0; ixya=0
DO WHILE(.TRUE.)
  read(indev, '(a256)') oneline
  IF(PE%Master) CALL message(outdev,FALSE, FALSE,oneline)
  IF(probe.eq.BANG) cycle;     IF(probe.eq.POUND) cycle  
  iya = iya + 1
  nFieldsLine = nfields(oneline)
  READ(oneline, *) (inpdat(i), i = 1, nFieldsLine)
  
  DO ixa = 1, nFieldsLine
    ixya = ixya + 1
    CrCspMap(ixya)=inpdat(ixa)
    IF(inpdat(ixa) .NE. 0) ndat = ndat + 1
  ENDDO
  IF(iya .EQ. nya) EXIT
ENDDO
CspFlag(2) = .TRUE.
END SUBROUTINE

SUBROUTINE SetnCrCspDat(i)
IMPLICIT NONE
INTEGER, INTENT(IN) :: i
nCrCspDat = MAX(i, nCrCspDat)
END SUBROUTINE

SUBROUTINE ReadCrCspDat(PE)
USE FILES,    ONLY : io15,                localfn
USE IOUTIL,   ONLY : OPENFILE,            TERMINATE
USE ALLOCS
#ifdef MPI_ENV
USE MPIComm_Mod, ONLY : MPIWaitTurn
#endif
IMPLICIT NONE

TYPE(PE_TYPE) :: PE

INTEGER :: i, j, k,  icase, iz, ig
INTEGER :: nfxr
REAL :: crpos
#ifdef MPI_ENV
CALL MPIWaitTurn(PE%MPI_Cmfd_Comm, PE%myCmfdRank, PE%nCmfdproc, .FALSE.)
#endif
DO i = 1, nCrCspDat
  IF(.NOT. CrCspDat(i)%luse)  CYCLE
  localfn = CrCspDat(i)%FN
  CALL OPENFILE(io15, .TRUE., .TRUE., .FALSE., localfn)
  READ(io15) CrCspDat(i)%ncrdat, CrCspDat(i)%icrbeg, CrCspDat(i)%icrend, CrCspDat(i)%nzdat, CrCspDat(i)%ng
  CALL DMALLOC0(CrCspDat(i)%h, 0, CrCspDat(i)%nzdat)
  CALL DMALLOC(CrCspDat(i)%hz, CrCspDat(i)%nzdat)
  CALL DMALLOC(CrCspDat(i)%crpos, CrCspDat(i)%ncrdat)
  CALL DMALLOC(CrCspDat(i)%icrpos, CrCspDat(i)%ncrdat)
  CALL DMALLOC(CrCspDat(i)%phi, CrCspDat(i)%nzdat, CrCspDat(i)%ng,CrCspDat(i)%ncrdat)
  CALL DMALLOC0(CrCspDat(i)%phi4th, 0, 4, 1, CrCspDat(i)%nzdat, 1, CrCspDat(i)%ng, 1, CrCspDat(i)%ncrdat)
  
  READ(io15)(CrCspDat(i)%h(j), j = 0, CrCspDat(i)%nzdat)
  DO k = 1, CrCspDat(i)%ncrdat
    READ(io15) icase, j, crpos
    CrCspDat(i)%crpos(icase)=crpos; CrCspDat(i)%iCrPos(icase) = j
    DO ig = 1, CrCspDat(i)%ng
      READ(io15) (CrCspDat(i)%Phi(iz, ig, icase), iz = 1, CrCspDat(i)%nzdat)
    ENDDO
  ENDDO
  DO k = 1, CrCspDat(i)%ncrdat
    READ(io15) icase, j, crpos
    DO ig = 1, CrCspDat(i)%ng
      READ(io15) (CrCspDat(i)%Phi4th(0:4, iz, ig, icase), iz = 1, CrCspDat(i)%nzdat)
    ENDDO
  ENDDO
  
  DO iz = 1, CrCspDat(i)%nzdat
    CrCspDat(i)%hz(iz) = CrCspDat(i)%h(iz) - CrCspDat(i)%h(iz-1)
  ENDDO
  
  ALLOCATE(CrCspDat(i)%IntSpec(CrCspDat(i)%ncrdat))
  DO iz = 1, CrCspDat(i)%ncrdat
    READ(io15) nfxr
    CrCspDat(i)%IntSpec(iz)%nFxr = nFxr
    CALL DMALLOC(CrCspDat(i)%IntSpec(iz)%PhiDist, nfxr, 2, ng)
    DO ig = 1, ng
      READ(io15) (CrCspDat(i)%IntSpec(iz)%PhiDist(k, 1, ig), k=1, nFxr) !In
      READ(io15) (CrCspDat(i)%IntSpec(iz)%PhiDist(k, 2, ig), k=1, nFxr) !Out
    ENDDO
  ENDDO
  CLOSE(io15)
ENDDO
#ifdef MPI_ENV
CALL MPIWaitTurn(PE%MPI_Cmfd_Comm, PE%myCmfdRank, PE%nCmfdproc, .TRUE.)
#endif
END SUBROUTINE

SUBROUTINE CrCspH2OUpdate(myFxr, nh, nh2o)
IMPLICIT NONE
TYPE(FxrInfo_Type) :: myFxr
REAL :: nh, nh2o
REAL :: nh_modi, nh2o_modi
INTEGER :: k
INTEGER :: niso
DO k = 1, 2
  IF(myFxr%CspFxr%lh2o(k)) THEN
    nh_modi = myFxr%CspFxr%h2ofrac(k) * nh; nh2o_modi = myFxr%CspFxr%h2ofrac(k) * nh2o
    niso = myFxr%CspFxr%niso(k)
    CALL H2ODensityUpdate(nH2O_modi, nH_modi, niso, myFxr%CspFxr%IsoList(1:niso,k), myFxr%CspFxr%pnum(1:niso, k))
  ENDIF  
ENDDO
END SUBROUTINE

SUBROUTINE CspFxrBoronUpdate(myFxr, BoronPpm)
USE PARAM
USE TypeDef,       ONLY : Fxrinfo_type,   CspFxr_Type
USE NuclidMap_mod, ONLY : AtomicWeight
USE Boron_mod,     ONLY : b10frac
IMPLICIT NONE

TYPE(FxrInfo_Type) :: myFxr
REAL :: Boronppm

INTEGER :: id, i, j, k, niso
REAL :: rho, aw
REAL :: tp1, tp2
LOGICAL :: lBoronIstope

REAL, POINTER :: pnum(:)
INTEGER, POINTER :: idiso(:)
!
LOGICAL :: lBoronAdd(2)
real :: wt
DO j = 1, 2
  lBoronAdd(j) = .FALSE.
  IF(.NOT. myFxr%cspFxr%lh2o(j)) CYCLE
  idiso => myFxr%CspFxr%isolist(:, j)
  pnum => myFxr%CspFxr%pnum(:, j)
  niso = myFxr%CspFxr%niso(j)
  rho = 0
  DO i = 1, niso
    id = idiso(i)/1000
    IF(id .EQ. 1 .OR. id .EQ. 8) THEN
      aw = AtomicWeight(idiso(i))
      rho = rho + pnum(i) * aw / AVOGADRO
    ENDIF
  ENDDO
  lBoronIstope = FALSE
  DO i = 1, niso
    !id = myFxr%idiso(i)/1000
    id = idiso(i)
    SELECT CASE(id)
    CASE(5010)
      aw = AtomicWeight(5010)
      pnum(i) = rho * boronppm * epsm6 * avogadro / aw * b10frac
      lBoronIstope = TRUE    
    CASE(5011)
      aw = AtomicWeight(5011)
      pnum(i) = rho * boronppm * epsm6 * avogadro / aw * (1._8-b10frac)
      lBoronIstope = TRUE
    END SELECT
!    IF(id .EQ. 5000) THEN
!      pnum(i) = rho * boronppm * epsm6 * avogadro / awboron
!      lBoronIstope = TRUE
!#define natb
!#ifndef natb        
!    ELSEIF(id .EQ. 5010) THEN
!      pnum(i) = rho * boronppm * epsm6 * avogadro / awboron * b10frac
!      lBoronIstope = TRUE      
!    ELSEIF(id .EQ. 5011) THEN
!      pnum(i) = rho * boronppm * epsm6 * avogadro / awboron * (1._8-b10frac)
!       lBoronIstope = TRUE      
!#endif
!    ENDIF
  ENDDO
  lBoronAdd(j) = .NOT. lBoronIstope
  IF(.NOT. lBoronIstope) THEN
!#ifdef natb
!    niso = niso + 1; i = niso
!    idiso(i) = 5000
!    pnum(i) = rho * boronppm * epsm6 * avogadro / awboron
!    DO k = 1, myFxr%niso
!      IF(myFxr%idiso(k) .NE. 5000) CYCLE
!      myFxr%CspFxr%isomap(i, j) = k
!      EXIT
!    ENDDO
!#else
    niso = niso + 2; i = niso-1
    idiso(i) = 5010;  aw = AtomicWeight(idiso(i))
    pnum(i) = rho * boronppm * epsm6 * avogadro / awboron * b10frac
    DO k = 1, myFxr%niso
      IF(myFxr%idiso(k) .NE. 5010) CYCLE
      myFxr%CspFxr%isomap(i, j) = k
      EXIT
    ENDDO      
    i = niso
    idiso(i) = 5011;  aw = AtomicWeight(idiso(i))
    pnum(i) = rho * boronppm * epsm6 * avogadro / awboron * (1._8 - b10frac)
    DO k = 1, myFxr%niso
      IF(myFxr%idiso(k) .NE. 5011) CYCLE
      myFxr%CspFxr%isomap(i, j) = k
      EXIT
    ENDDO
!#endif
    myFxr%CspFxr%niso(j) = niso
  ENDIF
ENDDO
!Update Volume Weighted Number Density

myFxr%idiso(1:myFxr%niso) = 0; myFxr%pnum(1:myFxr%niso) = 0
niso = 0
DO j = 1, 2
  wt = myFxr%CspFxr%vwt(j)
  DO i = 1, myFxr%CspFxr%niso(j)
    k = myFxr%CspFxr%isomap(i, j)
    niso = max(niso, k)
    myFxr%idiso(k) = myFxr%CspFxr%isolist(i, j)
    myFxr%pnum(k) = myFxr%pnum(k) + wt * myFxr%CspFxr%pnum(i, j) 
  ENDDO
ENDDO
myFxr%niso = niso

END SUBROUTINE


!SUBROUTINE UpdtCspIsotope(myFxr)
!IMPLICIT NONE
!TYPE(FxrInto_Type) :: myFxr
!
!
!END SUBROUTINE

END MODULE