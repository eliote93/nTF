#include <defines.h>

SUBROUTINE CspGenOut(io, Core, FmInfo, CmInfo, GroupInfo, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,      ONLY : CoreInfo_Type,      CmInfo_Type,       FmInfo_Type,         &
                         PE_TYPE,            GroupInfo_Type,    PinXS_TYPE
USE CNTL,         ONLY : nTracerCntl_Type,   OutpCntl_Type
USE CNTLROD_mod,  ONLY : CalCrFluxDip
USE FILES,        ONLY : LOCALFN,            CaseId
USE IOUTIL,       ONLY : OpenFile,           FnTrim

#ifdef MPI_ENV
USE MPIComm_Mod,  ONLY : MPIWaitTurn
#endif
IMPLICIT NONE
INTEGER :: io
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

CHARACTER(256) :: PinXsFn
TYPE(PinXS_Type), POINTER :: PinXS(:, :)
!
INTEGER :: nxy, nz, myzb, myze, ng
INTEGER :: ixy, iz, ig, igb, ige
INTEGER :: i, j
INTEGER :: bankid

nxy = Core%nxy; nz = Core%nz; ng = GroupInfo%ng
myzb = PE%myzb; myze = PE%myze
PinXS => CmInfo%PinXS 
BankId = nTracerCntl%OutpCntl%CspGenBankId


#ifdef MPI_ENV
CALL MPIWaitTurn(PE%MPI_CMFD_comm, PE%myCmfdRank, PE%nCmfdProc, .FALSE.)
#endif

DO iz = 1, nz
  IF((iz .LT. myzb) .OR. (iz .GT. myze)) CYCLE
  !localfn = trim(caseid)//'.pinxs'
  IF(iz .LT. 10) THEN
    WRITE(PinXsFn, '(A85, A10,I3, A10)') CaseID,'_pln00',iz, '.pinxs'
  ELSEIF(iz .LT. 100) THEN
    WRITE(PinXsFn, '(A85, A10,I3, A10)') CaseID,'_pln0',iz, '.pinxs'
  ELSE 
    WRITE(PinXsFn, '(A85, A10,I3, A10)') CaseID,'_pln',iz, '.pinxs'
  ENDIF
  CALL FnTrim(PinXsFn)
  CALL OpenFile(io, FALSE, TRUE, FALSE, PinXsFn)
  !Writing Header Files
  WRITE(io) nxy, ng, iz 
  DO ig = 1, ng
    WRITE(io) GroupInfo%InScatRange(1:2, ig), GroupInfo%OutScatRange(1:2, ig)
  ENDDO
  DO ixy = 1, nxy
    WRITE(io) PinXs(ixy, iz)%XSD(1:ng)
    WRITE(io) PinXs(ixy, iz)%XST(1:ng)
    WRITE(io) PinXs(ixy, iz)%XSTR(1:ng)
    WRITE(io) PinXs(ixy, iz)%XSR(1:ng)
    WRITE(io) PinXs(ixy, iz)%XSNF(1:ng)
    WRITE(io) PinXs(ixy, iz)%XSKF(1:ng)
    WRITE(io) PinXs(ixy, iz)%CHI(1:ng)
    DO j = 1, 4
      WRITE(io) PinXs(ixy, iz)%DTIL(j, 1:ng)
    ENDDO
    DO j = 1, 4
      WRITE(io) PinXs(ixy, iz)%DHAT(j, 1:ng)
    ENDDO
    DO ig = 1, ng
      igb = PinXS(ixy, iz)%xss(ig)%ib; ige = PinXs(ixy, iz)%xss(ig)%ie
      WRITE(io) igb, ige
      WRITE(io) PinXS(ixy, iz)%xss(ig)%WithInGroupScat  
      WRITE(io) PinXs(ixy, iz)%xss(ig)%from(igb:ige)
    ENDDO
  ENDDO
  CALL CalCrFluxDip(io, Core, BankId, iz, FmInfo, CmInfo, GroupInfo, nTracerCntl, PE)
  CLOSE(io)
ENDDO

#ifdef MPI_ENV
CALL MPIWaitTurn(PE%MPI_CMFD_comm, PE%myCmfdRank, PE%nCmfdProc, .TRUE.)
#endif


END SUBROUTINE
  
  
SUBROUTINE PrintPinXs(io, Core, CmInfo, GroupInfo, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,      ONLY : CoreInfo_Type,  CmInfo_Type,               &
                         PE_TYPE,        GroupInfo_Type, PinXS_TYPE
USE CNTL,         ONLY : nTracerCntl_Type,           OutpCntl_Type
USE FILES,        ONLY : LOCALFN
USE IOUTIL,       ONLY : OpenFile
#ifdef MPI_ENV
USE MPIComm_Mod,  ONLY : MPIWaitTurn
#endif
IMPLICIT NONE
!
TYPE(CoreInfo_Type) :: Core
TYPE(CmInfo_Type) :: CmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE
INTEGER :: io
!
CHARACTER(256) :: PinXsFn
TYPE(PinXS_Type), POINTER :: PinXS(:, :)
!
INTEGER :: nxy, nz, myzb, myze, ng
INTEGER :: ixy, iz, ig, igb, ige
INTEGER :: i, j
!
IF(nTracerCntl%OutpCntl%npinxsout .EQ. 0) RETURN
!REAL :: xst(ng), xsa(ng), xsnf(ng), xskf(ng), chi(ng), xss(ng, ng)

nxy = Core%nxy; nz = Core%nz; ng = GroupInfo%ng
myzb = PE%myzb; myze = PE%myze
PinXS => CmInfo%PinXS 
#ifdef MPI_ENV
CALL MPIWaitTurn(PE%MPI_CMFD_comm, PE%myCmfdRank, PE%nCmfdProc, .FALSE.)
#endif

DO i = 1, nTracerCntl%OutpCntl%npinxsout
  iz = nTracerCntl%OutpCntl%PinXsOutList(i)
  IF((iz .LT. myzb) .OR. (iz .GT. myze)) CYCLE
  PinXsFn = nTracerCntl%OutpCntl%PinXsFn(i)
  localfn = trim(PinXsFn)
  CALL OpenFile(io, FALSE, TRUE, FALSE, LocalFn)
  !Writing Header Files
  WRITE(io) nxy, ng, iz 
  !Write Group Information
!  WRITE(io) ng, GroupInfo%lUpScat, GroupInfo%UpscatRange
  DO ig = 1, ng
    WRITE(io) GroupInfo%InScatRange(1:2, ig), GroupInfo%OutScatRange(1:2, ig)
  ENDDO
  DO ixy = 1, nxy
    WRITE(io) PinXs(ixy, iz)%XSD(1:ng)
    WRITE(io) PinXs(ixy, iz)%XST(1:ng)
    WRITE(io) PinXs(ixy, iz)%XSTR(1:ng)
    WRITE(io) PinXs(ixy, iz)%XSR(1:ng)
    WRITE(io) PinXs(ixy, iz)%XSNF(1:ng)
    WRITE(io) PinXs(ixy, iz)%XSKF(1:ng)
    WRITE(io) PinXs(ixy, iz)%CHI(1:ng)
    DO j = 1, 4
      WRITE(io) PinXs(ixy, iz)%DTIL(j, 1:ng)
    ENDDO
    DO j = 1, 4
      WRITE(io) PinXs(ixy, iz)%DHAT(j, 1:ng)
    ENDDO
    DO ig = 1, ng
      igb = PinXS(ixy, iz)%xss(ig)%ib; ige = PinXs(ixy, iz)%xss(ig)%ie
      WRITE(io) igb, ige
      WRITE(io) PinXS(ixy, iz)%xss(ig)%WithInGroupScat  
      WRITE(io) PinXs(ixy, iz)%xss(ig)%from(igb:ige)
    ENDDO
  ENDDO
  
  
  CLOSE(io)
ENDDO

#ifdef MPI_ENV
CALL MPIWaitTurn(PE%MPI_CMFD_comm, PE%myCmfdRank, PE%nCmfdProc, .TRUE.)
#endif
END SUBROUTINE
  
SUBROUTINE MakeMGXS()
USE PARAM
USE geom
USE Core_mod
USE PE_Mod
USE MOC_Mod,          ONLY : FxrAvgPhi
USE BasicOperation
USE MacXsLib_Mod
IMPLICIT NONE
INTEGER :: ixy, iz, icel, ig, ig2, j, i, l, ifxr, igb, ige
INTEGER :: FxrIdxSt, nFxrLocal
REAL :: phifxr(100), localphi
REAL :: xsa(100,3), xsnf(100, 3), chi(100,3), xss(100,100, 3),xstr(100,3),xskf(100,3)
REAL :: RR(100,3), phisum(3)
TYPE(FxrInfo_Type), POINTER :: myFxr
TYPE(XsMac_Type) :: XsMac
INTEGER :: norg, iresogrpbeg, iresogrpend, nchi

ixy =1; iz=1
FxrIdxSt = Core%Pin(ixy)%FxrIdxSt; icel = Core%Pin(ixy)%Cell(iz); nFxrLocal = Core%CellInfo(icel)%nFxr
norg = GroupInfo%norg; nChi = GroupInfo%nChi
iResoGrpBeg = GroupInfo%nofg + 1
iResoGrpEnd = GroupInfo%nofg + GroupInfo%norg 
DO j= 1,3
  xsa(:,j)=0; xsnf(:,j)=0; chi(:,j)=0; xss(:, :,j)=0;xstr(:,j)=0;
ENDDO
DO ig = 1, ng
  RR = 0; phisum=0
  DO j = 1, nFxrLocal
    ifxr = FxrIdxSt + j -1
    phifxr(1:ng) = FxrAvgPhi(Core, FmInfo%Fxr, FmInfo%Phis, ixy, j, iz, ng, PE)
    myFxr => FmInfo%FXR(ifxr, iz)
    CALL MacXsBase(XsMac, myFxr, ig, ig, ng, 1._8, FALSE, TRUE, TRUE)
        
    CALL MacXsScatMatrix(XsMac, myFxr, ig, ig, ng, GroupInfo, FALSE, TRUE)
    !Self-Sheilding Effect
    IF(myFxr%lres) THEN
      IF ((ig.ge.iResoGrpBeg).and.(ig.le.iResoGrpEnd)) THEN
        XsMac%XsMacA(ig) = XsMac%XsMacA(ig) * myFxr%FresoA(ig)
        XsMac%XsMacNf(ig) = XsMac%XsMacNf(ig) * myFxr%FresoNF(ig)
        XsMac%XsMacKf(ig) = XsMac%XsMacKf(ig) * myFxr%FresokF(ig)
        ENDIF
    ENDIF
    CALL CP_CA(XsMac%CHI, 0._8, ng)
    IF(myFxr%lDepl) THEN
      CALL CP_VA(XsMac%CHI(1:nCHI), myFxr%CHI(1:nCHI), nCHI)
      XsMac%lFuel = TRUE
    ENDIF
    XsMac%XsMacTr(ig) = XsMac%XsMacA(ig) + XsMac%XsMacStr(ig)
    XsMac%XsMacT(ig) = XsMac%XsMacA(ig) + XsMac%XsMacS(ig)  
    IF(myFxr%lfuel) THEN
      l=1
    ELSEIF(myFxr%lH2o) THEN
      l=2
    ELSE
      l=3  
    ENDIF
    igb = GroupInfo%OutScatRange(1, ig)
    ige = GroupInfo%OutScatRange(2, ig)    
    localphi = phifxr(ig)*myFxr%area
    phisum(l) =phisum(l) + localphi
    RR(1,l) = RR(1,l) + localphi * XsMac%Xsmaca(ig)
    RR(2,l) = RR(2,l) + localphi * XsMac%Xsmactr(ig)
    RR(3,l) = RR(3,l) + localphi * XsMac%Xsmacnf(ig)
    RR(4,l) = RR(4,l) + localphi * XsMac%Xsmackf(ig)
    DO ig2 = igb, ige
      RR(4 + ig2, l) = RR(4 + ig2, l) + localphi * XsMac%Xsmacsm(ig, ig2)
    ENDDO    
    CHI(ig, l) = XSMAC%CHI(ig)
  ENDDO
  DO l = 1, 3
    CALL MULTI_CA(1._8/phisum(l), RR(1:ng+4,l), ng + 4)
    XSA(ig, l) = RR(1, l)
    XSTR(ig, l) = RR(2,l)
    XSNF(ig, l) = RR(3,l); XSKF(ig, l) = RR(4,l);
    XSS(1:ng, ig, l) = RR(5:4+ng, l)
  ENDDO
ENDDO
DO l = 1, 3
  DO ig = 1, ng
    WRITE(112, '(200e15.7)')  XSTR(ig, l), XSA(ig, l), XSNF(ig, l), XSKF(ig, l), chi(ig, l)
  ENDDO
  DO ig =1, ng
    WRITE(112, '(200e15.7)') XSS(1:ng, ig, l)
  ENDDO
  WRITE(112, '(200e15.7)')
ENDDO
STOP
CONTINUE
 !
END SUBROUTINE

!SUBROUTINE PrintPinXS(io, Core, CmInfo, ng, PE)
!USE PARAM
!USE TYPEDEF,   ONLY : CoreInfo_Type, CmInfo_Type, PE_TYPE,       &
!                      PinXS_TYPE
!USE FILES,  ONLY : LOCALFN, FILENAME, caseid
!USE IOUTIL, ONLY : OpenFile
!IMPLICIT NONE
!
!TYPE(CoreInfo_Type) :: Core
!TYPE(CmInfo_Type) :: CmInfo
!TYPE(PE_TYPE) :: PE
!INTEGER :: io, ng
!
!TYPE(PinXS_Type), POINTER :: PinXS(:, :)
!
!INTEGER :: nxy, myzb, myze
!INTEGER :: ixy, iz, ig, igb, ige
!
!REAL :: xst(ng), xsa(ng), xsnf(ng), xskf(ng), chi(ng), xss(ng, ng)
!
!nxy = Core%nxy
!myzb = PE%myzb; myze = PE%myze
!PinXS => CmInfo%PinXS
!
!localfn=trim(caseid)//'.pinxs'
!CALL openfile(io,FALSE,FALSE, FALSE, localfn)
!DO iz = myzb, myze
!  DO ixy = 1, nxy
!    WRITE(io, '(A7, 2I5)') 'PinXS', iz, ixy 
!    xst = PinXS(ixy, iz)%xst; xsnf = PinXS(ixy, iz)%xsnf
!    xskf = PinXS(ixy, iz)%xskf; chi = PinXS(ixy, iz)%chi
!    xss = 0
!    DO ig = 1, ng
!      igb = PinXS(ixy, iz)%xss(ig)%ib; ige = PinXS(ixy, iz)%xss(ig)%ie
!      !PinXS(ixy, iz)%xss(ig)%WithInGroupScat
!      xss(ig, igb:ige) = PinXS(ixy, iz)%xss(ig)%from(igb:ige)
!      xss(ig, ig) = PinXS(ixy, iz)%xss(ig)%WithInGroupScat
!    ENDDO
!   
!    DO ig = 1, ng
!      xsa(ig) = PinXS(ixy, iz)%xsr(ig) + xss(ig, ig)
!      xsa(ig) = xsa(ig) - SUM(xss(1:ng, ig))
!    ENDDO
!    
!    DO ig = 1, ng
!      WRITE(io, '(200e12.5)') xst(ig), xsa(ig), xsnf(ig), xskf(ig), chi(ig)
!    ENDDO
!    DO ig = 1, ng
!      WRITE(io, '(200e12.5)')  xss(1:ng, ig)
!    ENDDO
!  ENDDO
!ENDDO
!CLOSE(io)
!END SUBROUTINE
!
!SUBROUTINE WritePinXS(io, PinXS, ixy1, ixy2, iz1, iz2, ng)
!USE PARAM
!USE TYPEDEF,          ONLY : PinXS_Type
!IMPLICIT NONE
!INTEGER :: io, ixy1, ixy2, iz1, iz2, ng
!TYPE(PinXS_Type) :: PinXS(ixy1:ixy2, iz1:iz2)
!INTEGER :: iz, ig, ixy, igb, ige
!
!REAL :: xst(ng), xsa(ng), xsnf(ng), xskf(ng), chi(ng), xss(ng, ng)
!REAL :: xsd(ng), xsr(ng)
!DO iz = iz1, iz2
!  DO ixy = ixy1, ixy2
!    WRITE(io, '(A7, 2I5)') 'PinXS', iz, ixy 
!    xst = PinXS(ixy, iz)%xst; xsnf = PinXS(ixy, iz)%xsnf
!    xskf = PinXS(ixy, iz)%xskf; chi = PinXS(ixy, iz)%chi
!    xsr = PinXs(ixy, iz)%xsr
!    xsd = PinXs(ixy, iz)%XSD
!    xss = 0
!    DO ig = 1, ng
!      igb = PinXS(ixy, iz)%xss(ig)%ib; ige = PinXS(ixy, iz)%xss(ig)%ie
!      !PinXS(ixy, iz)%xss(ig)%WithInGroupScat
!      xss(ig, igb:ige) = PinXS(ixy, iz)%xss(ig)%from(igb:ige)
!      xss(ig, ig) = PinXS(ixy, iz)%xss(ig)%WithInGroupScat
!    ENDDO
!   
!    DO ig = 1, ng
!      xsa(ig) = PinXS(ixy, iz)%xsr(ig) + xss(ig, ig)
!      xsa(ig) = xsa(ig) - SUM(xss(1:ng, ig))
!    ENDDO
!    
!    DO ig = 1, ng
!      !WRITE(io, '(200e12.5)') xst(ig), xsa(ig), xsnf(ig), xskf(ig), chi(ig)
!      WRITE(io, '(200e14.5)') xsd(ig), xsr(ig), xsnf(ig), xskf(ig), chi(ig)
!    ENDDO
!    DO ig = 1, ng
!      WRITE(io, '(200e14.5)')  xss(1:ng, ig)
!    ENDDO
!  ENDDO
!ENDDO
!
!END SUBROUTINE