SUBROUTINE DcplLogOut(iRefPln, itemp)
USE PARAM
USE TYPEDEF
USE GEOM,            ONLY : Core, ng
USE RAYS,            ONLY : RayInfo
USE Core_mod,        ONLY : GroupInfo        ,CmInfo            ,FmInfo                 &
                           ,THInfo           ,Eigv              ,GcGroupInfo
USE DcplCore_Mod,    ONLY : DcplInfo         ,DcplFmInfo        ,DcplCmInfo             &
                           ,DcplThInfo       ,DcplFxr
USE PE_MOD,          ONLY : PE               ,DcplPE
USE CNTL,            ONLY : nTracerCntl      ,DcplControl
USE SubGrp_mod,      ONLY : SubGrpFsp
USE DcplXsGen_Mod,   ONLY : RefPlnXsGeneration ,DcplConvChk     ,ThXsGeneration
USE DcplTH_Mod,      ONLY : DcplThUpdate
USE ItrCNTL_mod,     ONLY : ItrCntl, DcplItrCntl
USE FILES,           ONLY : io8
USE IOUTIL,          ONLY : message
IMPLICIT NONE
INTEGER :: iRefPln

INTEGER :: ipln, ixy, iz, nxy, ifxr,ifsr, niso
INTEGER:: FxrIdxSt, FsrIdxSt, icel, nlocalFxr,io
INTEGER :: i, j, k, ig, itemp
INTEGER :: idum(100)
REAL :: rdum(0:100)

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :), myFxr
TYPE(Pin_Type), POINTER :: Pin(:)

Fxr => DcplInfo%DcplFmInfo(itemp, iRefPln)%Fxr
ipln = DcplPE(iRefPln)%myzb
iz =ipln
nxy = Core%nxy
Pin => Core%Pin
io = 100+irefpln
OPEN(UNIT =io, FILE = 'Fxr.out', FORM='UNFORMATTED', status ='replace', position='append')
DO ixy = 1, nxy
  FsrIdxSt = Pin(ixy)%FsrIdxSt; FxrIdxSt = Pin(ixy)%FxrIdxSt
  icel = Pin(ixy)%Cell(iz); nlocalFxr = Core%CellInfo(icel)%nFxr
  DO j =1, nLocalFxr
    ifxr = FxrIdxSt + j -1
    myFxr => Fxr(ifxr, iz)
    niso=myFxr%niso
    WRITE(io) niso
    WRITE(io) (myFxr%idiso(k), k = 1, niso)
    WRITE(io) (myFxr%pnum(k), k=1,niso)
    WRITE(io) myFxr%temp
  ENDDO
ENDDO
WRITE(io) DcplInfo%eigv(itemp, iRefPln)
DO ixy = 1, nxy
  FsrIdxSt = Pin(ixy)%FsrIdxSt; FxrIdxSt = Pin(ixy)%FxrIdxSt
  icel = Pin(ixy)%Cell(iz); nlocalFxr = Core%CellInfo(icel)%nFxr
  DO j =1, Core%CellInfo(icel)%nFsr
    iFsr = FsrIdxSt + j - 1
    WRITE(io) (DcplInfo%DcplFmInfo(itemp, iRefPln)%phis(iFsr, ipln, ig), ig=1, ng)
    WRITE(io) DcplInfo%DcplFmInfo(itemp, iRefPln)%psi(iFsr, ipln),DcplInfo%DcplFmInfo(1, iRefPln)%psid(iFsr, ipln)
  ENDDO
ENDDO
DO i = 1, RayInfo%nPolarAngle
  DO j = 1, RayInfo%nPolarAngle
    WRITE(io) (DcplInfo%DcplFmInfo(itemp, iRefPln)%PhiAngIn(j, i ,ipln, ig), ig = 1, ng)
  ENDDO
ENDDO

DO ixy = 1, nxy
  DO i =1,4 
    WRITE(io) (DcplInfo%DcplCmInfo(itemp, iRefPln)%PinXS(ixy, ipln)%pdhat(i,ig), ig =1, ng)
  ENDDO
ENDDO
!CALL CP_VA(PhiAngin1g, FmInfo%PhiAngin(:,: ,iz, ig), RayInfo%nPolarAngle, RayInfo%nPhiAngSv)
CLOSE(io)
PAUSE
!OPEN(UNIT =io, FILE = 'Fxr.out', FORM='UNFORMATTED', status ='OLD',action='read')
!DO ixy = 1, nxy
!  FsrIdxSt = Pin(ixy)%FsrIdxSt; FxrIdxSt = Pin(ixy)%FxrIdxSt
!  icel = Pin(ixy)%Cell(iz); nlocalFxr = Core%CellInfo(icel)%nFxr
!  DO j =1, nLocalFxr
!    ifxr = FxrIdxSt + j -1
!    myFxr => Fxr(ifxr, iz)
!    niso=myFxr%niso
!    READ(io) niso
!    READ(io) idum(2:niso+1)
!    READ(io) rdum(1:niso)
!    READ(io) rdum(0)
!    CONTINUE
!  ENDDO
!ENDDO
!CLOSE(io)
END SUBROUTINE

SUBROUTINE DcplLogOut0(iRefPln)
USE PARAM
USE TYPEDEF
USE GEOM,            ONLY : Core, ng
USE RAYS,            ONLY : RayInfo
USE Core_mod,        ONLY : GroupInfo        ,CmInfo            ,FmInfo                 &
                           ,THInfo           ,Eigv              ,GcGroupInfo
USE DcplCore_Mod,    ONLY : DcplInfo         ,DcplFmInfo        ,DcplCmInfo             &
                           ,DcplThInfo       ,DcplFxr
USE PE_MOD,          ONLY : PE               ,DcplPE
USE CNTL,            ONLY : nTracerCntl      ,DcplControl
USE SubGrp_mod,      ONLY : SubGrpFsp
USE DcplXsGen_Mod,   ONLY : RefPlnXsGeneration ,DcplConvChk     ,ThXsGeneration
USE DcplTH_Mod,      ONLY : DcplThUpdate
USE ItrCNTL_mod,     ONLY : ItrCntl, DcplItrCntl
USE FILES,           ONLY : io8
USE IOUTIL,          ONLY : message
IMPLICIT NONE
INTEGER :: iRefPln

INTEGER :: ipln, ixy, iz, nxy, ifxr,ifsr, niso
INTEGER:: FxrIdxSt, FsrIdxSt, icel, nlocalFxr,io
INTEGER :: i, j, k, ig
INTEGER :: idum(100)
REAL :: rdum(0:100)

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :), myFxr
TYPE(Pin_Type), POINTER :: Pin(:)

Fxr => DcplInfo%DcplFmInfo(1, iRefPln)%Fxr
ipln = DcplPE(iRefPln)%myzb
iz =ipln
nxy = Core%nxy
Pin => Core%Pin
io = 100+irefpln
OPEN(UNIT =io, FILE = 'Fxr.out', FORM='UNFORMATTED', status ='replace', position='append')
DO ixy = 1, nxy
  FsrIdxSt = Pin(ixy)%FsrIdxSt; FxrIdxSt = Pin(ixy)%FxrIdxSt
  icel = Pin(ixy)%Cell(iz); nlocalFxr = Core%CellInfo(icel)%nFxr
  DO j =1, nLocalFxr
    ifxr = FxrIdxSt + j -1
    myFxr => Fxr(ifxr, iz)
    niso=myFxr%niso
    READ(io) niso
    READ(io) (myFxr%idiso(k), k = 1, niso)
    READ(io) (myFxr%pnum(k), k=1,niso)
    READ(io) myFxr%temp
  ENDDO
ENDDO
READ(io) DcplInfo%eigv(1, iRefPln)
DO ixy = 1, nxy
  FsrIdxSt = Pin(ixy)%FsrIdxSt; FxrIdxSt = Pin(ixy)%FxrIdxSt
  icel = Pin(ixy)%Cell(iz); nlocalFxr = Core%CellInfo(icel)%nFxr
  DO j =1, Core%CellInfo(icel)%nFsr
    iFsr = FsrIdxSt + j - 1
    READ(io) (DcplInfo%DcplFmInfo(1, iRefPln)%phis(iFsr, ipln, ig), ig=1, ng)
    READ(io) DcplInfo%DcplFmInfo(1, iRefPln)%psi(iFsr, ipln),DcplInfo%DcplFmInfo(1, iRefPln)%psid(iFsr, ipln)
  ENDDO
ENDDO
DO i = 1, RayInfo%nPolarAngle
  DO j = 1, RayInfo%nPolarAngle
    READ(io) (DcplInfo%DcplFmInfo(1, iRefPln)%PhiAngIn(j, i ,ipln, ig), ig = 1, ng)
  ENDDO
ENDDO

DO ixy = 1, nxy
  DO i =1,4 
    READ(io) (DcplInfo%DcplCmInfo(1, iRefPln)%PinXS(ixy, ipln)%pdhat(i,ig), ig =1, ng)
  ENDDO
ENDDO
!CALL CP_VA(PhiAngin1g, FmInfo%PhiAngin(:,: ,iz, ig), RayInfo%nPolarAngle, RayInfo%nPhiAngSv)
CLOSE(io)
PAUSE
!OPEN(UNIT =io, FILE = 'Fxr.out', FORM='UNFORMATTED', status ='OLD',action='read')
!DO ixy = 1, nxy
!  FsrIdxSt = Pin(ixy)%FsrIdxSt; FxrIdxSt = Pin(ixy)%FxrIdxSt
!  icel = Pin(ixy)%Cell(iz); nlocalFxr = Core%CellInfo(icel)%nFxr
!  DO j =1, nLocalFxr
!    ifxr = FxrIdxSt + j -1
!    myFxr => Fxr(ifxr, iz)
!    niso=myFxr%niso
!    READ(io) niso
!    READ(io) idum(2:niso+1)
!    READ(io) rdum(1:niso)
!    READ(io) rdum(0)
!    CONTINUE
!  ENDDO
!ENDDO
!CLOSE(io)
END SUBROUTINE

SUBROUTINE DcplLogOutRead()
USE PARAM
USE TYPEDEF
USE GEOM,            ONLY : Core, ng
USE RAYS,            ONLY : RayInfo
USE Core_mod,        ONLY : GroupInfo        ,CmInfo            ,FmInfo                 &
                           ,THInfo           ,Eigv              ,GcGroupInfo
USE DcplCore_Mod,    ONLY : DcplInfo         ,DcplFmInfo        ,DcplCmInfo             &
                           ,DcplThInfo       ,DcplFxr
USE PE_MOD,          ONLY : PE               ,DcplPE
USE CNTL,            ONLY : nTracerCntl      ,DcplControl
USE SubGrp_mod,      ONLY : SubGrpFsp
USE DcplXsGen_Mod,   ONLY : RefPlnXsGeneration ,DcplConvChk     ,ThXsGeneration
USE DcplTH_Mod,      ONLY : DcplThUpdate
USE ItrCNTL_mod,     ONLY : ItrCntl, DcplItrCntl
USE FILES,           ONLY : io8
USE IOUTIL,          ONLY : message
IMPLICIT NONE
INTEGER :: iRefPln

INTEGER :: ipln, ixy, iz, nxy, ifxr,ifsr, niso
INTEGER:: FxrIdxSt, FsrIdxSt, icel, nlocalFxr,io
INTEGER :: i, j, k, ig
INTEGER :: idum(100)
REAL :: rdum(0:100), minv, maxv

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :), myFxr
TYPE(Pin_Type), POINTER :: Pin(:)

Fxr =>FmInfo%Fxr
iz =1; ipln = 1
nxy = Core%nxy
Pin => Core%Pin
io = 100
OPEN(UNIT =io, FILE = 'Fxr.out', FORM='UNFORMATTED', status ='OLD',action='read')
DO ixy = 1, nxy
  FsrIdxSt = Pin(ixy)%FsrIdxSt; FxrIdxSt = Pin(ixy)%FxrIdxSt
  icel = Pin(ixy)%Cell(iz); nlocalFxr = Core%CellInfo(icel)%nFxr
  DO j =1, nLocalFxr
    ifxr = FxrIdxSt + j -1
    myFxr => Fxr(ifxr, iz)
    niso=myFxr%niso
    READ(io) niso
    READ(io) (myFxr%idiso(k), k = 1, niso)
    READ(io) (myFxr%pnum(k), k=1,niso)
    READ(io) myFxr%temp
  ENDDO
ENDDO
READ(io) eigv
minv = 000000000000
maxv = -1000000000
DO ixy = 1, nxy
  FsrIdxSt = Pin(ixy)%FsrIdxSt; FxrIdxSt = Pin(ixy)%FxrIdxSt
  icel = Pin(ixy)%Cell(iz); nlocalFxr = Core%CellInfo(icel)%nFxr
  DO j =1, Core%CellInfo(icel)%nFsr
    iFsr = FsrIdxSt + j - 1
    READ(io) (FmInfo%phis(iFsr, ipln, ig), ig=1, ng)
    READ(io) FmInfo%psi(iFsr, ipln), FmInfo%psid(iFsr, ipln)
    !FmInfo%phis(iFsr, ipln, :)=1;
   ! max
  ENDDO
ENDDO
DO i = 1, RayInfo%nPolarAngle
  DO j = 1, RayInfo%nPolarAngle
    READ(io) (FmInfo%PhiAngIn(j, i ,ipln, ig), ig = 1, ng)
  ENDDO
ENDDO

DO ixy = 1, nxy
  DO i =1,4 
    READ(io) (CmInfo%PinXS(ixy, ipln)%pdhat(i,ig), ig =1, ng)
  ENDDO
ENDDO
!CALL CP_VA(PhiAngin1g, FmInfo%PhiAngin(:,: ,iz, ig), RayInfo%nPolarAngle, RayInfo%nPhiAngSv)
CLOSE(io)
!PAUSE

END SUBROUTINE