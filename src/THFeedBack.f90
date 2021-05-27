#include <defines.h>
SUBROUTINE THFeedBack(Core, CmInfo, FmInfo, ThInfo, nTRACERCntl, PE)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,        CMInfo_Type,      FmInfo_Type,     &
                           ThInfo_Type,          ThVar_TYpe,       ThOpt_Type,      &
                           PE_Type,                                                 &
                           FxrInfo_Type,         Pin_Type,        CoolantTH_Type,   &
                           ThVar_Type,           THOpt_Type,       FuelTh_Type,     &
                           Cell_Type,            THCell_Type,      AsyInfo_Type,    &
                           Asy_Type
                           
USE TH_Mod,         ONLY : ThVar,                ThOpt,                             &
                           CalRelPower,          SteadyCoolantTH,              &
                           Cal_RefFuelTemp,      SteadyFuelConduction
USE Boron_mod,      ONLY : SetBoronCoolant, SetBoronCoolantCTF
USE CrCsp_Mod,      ONLY : CrCspH2OUpdate
USE CNTL,           ONLY : nTracerCntl_Type
USE FILES,          ONLY : IO8
USE IOUTIL,         ONLY : message,              terminate
use timer,          ONLY : nTracer_dclock,       TimeChk
IMPLICIT NONE

INCLUDE 'mpif.h'
TYPE(CoreInfo_Type) :: Core
TYPE(CMInfo_Type) :: CmInfo
TYPE(FMInfo_Type) :: FmInfo
TYPE(ThInfo_Type) :: ThInfo
TYPE(nTracerCntl_Type) :: nTRACERCntl
TYPE(PE_Type) :: PE

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Fxrinfo_type), POINTER :: Fxr(:,:)
TYPE(FuelTh_Type), POINTER :: FuelTH(:)      
TYPE(CoolantTH_Type), POINTER :: CoolantTH(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(Asy_Type), POINTER :: Asy(:)
TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
TYPE(THCell_Type), POINTER :: ThCell
REAL, POINTER :: DenCool(:, :), TCool(:, :), Tfvol(:, :, :)
REAL, POINTER :: hz(:)

TYPE(Fxrinfo_type), POINTER :: myFxr
TYPE(Cell_Type), POINTER :: myCell

REAL :: powlin, PowLv

INTEGER :: myzb, myze, xyb, xye   !--- CNJ Edit : Domain Decomposition + MPI
INTEGER :: nzth, nxy, nbd, niso, npr1, npr2
INTEGER :: FxrIdxSt, nlocalFxr
INTEGER :: CLDReg(2), FuelReg(2), CoolReg(2)
INTEGER :: i, j, k, l, iz, ixy, icel, ifxr, ir
INTEGER :: iasy, AsyType
REAL :: rho, nh2o, nh, BulkTempK, TempK
REAL :: nh_modi, nh2o_modi
REAL :: DenIn, Tin
LOGICAL :: ChkList(nMaxFsr)

REAL :: buf0(6), buf(6)
REAL :: ncvolsum, tdsum, tdvolsum, densum, denvolsum, tmsum, tmvolsum
REAL :: tfxrsum, volfxr, volfxrsum
INTEGER :: ncell, ierr, nfsrinfxr
INTEGER :: ireg

IF(nTRACERCntl%libtyp .EQ. 11) THEN
  CALL Neacrp_thfeedback(Core, CmInfo, FmInfo, THInfo, nTracerCntl, PE) 
  RETURN
END IF
IF(nTRACERCntl%lBenchXS) RETURN
IF(.NOT. nTRACERCntl%lFeedBack) RETURN

myzb = PE%myzb; myze = PE%myze
nxy = Core%nxy
Pin => Core%Pin; Fxr => FmInfo%Fxr
Asy => Core%Asy; AsyInfo => Core%AsyInfo
CellInfo => Core%CellInfo

!--- CNJ Edit : Domain Decomposition + MPI
xyb = PE%myPinBeg; xye = PE%myPinEnd
IF (PE%RTMASTER) THEN
  xyb = 1; xye = nxy
ENDIF

DenCool => THInfo%DenCool; Tcool => ThInfo%Tcool
Tfvol => ThInfo%Tfvol
DenIn = ThInfo%DenIn; Tin = ThInfo%Tin


FuelTH => ThInfo%FuelTH; CoolantTH => ThInfo%CoolantTH
powlin = ThInfo%PowLin; PowLv = ThInfo%PowLv

nzth = ThVar%nzth; hz => ThVar%hz
npr2 = ThVar%npr2
npr1 = ThVar%npr1   !Gap
!Bottom and Top Reflector of Fuel Cell
DO iz = myzb, myze
  IF(Core%lFuelPlane(iz)) CYCLE
  DO ixy = xyb, xye
    IF(.NOT. Pin(ixy)%lfuel) CYCLE
    rho = DenCool(iz, ixy) * epsm3 ; BulkTempK = Tcool(iz, ixy) + CKELVIN
    nh2o = rho / awh2o * avogadro; nh = 2. * nh2o
    icel = Pin(ixy)%cell(iz); myCell => CellInfo(icel)
    FxrIdxSt = Pin(ixy)%FxrIdxSt; nlocalFxr = myCell%nFxr
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j - 1; myFxr => Fxr(ifxr, iz)
      myFxr%temp = BulkTempK; niso = myFxr%niso
      CALL ModFxrDensityFeedback(myFxr, nh, nh2o, BulkTempK) 
    ENDDO
  ENDDO
ENDDO

!Fuel Pin for Fuel Region
tmvolsum = 0
denvolsum = 0
tdvolsum = 0
ncvolsum = 0
tfxrsum = 0
volfxrsum = 0
DO iz = myzb, myze
  IF(.NOT. Core%lFuelPlane(iz)) CYCLE
  tmsum = 0
  densum = 0
  tdsum = 0
  ncell = 0
  DO ixy = xyb, xye
    IF(.NOT. Pin(ixy)%lfuel) CYCLE
    !Water Density Information for Moderator Region
    rho = DenCool(iz, ixy) * epsm3 ; BulkTempK = Tcool(iz, ixy) + CKELVIN
    nh2o = rho / awh2o * avogadro; nh = 2. * nh2o
    !
    icel = Pin(ixy)%cell(iz); myCell => CellInfo(icel)
    ThCell => myCell%ThCell
    FxrIdxSt = Pin(ixy)%FxrIdxSt; nlocalFxr = myCell%nFxr
    
    ! Average Information Monitoring
    ncell = ncell + 1
    tdsum = tdsum + FuelTH(ixy)%tfuel(ThVar%npr5, iz)
    densum = densum + rho
    tmsum = tmsum + Tcool(iz, ixy)

    !Initialize the checklist
    ChkList(1:nLocalFxr) = .FALSE.
    
    !Region Information
    CLDReg = myCell%ThCell%CldReg; FuelReg = myCell%ThCell%FuelReg
    CoolReg = myCell%ThCell%CoolReg
    
    !Coolant Region Temp and H2O density Update
    DO j = CoolReg(1), CoolReg(2)
      ChkList(j) = .TRUE.
      ifxr = FxrIdxSt + j - 1; myFxr => Fxr(ifxr, iz)
      myFxr%TEMP = BulkTempK
      IF(myFxr%lH2O) CALL ModFxrDensityFeedback(myFxr, nh, nh2o, BulkTempK) 
    ENDDO
    !CLAD Region Update
    DO j = CLDReg(1), CLDReg(2)
      ChkList(j) = .TRUE.
      ifxr = FxrIdxSt + j - 1; myFxr => Fxr(ifxr, iz)
      myFxr%Temp = MIN(3000._8+CKELVIN, Tfvol(npr2, iz, ixy))
    ENDDO
    !Gap Region Update (if it exists)
    DO j = CLDReg(2) + 1, FuelReg(1) -1, 1
      ChkList(j) = .TRUE.
      ifxr = FxrIdxSt + j - 1; myFxr => Fxr(ifxr, iz)
      myFxr%Temp = MIN(3000._8+CKELVIN, Tfvol(npr1, iz, ixy))      
    ENDDO
    !Fuel Region Temperature Update
    DO j = FuelReg(1), FuelReg(2)
      ChkList(j) = .TRUE.
      ifxr = FxrIdxSt + j - 1; myFxr => Fxr(ifxr, iz)
      TempK = 0
      DO ir = ThCell%FuelMapping(1, j), ThCell%FuelMapping(2, j)
        TempK = TempK + ThCell%Frac(ir, j) * tfvol(ir, iz, ixy)
      ENDDO
      myFxr%Temp = MIN(3000._8+CKELVIN, TempK)

      nfsrinfxr = myCell%nFsrInFxr(j)
      volfxr = 0
      DO k = 1, nfsrinfxr
        ireg = myCell%MapFxr2FsrIdx(k, j)
        volfxr = volfxr + myCell%vol(ireg)
      END DO
      volfxr = volfxr * hz(iz)
      tfxrsum = tfxrsum + volfxr * myFxr%Temp
      volfxrsum = volfxrsum + volfxr
    ENDDO
    !Grid Region or regions except for fuel, grid , air gap
    DO j = 1, nLocalFxr
      IF(ChkList(j)) CYCLE
      ifxr = FxrIdxSt + j - 1; myFxr => Fxr(ifxr, iz)
      myFxr%TEMP = BulkTempK
      CALL ModFxrDensityFeedback(myFxr, nh, nh2o, BulkTempK) 
    ENDDO
  ENDDO
  tdvolsum = tdvolsum + tdsum * hz(iz)
  tmvolsum = tmvolsum + tmsum * hz(iz)
  denvolsum = denvolsum + densum * hz(iz)
  ncvolsum = ncvolsum + real(ncell) * hz(iz)
ENDDO
!GUIDE Tube Information
DO iz =  myzb, myze
 ! IF(.NOT. Core%lFuelPlane(iz)) CYCLE
  DO ixy = xyb, xye
    IF(Pin(ixy)%lFuel) CYCLE
    iasy = Pin(ixy)%iasy; AsyType = Asy(iasy)%AsyType
    IF(.NOT. AsyInfo(AsyType)%lfuel) CYCLE
    !IF(.NOT. Pin(ixy)%lGT) CYCLE
    !Node Average Temperature and Density
#ifdef GTTemp_old
    nbd = Pin(ixy)%nBd; l = 0
    !Tempk = Tcool(iz, ixy); !rho = DenCool(iz, ixy)
    Tempk = 0; rho = 0
    DO j = 1, nbd
      i = Pin(ixy)%NeighIdx(j)
      IF(i .lt. 1 .AND. .NOT. Pin(i)%lfuel) CYCLE
      l = l + 1
      TempK = TempK + Tcool(iz, i);
      rho = rho + DenCool(iz, i);
    ENDDO
    IF(l .lt. 1) CALL terminate('wrong GT feedback') 
    
    rho = rho / REAL(l, 8); TempK = Tempk / REAL(l, 8)
    Tcool(iz, ixy) = TempK
#endif    
    Tempk = Tcool(iz, ixy); rho = DenCool(iz, ixy)
    !H2O Number Density Information
    rho = rho * epsm3 ; TempK = TempK + CKELVIN
    nh2o = rho / awh2o * avogadro; nh = 2. * nh2o
    !
    icel = Pin(ixy)%cell(iz); myCell => CellInfo(icel)
    ThCell => myCell%ThCell
    FxrIdxSt = Pin(ixy)%FxrIdxSt; nlocalFxr = myCell%nFxr    
    DO j = 1, nlocalFxr
      ifxr = FxrIdxSt + j - 1; myFxr => Fxr(ifxr, iz)
      myFxr%TEMP = TempK
      IF(myFxr%lh2o) CALL ModFxrDensityFeedback(myFxr, nh, nh2o, TempK) 
    ENDDO
  ENDDO
ENDDO
!Radial Reflector Information
DO iz = myzb, myze
 ! IF(.NOT. Core%lFuelPlane(iz)) CYCLE
  DO ixy = xyb, xye
    IF(Pin(ixy)%lFuel) CYCLE
    iasy = Pin(ixy)%iasy; AsyType = Asy(iasy)%AsyType
    IF(AsyInfo(AsyType)%lfuel) CYCLE
    rho = DenIn * epsm3 ; TempK = Tin + CKELVIN
    nh2o = rho / awh2o * avogadro; nh = 2. * nh2o    
    !
    icel = Pin(ixy)%cell(iz); myCell => CellInfo(icel)
    ThCell => myCell%ThCell
    FxrIdxSt = Pin(ixy)%FxrIdxSt; nlocalFxr = myCell%nFxr    
    DO j = 1, nlocalFxr
      ifxr = FxrIdxSt + j - 1; myFxr => Fxr(ifxr, iz)
      myFxr%TEMP = TempK
      IF(myFxr%lh2o) CALL ModFxrDensityFeedback(myFxr, nh, nh2o, TempK)
    ENDDO    
  ENDDO
ENDDO

IF (nTracerCntl%lborontrack) THEN
  mesg = '    The non-uniform boron concentration is considered by boron tracking model'
  IF(PE%master) CALL message(io8,FALSE, TRUE, mesg)  
	CALL SetBoronCoolantCTF(Core, Fxr,ThInfo, myzb, myze)
ELSEIF(nTracerCntl%lInitBoron) THEN
  mesg = '    The boron distrubution is uniform'
  IF(PE%master) CALL message(io8,FALSE, TRUE, mesg)
  CALL SetBoronCoolant(Core, Fxr, nTracerCntl%BoronPPM , myzb, myze)
ENDIF

#ifdef MPI_ENV
Buf0 =(/tdvolsum, denvolsum, tmvolsum, ncvolsum, tfxrsum, volfxrsum/)
CALL MPI_ALLREDUCE(Buf0, Buf, 6, MPI_DOUBLE_PRECISION, MPI_SUM, PE%MPI_CMFD_COMM, ierr)
tdvolsum = Buf(1)
denvolsum = Buf(2)
tmvolsum = Buf(3)
ncvolsum = Buf(4)
tfxrsum = Buf(5)
volfxrsum = Buf(6)
#endif
WRITE(mesg, '(2x,3(A12,es12.4))') 'Avg T.Dop:', tdvolsum/ncvolsum, 'Avg Rho:', denvolsum/ncvolsum, 'Avg T.mod:', Tmvolsum/ncvolsum
IF (PE%MASTER) CALL message(io8, .FALSE., .TRUE., mesg)
WRITE(mesg, '(2x,A16,es14.6)') 'Applied T.Dop:', tfxrsum/volfxrsum - CKELVIN
IF (PE%MASTER) CALL message(io8, .FALSE., .TRUE., mesg)

NULLIFY(Pin); NULLIFY(Fxr)
NULLIFY(FuelTH); NULLIFY(CoolantTH)
NULLIFY(CellInfo); NULLIFY(Asy)
NULLIFY(AsyInfo); NULLIFY(hz)
NULLIFY(ThCell)
NULLIFY(DenCool, TCool, Tfvol)
NULLIFY(myFxr, myCell)


END SUBROUTINE

!--- CNJ Edit : Unused
SUBROUTINE MTCH2Onumdensity(Fxr, nfxr, iz1, iz2, nTracerCntl)
USE PARAM
USE TYPEDEF, ONLY : FxrInfo_Type
USE CNTL, ONLY : nTracerCntl_Type
USE SteamTBL_mod,  ONLY : steamtbl
IMPLICIT NONE
TYPE(FxrInfo_Type) :: Fxr(nfxr, iz1:iz2)
TYPE(nTracerCntl_Type) :: nTracerCntl
INTEGER :: iz1, iz2, nFxr

REAL :: wt, wh, hin, wrho, wvin, wxin, wbetain, wkapain, wcpin
REAL :: pexit, nh2o, nh
INTEGER :: iz, ifxr
INTEGER :: niso

Pexit = nTracerCntl%Pexit
DO iz = iz1, iz2
  DO ifxr = 1, nFxr
    IF(.NOT. Fxr(ifxr, iz)%lh2o) CYCLE
    wt = nTracerCntl%TempInlet + CKELVIN + 3._8
    Fxr(ifxr, iz)%temp = wt
    CALL steamtbl(TRUE,pexit,wt,wh,wrho,wvin,wxin,wbetain,wkapain,wcpin)
    wrho = wrho* epsm3
    nh2o = wrho / awh2o * avogadro; nh = 2. * nh2o 
    niso = Fxr(ifxr, iz)%niso
    CALL H2oDensityUpdate(nh2o, nh, niso, Fxr(ifxr, iz)%idiso(1:niso), Fxr(ifxr, iz)%pnum(1:niso))
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE ModFxrDensityFeedback(myFxr, nh, nh2o, TempK)
USE PARAM
USE TYPEDEF,    ONLY : FxrInfo_Type
USE CrCsp_mod,  ONLY : CrCspH2OUpdate
IMPLICIT NONE
TYPE(FxrInfo_Type) :: myFxr
REAL :: nh, nh2o, TempK
REAL :: nh_modi, nh2o_modi

nh_modi = nh; nh2o_modi = nh2o
IF(myFxr%lMixtureMixing) THEN
  nh_modi = myFxr%h2ofrac * nh; nh2o_modi = myFxr%h2ofrac * nh2o
ENDIF        
myFxr%TEMP = TempK
CALL H2ODensityUpdate(nH2O_modi, nH_modi, myFxr%niso, myFxr%IdIso, myFxr%pnum)
!Cusping H2O Data Update
IF(myFxr%lMixtureMixing .AND. myFxr%lCrFxr .AND. myFxr%lCrCspFtn) THEN
  CALL CrCspH2OUpdate(myFxr, nh, nh2o)  !Control Rod Cusping
ENDIF  
END SUBROUTINE

SUBROUTINE H2ODensityUpdate(nh2o, nh, niso, idiso, pnum)
USE PARAM
IMPLICIT NONE
INTEGER :: niso, idiso(niso)
REAL :: nH2O, nH, pnum(niso)
INTEGER :: i, iprobe

DO i = 1, niso
  iprobe = idiso(i)/1000
  IF(iprobe .eq. 1) pnum(i) = nh
  IF(iprobe .eq. 8) pnum(i) = nh2o
ENDDO
END SUBROUTINE
  
SUBROUTINE Neacrp_thfeedback(Core, CmInfo, FmInfo, ThInfo, nTRACERCntl, PE)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,     CmInfo_Type,          FmInfo_Type,      &
                           ThInfo_Type,       ThVar_Type,           ThOpt_Type,       &
                           PE_Type,           Pin_Type,             Cell_Type,        &
                           FxrInfo_Type
USE TH_Mod,         ONLY : ThVar,             ThOpt
USE CNTL,           ONLY : nTracerCntl_Type
USE files,          ONLY : io8
USE ioutil,         ONLY : message
USE timer,          ONLY : nTracer_dclock,    TimeChk
IMPLICIT NONE

Include 'mpif.h'
TYPE(CoreInfo_Type) :: Core
TYPE(CmInfo_Type) :: CmInfo
TYPE(FmInfo_TYpe) :: FmInfo
TYPE(ThInfo_Type) :: ThInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
REAL, POINTER :: DenCool(:, :), Tcool(:, :), Tfvol(:, :, :)

REAL :: alpha = 0.7
REAL :: rho, BulkTempC, DopTempC
REAL :: hz
INTEGER :: FxrIdxSt, nlocalFxr
INTEGER :: myzb, myze, xyb, xye
INTEGER :: nxy, nsurf
INTEGER :: iz, ixy, icel, ifxr
INTEGER :: i

REAL :: buf0(4), buf(4)
REAL :: ncvolsum, tdsum, tdvolsum, densum, denvolsum, tmsum, tmvolsum
INTEGER :: ncell, ierr

myzb = PE%myzb
myze = PE%myze
nxy = Core%nxy
nsurf = ThVar%npr1

Pin => Core%Pin
Cell => Core%CellInfo
Fxr => FmInfo%Fxr
DenCool => THInfo%DenCool
Tcool => ThInfo%Tcool
Tfvol => ThInfo%Tfvol

xyb = PE%myPinBeg
xye = PE%myPinEnd
IF(PE%RTMASTER) THEN
  xyb = 1
  xye = nxy
END IF

! Bottom and Top Reflector of Fuel Pin
!DO iz = myzb, myze
!  IF(Core%lFuelPlane(iz)) CYCLE
!  DO ixy = xyb, xye
!    IF(.NOT. Pin(ixy)%lfuel) CYCLE
!    rho = DenCool(iz, ixy) * epsm3
!    BulkTempC = Tcool(iz, ixy)
!    icel = Pin(ixy)%cell(iz)
!    FxrIdxSt = Pin(ixy)%FxrIdxSt
!    nLocalFxr = Cell(icel)%nFxr
!    DO i = 1, nLocalFxr
!      ifxr = FxrIdxSt + i - 1
!      Fxr(ifxr, iz)%temp = BulkTempC + CKELVIN
!      Fxr(ifxr, iz)%rho = rho
!    END DO
!  END DO 
!END DO

! Fuel Cell
tmvolsum = 0
denvolsum = 0
tdvolsum = 0
ncvolsum = 0
!PRINT'(20f8.2)', ThInfo%FuelTH(1)%tfuel(1:ThVar%npr4, 2)
DO iz = myzb, myze
  IF(.NOT. Core%lFuelPlane(iz)) CYCLE
  hz = Core%hz(iz)
  tmsum = 0
  densum = 0
  tdsum = 0
  ncell = 0
  DO ixy = xyb, xye
    IF(.NOT. Pin(ixy)%lfuel) CYCLE
    rho = DenCool(iz, ixy) * epsm3
    BulkTempC = Tcool(iz, ixy)

    DopTempC = (1-alpha) * ThInfo%FuelTH(ixy)%tfuel(1,iz) + alpha * ThInfo%FuelTH(ixy)%tfuel(nsurf,iz)
!    DopTempC = 287
    
    icel = Pin(ixy)%cell(iz)
    FxrIdxSt = Pin(ixy)%FxrIdxSt
    nLocalFxr = Cell(icel)%nFxr
    
    ncell = ncell + 1
    tdsum = tdsum + DopTempC
    densum = densum + rho
    tmsum = tmsum + BulkTempC

    DO i = 1, nLocalFxr
      ifxr = FxrIdxSt + i - 1
      Fxr(ifxr, iz)%temp = BulkTempC + CKELVIN
      Fxr(ifxr, iz)%Doptemp = DopTempC + CKELVIN
      Fxr(ifxr, iz)%rho = rho
    END DO 
  END DO 
  tdvolsum = tdvolsum + tdsum * hz
  tmvolsum = tmvolsum + tmsum * hz
  denvolsum = denvolsum + densum * hz
  ncvolsum = ncvolsum + real(ncell) * hz
END DO 

#ifdef MPI_ENV
Buf0 =(/tdvolsum, denvolsum, tmvolsum, ncvolsum/)
CALL MPI_ALLREDUCE(Buf0, Buf, 4, MPI_DOUBLE_PRECISION, MPI_SUM, PE%MPI_CMFD_COMM, ierr)
tdvolsum = Buf(1)
denvolsum = Buf(2)
tmvolsum = Buf(3)
ncvolsum = Buf(4)
#endif
WRITE(mesg, '(2x,3(A12,es12.4))') 'Avg T.Dop:', tdvolsum/ncvolsum, 'Avg Rho:', denvolsum/ncvolsum, 'Avg T.mod:', Tmvolsum/ncvolsum
IF (PE%MASTER) CALL message(io8, .TRUE., .TRUE., mesg)

!DO iz = myzb, myze
!  IF(.NOT. Core%lFuelPlane(iz)) CYCLE
!  DO ixy = xyb, xye
!    IF(.NOT. Pin(ixy)%lfuel) CYCLE
!    
!    icel = Pin(ixy)%cell(iz)
!    FxrIdxSt = Pin(ixy)%FxrIdxSt
!    nLocalFxr = Cell(icel)%nFxr
!
!    DO i = 1, nLocalFxr
!      ifxr = FxrIdxSt + i - 1
!      !Fxr(ifxr, iz)%temp = tmvolsum/ncvolsum + CKELVIN
!      Fxr(ifxr, iz)%Doptemp = tdvolsum/ncvolsum + CKELVIN
!      !Fxr(ifxr, iz)%rho = denvolsum/ncvolsum
!    END DO 
!  END DO 
!END DO 

! Radial Reflector
!DO iz = myzb, myze
!  DO ixy = xyb, xye
!    IF(Pin(ixy)%lFuel) CYCLE
!    rho = DenCool(iz, ixy) * epsm3
!    BulkTempC = Tcool(iz, ixy)
!    icel = Pin(ixy)%cell(iz)
!    FxrIdxSt = Pin(ixy)%FxrIdxSt
!    nLocalFxr = Cell(icel)%nFxr
!    DO i = 1, nLocalFxr
!      ifxr = FxrIdxSt + i - 1
!      Fxr(ifxr, iz)%temp = BulkTempC + CKELVIN
!      Fxr(ifxr, iz)%rho = rho
!    END DO
!  END DO 
!END DO 

END SUBROUTINE

SUBROUTINE SaveFxrTemp(Fxr, ithstep, nFxr, myzb, myze)
USE TYPEDEF,    ONLY : FxrInfo_Type
USE TH_Mod,     ONLY : FxrTempSave
IMPLICIT NONE
TYPE(FxrInfo_Type), POINTER :: Fxr(:,:)
INTEGER :: ithstep, nFxr, myzb, myze

INTEGER ::  ifxr, iz

!$OMP PARALLEL 
DO iz = myzb, myze
  !$OMP DO SCHEDULE(GUIDED)
  DO ifxr = 1, nfxr
    IF(Fxr(ifxr, iz)%niso .EQ. 0) CYCLE
    FxrTempSave(ifxr, iz, ithstep) = Fxr(ifxr, iz)%temp
  END DO
  !$OMP END DO
END DO
!$OMP END PARALLEL

END SUBROUTINE
SUBROUTINE RecoverFxrTemp(Fxr, ithstep, nFxr, myzb, myze)
USE TYPEDEF,    ONLY : FxrInfo_Type
USE TH_Mod,     ONLY : FxrTempSave
IMPLICIT NONE
TYPE(FxrInfo_Type), POINTER :: Fxr(:,:)
INTEGER :: ithstep, nFxr, myzb, myze

INTEGER ::  ifxr, iz

!$OMP PARALLEL 
DO iz = myzb, myze
  !$OMP DO SCHEDULE(GUIDED)
  DO ifxr = 1, nfxr
    IF(Fxr(ifxr, iz)%niso .EQ. 0) CYCLE
    Fxr(ifxr, iz)%temp = FxrTempSave(ifxr, iz, ithstep)
  END DO
  !$OMP END DO
END DO
!$OMP END PARALLEL

END SUBROUTINE