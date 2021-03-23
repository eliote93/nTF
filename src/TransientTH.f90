#include <defines.h>
  
SUBROUTINE TransientTH(Core, CmInfo, FmInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE, delt_inp)
USE PARAM
USE TYPEDEF,                    ONLY : CoreInfo_Type,              FmInfo_Type,               CmInfo_Type,          &
                                       ThInfo_Type,                GroupInfo_Type,            TranCntl_Type,        &
                                       PE_TYPE,                                                                     &
                                       FuelTH_Type,                CoolantTH_Type,            Pin_Type,             &
                                       FxrInfo_Type,               PinXS_Type
USE CNTL,                       ONLY : nTracerCntl_Type
USE TH_Mod,                     ONLY : ThOpt,                      ThVar,                                           &
                                       CalRelPower,                SetPwShape,                tfcaltr,              &
                                       TranRadFuelCondFDM,         Grp_RelPower,              hGapArray, Tfcaltr_pwshape
USE FuelProperty_Mod,           ONLY : fhtcoef,                    INTGAPCOND
USE timer,                      ONLY : nTracer_dclock,             TimeChk
USE FILES,                      ONLY : IO8
USE IOUTIL,                     ONLY : message
USE TIMER,                      ONLY : nTracer_dclock,             TimeChk
#ifdef MPI_ENV
USE MPIComm_Mod,                ONLY : MPI_SYNC
#endif
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(CmInfo_Type) :: CmInfo
TYPE(FmInfo_Type) :: FmInfo
TYPE(ThInfo_Type) :: ThInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
REAL, OPTIONAL :: delt_inp

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Fxrinfo_type), POINTER :: Fxr(:,:)
TYPE(PinXS_Type), POINTER :: PinXS(:, :)
TYPE(FuelTh_Type), POINTER :: FuelTH(:)      
TYPE(CoolantTH_Type), POINTER :: CoolantTH(:)
REAL, POINTER :: RelPower(:, :)
REAL, POINTER :: PwShapeAsy(:,:,:), PwShapedAsy(:,:,:)
REAL, POINTER :: PwShape(:, :), PwShaped(:, :)

INTEGER :: nowstep
INTEGER :: ng, nxy, nzth, nr, nsurf, nz, nchannel, maxnxy
INTEGER :: nasy, iasy, iasytype, ixy0, igrp
INTEGER :: Comm
INTEGER :: ixy, iz, iter, ixya, ichtyp

REAL :: qeff_grp(4), qeffd_grp(4), vol_grp(4)
REAL :: Tbeg, Tend
REAL :: dz, qprime, qf, qfd, qflux, qeff, qeffd
REAL :: Tdopld, TdoplMax, TdoplMax1
REAL :: TBulk, Tsurf, Tout, ToutAvg, Tfcl_max, Tfavg_max, htcoeff , twall_max, tmod_max, pwmax
REAL :: TBulkd, htcoeffd
REAL :: Deq, Fracdf, afp
REAL :: Fracdc, qc, acf
REAL :: PEXIT, BTEMP
REAL :: hgapavg, nhgap

LOGICAL :: Master, Slave

IF(nTracerCntl%libtyp .NE. 11) THEN 
IF(nTRACERCntl%lBenchXs) RETURN
END IF

IF(PRESENT(delt_inp)) THEN
  ThVar%DelT = delt_inp
ELSE
  NowStep = TranCntl%NowStep
  ThVar%DelT = TranCntl%DelT(NowStep)
END IF
Pin => Core%Pin;                  Fxr => FmInfo%FXR
PinXs => CmInfo%PinXS;            RelPower => ThInfo%RelPower
FuelTH => ThInfo%FuelTH;          CoolantTH => ThInfo%CoolantTH

Tbeg = nTracer_dclock(false, false)

IF(nTracerCntl%ThCh_mod .GE. 1 .AND. .NOT. nTracerCntl%lthch_tf) THEN 
  maxnxy = 0
  DO iasy = 1, Core%nAsyType
    maxnxy = MAX(maxnxy, Core%AsyInfo(iasy)%nxy)
  END DO 
  ALLOCATE(PwShapeAsy(ThVar%npr5, ThVar%nzth, maxnxy), PwShapedAsy(ThVar%npr5, ThVar%nzth, maxnxy))
ELSE
  ALLOCATE(PwShape(ThVar%npr5, ThVar%nzth), PwShaped(ThVar%npr5, ThVar%nzth))
END IF


Master = PE%CmfdMaster; Slave = .NOT. Master
Comm = PE%MPI_CMFD_COMM
IF(.NOT. TranCntl%ladptt) THEN
  WRITE(MESG, '(A)') 'Performing T/H Calculation...'
  IF(Master) CALL MESSAGE(io8,TRUE,TRUE,'Performing T/H Calculation...')
END IF

PEXIT = ThInfo%PExit; Btemp = ThVar%BoilingTemp  
nxy = Core%nxy; nzth = ThVar%nzth; ng = GroupInfo%ng 

nr = ThVar%npr5;nsurf = Thvar%npr4; nz = ThVar%nzth

fracdf = THVar%fracdf; DEQ =ThVar%Deq; afp = ThVar%afp
fracdc = THVar%fracdc; acf = ThVar%acf

CALL CalRelPower(Core, CmInfo, RelPower, GroupInfo%ng, nTracerCntl, PE, .TRUE.)
pwmax = 0
DO ixy = 1, nxy
  IF (.NOT. CoolantTH(ixy)%lfuel) CYCLE
  pwmax = max(pwmax, RelPower(0,ixy))
END DO 
! iasytype = Core%Asy(1)%AsyType
! DO ixy0 = 1, Core%AsyInfo(iasytype)%nxy
!   ixy = core%asy(1)%GlobalPinIdx(ixy0)
!   write(123, '(i6, f10.6)') ixy, RelPower(1,ixy)
! END DO 
! STOP

nchannel = 0; ToutAvg = 0; Tfcl_max = 0; Tfavg_max = 0; Twall_max = 0
tmod_max = 0; 
IF(nTracerCntl%ThCh_mod .GE. 1 .AND. .NOT. nTracerCntl%lthch_tf) THEN 
  nasy = Core%nxya
  DO iasy = 1, nasy
    iasytype = Core%Asy(iasy)%AsyType
    Deq = ThVar%deq
    acf = ThVar%acf
    IF(nTracerCntl%lThChConf) THEN 
      ichtyp = Core%THChMap(iasy)
      IF(ichtyp .NE. 0) THEN 
        Deq = ThVar%ThCh(ichtyp)%Deq 
        acf = ThVar%ThCh(ichtyp)%acf
      END IF
    END IF
    DO ixy0 = 1, Core%AsyInfo(iasytype)%nxy
      ixy = Core%Asy(iasy)%GlobalPinIdx(ixy0)
      CALL SetPwShape(PwShapeAsy(:,:,ixy0), Core, FmInfo%Fxr, FmInfo%Power, ixy, nzth, ThVar%npr5, Thvar, ThOpt, PE)
      CALL SetPwShape(PwshapedAsy(:,:,ixy0), Core, FmInfo%Fxr, FmInfo%TranPower, ixy, nzth, ThVar%npr5, Thvar, ThOpt, PE)
    END DO 
    DO iz = 1, ThVar%nzth
      DO iter = 1, 5
        qeff_grp = 0.
        qeffd_grp = 0.
        vol_grp = 0.
        DO ixy0 = 1, Core%AsyInfo(iasytype)%nxy
          igrp = Core%AsyInfo(iasytype)%ThChGrp(ixy0)
          ixy = Core%Asy(iasy)%GlobalPinIdx(ixy0)
          IF(.NOT. CoolantTH(ixy)%lfuel) CYCLE
          IF(iz .EQ.1 .AND. iter .EQ. 1) nchannel = nchannel + 1
          dz = ThVar%hz(iz)
          qprime = ThInfo%PowLv * ThInfo%PowLin * RelPower(iz, ixy)
          qf = fracdf * qprime / afp
          qfd = FuelTH(ixy)%qvold(iz)
          qc = fracdc * qprime / acf
          IF(ThOpt%hGapModel .EQ. 2) THEN
            ThOpt%hgap = INTGAPCOND(qprime)
          END IF
          Tbulk = CoolantTH(ixy)%Tcool(iz)
          Tbulkd = CoolantTH(ixy)%Tcoold(iz)
          htcoeff = fhtcoef(Tbulk, deq, CoolantTH(ixy)%RhoU(iz))                  !Heat Transfer Coefficient
          htcoeffd = fhtcoef(Tbulkd, deq, CoolantTH(ixy)%RhoUd(iz))               !Heat Transfer Coefficient
          ! CALL tfcaltr(FuelTH(ixy)%tfvol(1:nr-3,iz), FuelTH(ixy)%tfuel(1:nr, iz), FuelTH(ixy)%tfueld(1:nr, iz), Tbulk, htcoeff, &
          !              qf, qfd, nr, ThVar, ThOpt, hGapArray(iz, ixy))
          CALL Tfcaltr_pwshape(FuelTH(ixy)%tfvol(1:nr-3,iz), FuelTH(ixy)%tfuel(1:nr, iz), FuelTH(ixy)%tfueld(1:nr, iz), Tbulk, htcoeff, &
                       qf, qfd,PwShapeAsy(:,iz,ixy0), PwShapedAsy(:,iz,ixy0), nr, ThVar, ThOpt, hGapArray(iz, ixy))
          
          Tsurf = FuelTh(ixy)%tfuel(nsurf, iz)
          qflux = htcoeff * (Tsurf - Tbulk)
          qeff = qflux * ThVar%zetap + qc
          qeffd = CoolantTH(ixy)%qeffd(iz)
          qeff_grp(igrp) = qeff_grp(igrp) + Core%PinVol(ixy,iz) * qeff 
          qeffd_grp(igrp) = qeffd_grp(igrp) + Core%PinVol(ixy,iz) * qeffd 
          vol_grp(igrp) = vol_grp(igrp) + Core%PinVol(ixy,iz)
          !Maximum Fuel Temperature
          IF(iter .EQ. 5) THEN 
            Twall_max = MAX(FuelTh(ixy)%tfvol(ThVar%npr2, iz), Twall_max)
            Tfcl_max = MAX(FuelTh(ixy)%tfuel(1, iz), Tfcl_max)
            Tfavg_max = MAX(FuelTh(ixy)%tfuel(THVar%npr5, iz), Tfavg_max)
            !Save variables
            CoolantTH(ixy)%qeff(iz) = qeff
            FuelTH(ixy)%htcoef(iz) = htcoeff
            FuelTH(ixy)%qvol(iz) = qf
          END IF
        END DO 

        DO igrp = 1, Core%AsyInfo(iasytype)%nThChGrp
          IF(vol_grp(igrp) .LT. 1.0E-5) CYCLE 
          qeff_grp(igrp) = qeff_grp(igrp) / vol_grp(igrp)
          qeffd_grp(igrp) = qeffd_grp(igrp) / vol_grp(igrp)
        END DO 

        DO ixy0 = 1, Core%AsyInfo(iasytype)%nxy
          igrp = Core%AsyInfo(iasytype)%ThChGrp(ixy0)
          ixy = Core%Asy(iasy)%GlobalPinIdx(ixy0)
          IF(.NOT. CoolantTH(ixy)%lfuel) CYCLE
          IF(nTracerCntl%lThChConf) THEN
            CALL TransientCoolantTH_ThCh(Core, qeff_grp(igrp), qeffd_grp(igrp), ThInfo%PExit, CoolantTH(ixy), iz, THVar, ThOpt, ixy)
          ELSE
            CALL TransientCoolantTH(qeff_grp(igrp), qeffd_grp(igrp), ThInfo%PExit, CoolantTH(ixy), iz, THVar, ThOpt)
          END IF
          !Maximum Fuel Temperature
          IF(iter .EQ. 5) THEN 
            tmod_max = max(coolantth(ixy)%tcool(iz), tmod_max)
            IF(iz .EQ. ThVar%nzth) THEN 
              Tout = CoolantTH(ixy)%TCoolInOut(2, nz) 
              ToutAvg = ToutAvg + Tout
            END IF
          END IF
        END DO 
      END DO 
    END DO 
  END DO 
ELSE
IF(nTracerCntl%ThCh_mod .GE. 1) CALL Grp_RelPower(Core, CmInfo, RelPower, RelPower, ng, nTracerCntl, PE, .TRUE.)
DO ixy = 1, nxy
  IF(.NOT. CoolantTH(ixy)%lfuel) CYCLE
  Deq = ThVar%Deq
  acf = ThVar%acf
  IF(nTracerCntl%lThChConf) THEN 
    ixya = Pin(ixy)%iasy
    ichtyp = Core%ThChMap(ixya)
    IF(ichtyp .NE. 0) THEN
      Deq = ThVar%ThCh(ichtyp)%Deq
      acf = ThVar%ThCh(ichtyp)%acf
    END IF
  END IF
  
  nchannel = nchannel + 1
  !Set Power Shape
  CALL SetPwShape(Pwshape, Core, FmInfo%Fxr, FmInfo%Power, ixy, nzth, ThVar%npr5, Thvar, ThOpt, PE)
  CALL SetPwShape(Pwshaped, Core, FmInfo%Fxr, FmInfo%TranPower, ixy, nzth, ThVar%npr5, Thvar, ThOpt, PE)
  DO iz = 1, ThVar%nzth
    dz = ThVar%hz(iz) 
    qprime = ThInfo%PowLv * ThInfo%PowLin * RelPower(iz, ixy)
    qf = fracdf * qprime / afp; qfd = FuelTH(ixy)%qvold(iz)
    qc = fracdc * qprime / acf
    IF(ThOpt%hGapModel .EQ. 2) THEN
      ThOpt%hgap = INTGAPCOND(qprime)
    END IF

    !IF(PE%MASTER) PRINT*, ixy, iz, RelPower(iz, ixy), qf 
    DO iter = 1, 5 
      Tbulkd = CoolantTH(ixy)%Tcoold(iz); Tbulk = CoolantTH(ixy)%tcool(iz)                                         !Water Bulk Temperature
      htcoeff = fhtcoef(Tbulk, deq, CoolantTH(ixy)%RhoU(iz))                  !Heat Transfer Coefficient
      htcoeffd = fhtcoef(Tbulkd, deq, CoolantTH(ixy)%RhoUd(iz))               !Heat Transfer Coefficient
      IF(.NOT. nTracerCntl%LIHP) THEN
        Pwshape = 1; Pwshaped = 1
      ENDIF
      IF(nTracerCntl%lFuelCondFDM) THEN
        CALL TranRadFuelCondFDM(FuelTh(ixy)%tfvol(1:nr-3, iz), FuelTh(ixy)%tfuel(1:nr, iz), FuelTh(ixy)%tfueld(1:nr, iz),               &
                              Tbulk, Tbulkd, htcoeff, htcoeffd, qf, qfd, Pwshape(:, iz), Pwshaped(:, iz), nr, ThVar, ThOpt)
      ELSE
        !CALL tfcaltr(FuelTh(ixy)%tfvol(1:nr-3, iz), FuelTh(ixy)%tfuel(1:nr, iz), FuelTh(ixy)%tfueld(1:nr, iz), Tbulk, htcoeff, qf, qfd, nr, ThVar, ThOpt, hGapArray(iz, ixy))
        CALL Tfcaltr_pwshape(FuelTh(ixy)%tfvol(1:nr-3, iz), FuelTh(ixy)%tfuel(1:nr, iz), FuelTh(ixy)%tfueld(1:nr, iz), Tbulk, htcoeff, &
                      qf, qfd, PwShape(:,iz), PwShaped(:,iz), nr, ThVar, ThOpt, hGapArray(iz, ixy))

      ENDIF
      Tsurf =  FuelTh(ixy)%tfuel(nsurf, iz)
      qflux = htcoeff * (Tsurf - Tbulk)     !Heat Flux from Clad to Coolant
      qeff = qflux * ThVar%zetap + qc; qeffd = CoolantTH(ixy)%qeffd(iz)
      !IF(Core%lfuelplane(iz)) qeff = 0
      IF(nTracerCntl%lThChConf) THEN
        CALL TransientCoolantTH_ThCh(Core, qeff, qeffd, ThInfo%PExit, CoolantTH(ixy), iz, THVar, ThOpt, ixy)
      ELSE
        CALL TransientCoolantTH(qeff, qeffd, ThInfo%PExit, CoolantTH(ixy), iz, THVar, ThOpt)
      END IF
    ENDDO
      
    !Maximum Fuel Temperature
    Twall_max = MAX(FuelTh(ixy)%tfvol(ThVar%npr2, iz), Twall_max)
    Tfcl_max = MAX(FuelTh(ixy)%tfuel(1, iz), Tfcl_max)
    Tfavg_max = MAX(FuelTh(ixy)%tfuel(THVar%npr5, iz), Tfavg_max)
    tmod_max = max(coolantth(ixy)%tcool(iz), tmod_max)
    !Save variables
    CoolantTH(ixy)%qeff(iz) = qeff
    FuelTH(ixy)%htcoef(iz) = htcoeff
    FuelTH(ixy)%qvol(iz) = qf
  ENDDO
  Tout = CoolantTH(ixy)%TCoolInOut(2, nz) 
  ToutAvg = ToutAvg + Tout
ENDDO
END IF
!Assgin Temperature for Guide Tube and Assembly Gap
CALL SetGTnGapCoolTemp(Core, ThInfo, PE)
!Average Outlet Temperature
ToutAvg = ToutAvg / nChannel

!Doppler Temperature Update
TdoplMax = 0; nchannel = 0
DO ixy = 1, nxy
  IF(.NOT. CoolantTH(ixy)%lfuel) CYCLE
  nchannel = nchannel + 1
  DO iz = 1, nz
    Tdopld = ThInfo%Tdop(iz, ixy)
    ThInfo%Tdop(iz, ixY) = SQRT(CKELVIN + FuelTH(ixy)%tfueL(nr, iz))
    TdoplMax = max(TdoplMax, ABS(1._8 - tdopLd/THInfo%Tdop(iz, ixy)))
  ENDDO
ENDDO

ThInfo%Tfmax = Tfcl_max; ThInfo%TdopChg = Tdoplmax
ThInfo%Tfcl_max = Twall_max; ThInfo%Tfavg_max = Tfavg_max
!ThInfo%TModoutAvg = ToutAvg
ThInfo%TModoutAvg = tmod_max
ThInfo%pwpeakf = pwmax

hgapavg = 0.
nhgap = 0
DO ixy = 1, nxy
  IF(.NOT. CoolantTH(ixy)%lfuel) CYCLE
  DO iz = 1, nzth
    nhgap = nhgap + 1
    hgapavg = hgapavg + hGapArray(iz, ixy)
  END DO 
END DO 
hgapavg = hgapavg / nhgap

WRITE(mesg,'(2x,a,es14.6)') 'Average Gap Conductance: ', hgapavg 
IF(master) CALL message(io8,false,TRUE,mesg)

WRITE(mesg,601) Tfcl_max, Tfavg_max, toutavg,tdoplmax, FALSE
IF(master) CALL message(io8,false,TRUE,mesg)
601 format(4x,"Max Tf=",f6.1, "/" f6.1, "C, Avg Outlet Temp",f7.2,"C", ", Dop. Change=",1p,e9.2,l2,i3)
    
IF(nTracerCntl%ThCh_mod .GE. 1 .AND. .NOT. nTracerCntl%lthch_tf) THEN 
  DEALLOCATE(PwShapeAsy, PwShapedAsy)
ELSE
  DEALLOCATE(PwShape, PwShaped)
END IF
    
#ifdef MPI_ENV
CALL MPI_SYNC(Comm)
#endif

END SUBROUTINE

SUBROUTINE SaveTranTHsol(Core, ThInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,                     ONLY : CoreInfo_Type,           ThInfo_Type,             TranCntl_Type,      &
                                        PE_Type,                                                              &
                                        CoolantTH_Type,          FuelTh_Type
USE CNTL,                        ONLY : nTracerCntl_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(ThInfo_Type) :: ThInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

TYPE(CoolantTH_Type), POINTER :: CoolantTH(:)
TYPE(FuelTH_Type), POINTER :: FuelTH(:)


INTEGER :: nz, nxy
INTEGER :: ixy

CoolantTH => ThInfo%CoolantTH; FuelTH => ThInfo%FuelTH
nz = Core%nz; nxy = Core%nxy

DO ixy = 1, nxy
  IF(.NOT. CoolantTH(ixy)%lfuel) CYCLE
  CoolantTH(ixy)%TCoold(1:nz) = CoolantTH(ixy)%TCool(1:nz)
  CoolantTH(ixy)%DenCoold(1:nz) = CoolantTH(ixy)%DenCool(1:nz)
  CoolantTH(ixy)%hCoold(1:nz) = CoolantTH(ixy)%hCool(1:nz)
  CoolantTH(ixy)%rhoud(0:nz) = CoolantTH(ixy)%rhou(0:nz)
  CoolantTH(ixy)%rhohud(0:nz) = CoolantTH(ixy)%rhohu(0:nz)
  CoolantTH(ixy)%ud(0:nz) = CoolantTH(ixy)%u(0:nz)
  CoolantTH(ixy)%qeffd(1:nz) = CoolantTH(ixy)%qeff(1:nz) 
  FuelTH(ixy)%qvold(1:nz) = FuelTH(ixy)%qvol(1:nz)
  FuelTH(ixy)%tfueld(:, 1:nz) = FuelTH(ixy)%tfuel(:, 1:nz)
ENDDO

NULLIFY(CoolantTH, FuelTH)

END SUBROUTINE

SUBROUTINE RecoverTranTHsol(Core, ThInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,                     ONLY : CoreInfo_Type,           ThInfo_Type,             TranCntl_Type,      &
                                        PE_Type,                                                              &
                                        CoolantTH_Type,          FuelTh_Type
USE CNTL,                        ONLY : nTracerCntl_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(ThInfo_Type) :: ThInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

TYPE(CoolantTH_Type), POINTER :: CoolantTH(:)
TYPE(FuelTH_Type), POINTER :: FuelTH(:)


INTEGER :: nz, nxy
INTEGER :: ixy

CoolantTH => ThInfo%CoolantTH; FuelTH => ThInfo%FuelTH
nz = Core%nz; nxy = Core%nxy

DO ixy = 1, nxy
  IF(.NOT. CoolantTH(ixy)%lfuel) CYCLE
  CoolantTH(ixy)%TCool(1:nz) = CoolantTH(ixy)%TCoold(1:nz)
  CoolantTH(ixy)%DenCool(1:nz) = CoolantTH(ixy)%DenCoold(1:nz)
  CoolantTH(ixy)%hCool(1:nz) = CoolantTH(ixy)%hCoold(1:nz)
  CoolantTH(ixy)%rhou(0:nz) = CoolantTH(ixy)%rhoud(0:nz)
  CoolantTH(ixy)%rhohu(0:nz) = CoolantTH(ixy)%rhohud(0:nz)
  CoolantTH(ixy)%u(0:nz) = CoolantTH(ixy)%ud(0:nz)
  CoolantTH(ixy)%qeff(1:nz) = CoolantTH(ixy)%qeffd(1:nz) 
  FuelTH(ixy)%qvol(1:nz) = FuelTH(ixy)%qvold(1:nz)
  FuelTH(ixy)%tfuel(:, 1:nz) = FuelTH(ixy)%tfueld(:, 1:nz)
ENDDO
NULLIFY(CoolantTH, FuelTH)

END SUBROUTINE

SUBROUTINE SaveBaseTH(Core, THInfo)
USE TYPEDEF,      ONLY : CoreInfo_Type,     ThInfo_Type,    CoolantTH_Type,   FuelTH_Type
USE Th_Mod,       ONLY : CoolantTHSave,     FuelTHSave
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(ThInfo_Type) :: THInfo
TYPE(CoolantTH_Type), POINTER :: CoolantTH(:)
TYPE(FuelTH_Type), POINTER :: FuelTH(:)

INTEGER :: nz, nxy
INTEGER :: ixy

CoolantTH => ThInfo%CoolantTH; FuelTH => ThInfo%FuelTH
nz = Core%nz; nxy = Core%nxy
DO ixy = 1, nxy
  IF(.NOT. CoolantTH(ixy)%lfuel) CYCLE
  CoolantTHSave(ixy)%Tcool(1:nz) = CoolantTH(ixy)%Tcool(1:nz)
  CoolantTHSave(ixy)%Dencool(1:nz) = CoolantTH(ixy)%Dencool(1:nz)
  CoolantTHSave(ixy)%hcool(1:nz) = CoolantTH(ixy)%hcool(1:nz)
  CoolantTHSave(ixy)%rhou(0:nz) = CoolantTH(ixy)%rhou(0:nz)
  CoolantTHSave(ixy)%rhohu(0:nz) = CoolantTH(ixy)%rhohu(0:nz)
  CoolantTHSave(ixy)%u(0:nz) = CoolantTH(ixy)%u(0:nz)
  CoolantTHSave(ixy)%qeff(0:nz) = CoolantTH(ixy)%qeff(0:nz)
  FuelTHSave(ixy)%qvol(1:nz) = FuelTH(ixy)%qvol(1:nz)
  FuelTHSave(ixy)%tfuel(:,1:nz) = FuelTH(ixy)%tfuel(:,1:nz)

  CoolantTHSave(ixy)%Tcoold(1:nz) = CoolantTH(ixy)%Tcoold(1:nz)
  CoolantTHSave(ixy)%Dencoold(1:nz) = CoolantTH(ixy)%Dencoold(1:nz)
  CoolantTHSave(ixy)%hcoold(1:nz) = CoolantTH(ixy)%hcoold(1:nz)
  CoolantTHSave(ixy)%rhoud(0:nz) = CoolantTH(ixy)%rhoud(0:nz)
  CoolantTHSave(ixy)%rhohud(0:nz) = CoolantTH(ixy)%rhohud(0:nz)
  CoolantTHSave(ixy)%ud(0:nz) = CoolantTH(ixy)%ud(0:nz)
  CoolantTHSave(ixy)%qeffd(0:nz) = CoolantTH(ixy)%qeffd(0:nz)
  FuelTHSave(ixy)%qvold(1:nz) = FuelTH(ixy)%qvold(1:nz)
  FuelTHSave(ixy)%tfueld(:,1:nz) = FuelTH(ixy)%tfueld(:,1:nz)
END DO 
NULLIFY(CoolantTH, FuelTH)

END SUBROUTINE

SUBROUTINE RecoverBaseTH(Core, ThInfo)
USE TYPEDEF,      ONLY : CoreInfo_Type,     ThInfo_Type,    CoolantTH_Type,   FuelTH_Type
USE Th_Mod,       ONLY : CoolantTHSave,     FuelTHSave
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(ThInfo_Type) :: THInfo
TYPE(CoolantTH_Type), POINTER :: CoolantTH(:)
TYPE(FuelTH_Type), POINTER :: FuelTH(:)

INTEGER :: nz, nxy
INTEGER :: ixy

CoolantTH => ThInfo%CoolantTH; FuelTH => ThInfo%FuelTH
nz = Core%nz; nxy = Core%nxy
DO ixy = 1, nxy
  IF(.NOT. CoolantTH(ixy)%lfuel) CYCLE
  CoolantTH(ixy)%Tcool(1:nz) = CoolantTHSave(ixy)%Tcool(1:nz)
  CoolantTH(ixy)%Dencool(1:nz) = CoolantTHSave(ixy)%Dencool(1:nz)
  CoolantTH(ixy)%hcool(1:nz) = CoolantTHSave(ixy)%hcool(1:nz)
  CoolantTH(ixy)%rhou(0:nz) = CoolantTHSave(ixy)%rhou(0:nz)
  CoolantTH(ixy)%rhohu(0:nz) = CoolantTHSave(ixy)%rhohu(0:nz)
  CoolantTH(ixy)%u(0:nz) = CoolantTHSave(ixy)%u(0:nz)
  CoolantTH(ixy)%qeff(0:nz) = CoolantTHSave(ixy)%qeff(0:nz)
  FuelTH(ixy)%qvol(1:nz) = FuelTHSave(ixy)%qvol(1:nz)
  FuelTH(ixy)%tfuel(:,1:nz) = FuelTHSave(ixy)%tfuel(:,1:nz)

  CoolantTH(ixy)%Tcoold(1:nz) = CoolantTHSave(ixy)%Tcoold(1:nz)
  CoolantTH(ixy)%Dencoold(1:nz) = CoolantTHSave(ixy)%Dencoold(1:nz)
  CoolantTH(ixy)%hcoold(1:nz) = CoolantTHSave(ixy)%hcoold(1:nz)
  CoolantTH(ixy)%rhoud(0:nz) = CoolantTHSave(ixy)%rhoud(0:nz)
  CoolantTH(ixy)%rhohud(0:nz) = CoolantTHSave(ixy)%rhohud(0:nz)
  CoolantTH(ixy)%ud(0:nz) = CoolantTHSave(ixy)%ud(0:nz)
  CoolantTH(ixy)%qeffd(0:nz) = CoolantTHSave(ixy)%qeffd(0:nz)
  FuelTH(ixy)%qvold(1:nz) = FuelTHSave(ixy)%qvold(1:nz)
  FuelTH(ixy)%tfueld(:,1:nz) = FuelTHSave(ixy)%tfueld(:,1:nz)
END DO 
NULLIFY(CoolantTH, FuelTH)

END SUBROUTINE

SUBROUTINE UpdtTdopErr(ThInfo, Tdoperr, TDopSave, nxy, nz)
USE TYPEDEF,     ONLY : ThInfo_Type,    CoolantTH_Type
IMPLICIT NONE
TYPE(ThInfo_Type) :: ThInfo
REAL :: TDopErr
REAL :: TDopSave(nz, nxy)
INTEGER :: nxy, nz

TYPE(CoolantTH_Type),POINTER :: CoolantTH(:)
REAL :: tmperr
INTEGER :: ixy, iz

CoolantTH => ThInfo%CoolantTH
TDopErr = 0.
DO ixy = 1, nxy 
  IF(.NOT. CoolantTH(ixy)%lfuel) CYCLE
  DO iz = 1, nz
    tmperr = ABS(1._8 - TDopSave(iz,ixy)/THInfo%Tdop(iz,ixy))
    TDopErr = max(TDopErr, tmperr)
    TDopSave(iz,ixy) = ThInfo%TDop(iz,ixy)
  END DO 
END DO 

END SUBROUTINE
  
!
!SUBROUTINE TransientTH_Test(ThInfo, ThVar, ThOpt)
!USE PARAM
!USE TYPEDEF,    ONLY : ThInfo_Type,          ThVar_Type,         ThOpt_Type,     &
!                       FuelTh_Type,          CoolantTH_Type
!USE FuelProperty_Mod, ONLY : fhtcoef
!IMPLICIT NONE
!TYPE(THInfo_Type) :: ThInfo
!TYPE(ThVar_Type) :: ThVar
!TYPE(ThOpt_Type) :: ThOpt
!TYPE(FuelTh_Type), POINTER :: FuelTh
!TYPE(CoolantTH_Type), POINTER :: CoolantTH
!
!REAL, POINTER :: RelPower(:, :)
!
!INTEGER :: i, iz, j
!INTEGER :: nr, nsurf, nz
!REAL :: dz, qf, qfd, qprime, deq
!REAL :: afp, fracdf
!REAL :: htcoeff, Tbulk, Tsurf, qflux, qeff, qeffd, Tink
!
!TinK = ThInfo%Tin + CKELVIN
!
!ThVar%DelT = 0.01
!nr = Thvar%npr5; nsurf = Thvar%npr4; nz = ThVar%nzth
!Deq = ThVar%Deq
!fracdf = THVar%fracdf
!DEQ =ThVar%Deq; afp = ThVar%afp
!FuelTH => THInfo%FuelTH(1); CoolantTH => THInfo%CoolantTH(1)
!RelPower => ThInfo%RelPower
!
!CALL InitCoolantVar(Tink, ThInfo%PExit, ThInfo%MdotFa, ThVar%acf, ThInfo%CoolantTh, ThVar%npr , 1 , THVar%nzth)
!CALL  InitTransientTH(THInfo, ThVar, 1, ThVar%nzth)
!CoolantTH%DenCOOL(:) = CoolantTH%DenCOOLD(:)
!CoolantTH => ThInfo%CoolantTH(1)
!CoolantTH%Tcool(:) = Tink - CKELVIN
!CoolantTH%TcoolInout(:, :) = Tink - CKELVIN
!FuelTh%tfuel =  Tink - CKELVIN
!FuelTH%qvold = 0; CoolantTH%qeffd = 0 
!
!ThInfo%RelPower(1,1)	= 0.000000000000000     
!ThInfo%RelPower(2,1)	= 0.000000000000000     
!ThInfo%RelPower(3,1)	= 0.365526675209038     
!ThInfo%RelPower(4,1)	= 0.617572551458803     
!ThInfo%RelPower(5,1)	= 0.846207965299850     
!ThInfo%RelPower(6,1)	= 1.02082671617635	    
!ThInfo%RelPower(7,1)	= 1.14795070062357	    
!ThInfo%RelPower(8,1)	= 1.23571192682966	    
!ThInfo%RelPower(9,1)	= 1.29190873137794	    
!ThInfo%RelPower(10,1)	= 1.32317321156984	    
!ThInfo%RelPower(11,1)	= 1.33468705026985	    
!ThInfo%RelPower(12,1)	= 1.33012948121759	    
!ThInfo%RelPower(13,1)	= 1.31172680025683	    
!ThInfo%RelPower(14,1)	= 1.28032551835021	    
!ThInfo%RelPower(15,1)	= 1.23545633287335	    
!ThInfo%RelPower(16,1)	= 1.17537972986876	    
!ThInfo%RelPower(17,1)	= 1.09713025071640	    
!ThInfo%RelPower(18,1)	= 0.996590036280530     
!ThInfo%RelPower(19,1)	= 0.868720693856457     
!ThInfo%RelPower(20,1)	= 0.708133686689330     
!ThInfo%RelPower(21,1)	= 0.510443211582325     
!ThInfo%RelPower(22,1)	= 0.302398729493319     
!ThInfo%RelPower(23,1)	= 0.000000000000000     
!ThInfo%RelPower(24,1)	= 0.000000000000000     
!
!DO i = 1, 1000
!  DO iz = 1, ThVar%nzth
!    dz = ThVar%hz(iz)
!    qprime = ThInfo%PowLv * ThInfo%PowLin * RelPower(iz, 1)
!    qf = fracdf * qprime / afp; qfd = FuelTH%qvold(iz)
!    DO j= 1, 5
!      Tbulk = CoolantTH%tcool(iz)
!      htcoeff = fhtcoef(Tbulk, deq, CoolantTH%RhoU(iz))
!      CALL tfcaltr(FuelTh%tfvol(1:nr-3, iz), FuelTh%tfuel(1:nr, iz),       &
!                   Tbulk, htcoeff, qf, qfd, nr, ThVar, ThOpt)
!      Tsurf =  FuelTh%tfuel(nsurf, iz)
!      qflux = htcoeff * (Tsurf - Tbulk)
!      qeff = qflux * ThVar%zetap
!      qeffd = CoolantTH%qeffd(iz)
!      CALL TransientCoolantTH(qeff, qeffd, ThInfo%PExit, CoolantTH, iz, THVar, ThOpt)
!    ENDDO
!    CoolantTH%qeff(iz) = qeff
!    FuelTH%htcoef(iz) = htcoeff
!    FuelTH%qvol(iz) = qf
!  ENDDO
!  CoolantTH%DenCoold(1:nz) = CoolantTH%DenCool(1:nz)
!  CoolantTH%hCoold(1:nz) = CoolantTH%hCool(1:nz)
!  CoolantTH%rhoud(0:nz) = CoolantTH%rhou(0:nz)
!  CoolantTH%rhohud(0:nz) = CoolantTH%rhohu(0:nz)
!  CoolantTH%ud(0:nz) = CoolantTH%u(0:nz)
!  CoolantTH%qeffd(1:nz) = CoolantTH%qeff(1:nz) 
!  FuelTH%qvold = FuelTH%qvol
!  IF(mod(i, 25) .EQ. 0) THEN
!    WRITE(99, '(F10.4, A5, 100F15.5)') ThVar%DelT * i, 'sec', (CoolantTH%Tcool(iz) - ThInfo%Tin, iz = 1, nz)
!    WRITE(100, '(F10.4, A5, 100F15.5)') ThVar%DelT * i, 'sec', (FuelTh%tfuel(1, iz), iz = 1, nz)
!  ENDIF
!ENDDO
!
!END SUBROUTINE
!
!SUBROUTINE TransientTH_Test0(ThInfo, ThVar, ThOpt)
!USE PARAM
!USE TYPEDEF,    ONLY : ThInfo_Type,  ThVar_Type,   ThOpt_Type,   FuelTh_Type
!USE FuelProperty_Mod, ONLY : fhtcoef
!IMPLICIT NONE
!TYPE(THInfo_Type) :: ThInfo
!TYPE(ThVar_Type) :: ThVar
!TYPE(ThOpt_Type) :: ThOpt
!TYPE(FuelTh_Type), POINTER :: FuelTh
!
!REAL, POINTER :: hflux(:), qvol(:), tcool(:), htcoef(:), tfuel(:, :), tfvol(:, :), RhoU(:)
!
!REAL :: TinK, qf, qfd, DEQ, afp
!INTEGER :: iz, nz, nr, i, ir
!ThVar%DelT = 0.005
!
!TinK = ThInfo%Tin + CKELVIN
!nr = Thvar%npr5
!FuelTH => THInfo%FuelTH(1)
!DEQ =ThVar%Deq; afp = ThVar%afp
!hflux => FuelTH%hflux; qvol => FuelTh%qvol
!tcool => FuelTH%tcool; htcoef => FuelTh%htcoef
!tfuel => FuelTh%tfuel; tfvol => FuelTH%tfvol
!RhoU => ThInfo%CoolantTH(1)%RhoU
!!
!!CALL InitCoolantVar(Tink, ThInfo%PExit, ThInfo%MdotFa, ThVar%acf, ThInfo%CoolantTh, ThVar%npr , 1 , THVar%nzth)
!!CALL  InitTransientTH(THInfo, ThVar, 1, ThVar%nzth)
!!CoolantTH%DenCOOL(:) = CoolantTH%DenCOOLD(:)
!!CoolantTH => ThInfo%CoolantTH(1)
!!CoolantTH%Tcool(:) = Tink - CKELVIN
!!CoolantTH%TcoolInout(:, :) = Tink - CKELVIN
!!
!iz = 3 
!qf = qvol(iz)
!qfd = qf
!htcoef(iz) = fhtcoef(tcool(iz), deq, RhoU(iz))
!Tfuel(:, iz) = THinfo%Tin
!DO i = 1, 5000
!  CALL tfcaltr(tfvol(1:nr-3, iz), tfuel(1:nr, iz), tcool(iz), htcoef(iz), qf, qfd, nr, ThVar, ThOpt)
!  IF(mod(i, 200) .EQ. 0) THEN
!    WRITE(99, '(F10.4, A5, 100F15.5)') ThVar%DelT * i, 'sec', ( tfuel(ir, iz), ir = 1, nr)
!  ENDIF
!ENDDO
!stop
!!CALL Tfcaltr(tfvol, Tfuel,Tcool, htCoef, Qf, Qfd, nr, ThVar, ThOpt)
!END SUBROUTINE
!  
!SUBROUTINE TransientTH_Test1(ThInfo, ThVar, ThOpt)
!USE PARAM
!USE TYPEDEF,    ONLY : ThInfo_Type,  ThVar_Type,   ThOpt_Type,   CoolantTH_Type
!IMPLICIT NONE
!TYPE(THInfo_Type) :: ThInfo
!TYPE(ThVar_Type) :: ThVar
!TYPE(ThOpt_Type) :: ThOpt
!TYPE(CoolantTH_Type), POINTER :: CoolantTH
!REAL :: TinK, qeff, qeffd
!INTEGER :: iz, nz , i 
!ThVar%DelT = 0.01
!
!TinK = ThInfo%Tin + CKELVIN
!nz = THvar%nzth
!CALL InitCoolantVar(Tink, ThInfo%PExit, ThInfo%MdotFa, ThVar%acf, ThInfo%CoolantTh, ThVar%npr , 1 , THVar%nzth)
!CALL  InitTransientTH(THInfo, ThVar, 1, ThVar%nzth)
!THInfo%CoolantTH(1)%DenCOOL(:) = THInfo%CoolantTH(1)%DenCOOLD(:)
!CoolantTH => ThInfo%CoolantTH(1)
!CoolantTH%Tcool(:) = Tink - CKELVIN
!CoolantTH%TcoolInout(:, :) = Tink - CKELVIN
!DO i = 1, 500
!  DO iz = 1, ThVar%nzth
!    qeff = CoolantTH%qeff(iz)
!    qeffd =  CoolantTH%qeff(iz) 
!    CALL TransientCoolantTH(qeff, qeffd, ThInfo%PExit, CoolantTH, iz, THVar, ThOpt)
!    CONTINUE  
!  ENDDO
!  
!  CoolantTH%DenCoold(1:nz) = CoolantTH%DenCool(1:nz)
!  CoolantTH%hCoold(1:nz) = CoolantTH%hCool(1:nz)
!  CoolantTH%rhoud(0:nz) = CoolantTH%rhou(0:nz)
!  CoolantTH%rhohud(0:nz) = CoolantTH%rhohu(0:nz)
!  CoolantTH%ud(0:nz) = CoolantTH%u(0:nz)
!  !CoolantTH%qeffd(1:nz) = CoolantTH%qeff(1:nz)
!  CONTINUE
!  IF(mod(i, 25) .EQ. 0) THEN
!    WRITE(99, '(F10.4, A5, 100F15.5)') ThVar%DelT * i, 'sec', (CoolantTH%Tcool(iz) - ThInfo%Tin, iz = 1, nz)
!  ENDIF
!ENDDO
!STOP
!END SUBROUTINE
  
SUBROUTINE TransientCoolantTH(qeff, qeffd, Pexit, CoolantTH, iz, ThVar, ThOpt)
USE PARAM
USE TYPEDEF,          ONLY : CoolantTH_Type,       ThVar_Type,       ThOpt_Type
USE SteamTBL_mod,     ONLY : steamtbl
IMPLICIT NONE
TYPE(CoolantTH_Type) :: CoolantTH                            ! One Channel Coolant TH in 
TYPE(ThVar_Type) :: ThVar
TYPE(ThOPT_Type) :: ThOpt                
REAL :: qeff, qeffd, PEXIT 
INTEGER :: iz

REAL, POINTER :: hz(:)

REAL :: cetac, cetacb, cetacr                                ! Ceta method
REAL :: qc, tout, rhoout, uout, uout0, uoutd, uavg
REAL :: td, hd, rhohd, rhod, rhouind, rhououtd, rhohuind, rhohuoutd              !TH variable at the current time step
REAL :: t, h, rho, rhoh, rhouin, rhouout, rhohuin, rhohuout         !TH variable at the next time step
REAL :: acf, afp, xi, zeta, zetap, FracDC, Fracdf

REAL :: delrhoud, delrhohud, hin, hout, houtd, hind
REAL :: tmpterm, sqterm, dzcetadt, delh, delh1, delh2, delh3
REAL :: alpha, beta, gamma, delta, a, b, c, biga, bigb

REAL :: wt, wh, wrho, wvin, wxin, wbetain, wkapain, wcpin
REAL :: hr, hl, rhor, rhol, dh, drhodh, drhohdh
REAL, PARAMETER :: eps = 1.e-4_8

cetac = 0.5

cetacb = 1._8 - cetac; cetacr = cetacb / cetac

acf = ThVar%acf; afp = ThVar%afp
xi = ThVar%Xi; zeta = ThVar%zeta
zetap = ThVar%zetap
FracDc = ThVar%FracDC; FracDf = ThVar%FracDf; 

hz => ThVar%hz

qc = qeff + cetacr * qeffd
dzcetadt = hz(iz) / (cetac * ThVar%DelT)
!set current time step variable

h = CoolantTH%hcool(iz); rho = CoolantTH%DenCool(iz)
rhouout = CoolantTH%rhou(iz); 
rhouin = CoolantTH%rhou(iz-1); rhohuin = CoolantTH%rhohu(iz-1)
uout0 = CoolantTH%u(iz); uoutd = CoolantTH%ud(iz)

hd = CoolantTH%hcoold(iz); rhod = CoolantTH%DenCoold(iz)
rhououtd = CoolantTH%rhoud(iz); rhouind = CoolantTH%rhoud(iz-1)
rhohuoutd = CoolantTH%rhohud(iz); rhohuind = CoolantTH%rhohud(iz-1)
Td = CoolantTH%Tcool(iz)
hin = rhohuin / rhouin
! Calculate Derivaties
hr = (1._8+eps)*hin; hl = (1._8-eps)*hin
dh =hr - hl
wh = hr
CALL steamtbl(FALSE,pexit,wt,wh,wrho,wvin,wxin,wbetain,wkapain,wcpin)
rhor = wrho
wh = hl
CALL steamtbl(FALSE,pexit,wt,wh,wrho,wvin,wxin,wbetain,wkapain,wcpin)
rhol = wrho
drhodh = (rhor - rhol) / dh
drhohdh = (rhor * hr - rhol * hl) / dh

delrhoud =   (rhououtd - rhouind)
delrhohud = (rhohuoutd - rhohuind) 
hin = rhohuin / rhouin
alpha = - dzcetadt * drhodh
beta = rhouin - cetacr * delrhoud
tmpterm = 2._8 * hd - hin
gamma = alpha * tmpterm + 2 * beta
delta = beta * tmpterm

a = 2._8 * alpha
b = gamma + dzcetadt * drhohdh
c = cetacr * delrhohud + delta - rhohuin - hz(iz) * qc
sqterm = sqrt(b*b - 4._8 * a * c)
delh1 = half*(-b+sqterm)/a; delh2=half*(-b-sqterm)/a

!Check enthaly change
houtd = rhohuoutd / rhououtd
biga = dzcetadt*drhohdh+alpha*houtd
bigb= hz(iz) * qc - rhouin * houtd + rhohuin + cetacr*(delrhoud*houtd-delrhohud)
delh3 = bigb / biga

!Node Average Enhalpy Calculation
IF(abs(delh1-delh3) .lt. abs(delh2-delh3)) THEN
  delh = delh1
ELSE
  delh = delh2
ENDIF

h = hd + delh
wh = h
CALL steamtbl(FALSE,pexit,wt,wh,wrho,wvin,wxin,wbetain,wkapain,wcpin)
t = wt - CKELVIN;   rho = wrho

rhohd = rhod * hd; rhoh = rho * h
rhouout = rhouin - dzcetadt * (rho - rhod) - cetacr * delrhoud
rhohuout = rhohuin - dzcetadt * (rhoh - rhohd) - cetacr * delrhohud + hz(iz) * qc
hout = rhohuout / rhouout

wh = hout
CALL steamtbl(FALSE,pexit,wt,wh,wrho,wvin,wxin,wbetain,wkapain,wcpin)
tout = wt - CKELVIN; rhoout = wrho
uout = rhouout / rhoout; uavg = cetac * (uout + cetacr * uout0)

IF(((uoutd-uout0) * (uout-uout0) .LT. -1.e-8*uavg * uavg) .AND. cetac .LT. 1.0_8) THEN
  uout = uavg
  rhouout = rhoout * uout
  rhohuout = rhouout * hout
ENDIF

CoolantTH%Tcool(iz) = T; CoolantTH%DenCool(iz) = rho
CoolantTH%hcool(iz) = h; CoolantTH%u(iz) = uout
CoolantTH%TCoolInOut(2, iz) = Tout; CoolantTH%TCoolInOut(1, iz+1) = Tout
CoolantTH%Rhou(iz) = rhouout; CoolantTH%rhohu(iz) = rhohuout

END SUBROUTINE

SUBROUTINE TransientCoolantTH_ThCh(Core, qeff, qeffd, Pexit, CoolantTH, iz, ThVar, ThOpt, ixy)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type,   CoolantTH_Type,       ThVar_Type,       ThOpt_Type
USE SteamTBL_mod,     ONLY : steamtbl
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(CoolantTH_Type) :: CoolantTH                            ! One Channel Coolant TH in 
TYPE(ThVar_Type) :: ThVar
TYPE(ThOPT_Type) :: ThOpt                
REAL :: qeff, qeffd, PEXIT 
INTEGER :: iz, ixy

REAL, POINTER :: hz(:)

REAL :: cetac, cetacb, cetacr                                ! Ceta method
REAL :: qc, tout, rhoout, uout, uout0, uoutd, uavg
REAL :: td, hd, rhohd, rhod, rhouind, rhououtd, rhohuind, rhohuoutd              !TH variable at the current time step
REAL :: t, h, rho, rhoh, rhouin, rhouout, rhohuin, rhohuout         !TH variable at the next time step
REAL :: acf, afp, xi, zeta, zetap, FracDC, Fracdf

REAL :: delrhoud, delrhohud, hin, hout, houtd, hind
REAL :: tmpterm, sqterm, dzcetadt, delh, delh1, delh2, delh3
REAL :: alpha, beta, gamma, delta, a, b, c, biga, bigb

REAL :: wt, wh, wrho, wvin, wxin, wbetain, wkapain, wcpin
REAL :: hr, hl, rhor, rhol, dh, drhodh, drhohdh
REAL, PARAMETER :: eps = 1.e-4_8
INTEGER :: ixya, ichtyp

cetac = 0.5
cetacb = 1._8 - cetac; cetacr = cetacb / cetac

ixya = Core%Pin(ixy)%iasy
ichtyp = Core%ThChMap(ixya)

IF(ichtyp .EQ. 0) THEN 
  acf = ThVar%acf
  xi = ThVar%xi
  zetap = ThVar%zetap
ELSE
  acf = ThVar%ThCh(ichtyp)%acf
  xi = ThVar%ThCh(ichtyp)%xi
  zetap = ThVar%ThCh(ichtyp)%zetap
END IF

afp = ThVar%afp
zeta = ThVar%zeta
FracDc = ThVar%FracDC; FracDf = ThVar%FracDf; 

hz => ThVar%hz

qc = qeff + cetacr * qeffd
dzcetadt = hz(iz) / (cetac * ThVar%DelT)
!set current time step variable

h = CoolantTH%hcool(iz); rho = CoolantTH%DenCool(iz)
rhouout = CoolantTH%rhou(iz); 
rhouin = CoolantTH%rhou(iz-1); rhohuin = CoolantTH%rhohu(iz-1)
uout0 = CoolantTH%u(iz); uoutd = CoolantTH%ud(iz)

hd = CoolantTH%hcoold(iz); rhod = CoolantTH%DenCoold(iz)
rhououtd = CoolantTH%rhoud(iz); rhouind = CoolantTH%rhoud(iz-1)
rhohuoutd = CoolantTH%rhohud(iz); rhohuind = CoolantTH%rhohud(iz-1)
Td = CoolantTH%Tcool(iz)
hin = rhohuin / rhouin
! Calculate Derivaties
hr = (1._8+eps)*hin; hl = (1._8-eps)*hin
dh =hr - hl
wh = hr
CALL steamtbl(FALSE,pexit,wt,wh,wrho,wvin,wxin,wbetain,wkapain,wcpin)
rhor = wrho
wh = hl
CALL steamtbl(FALSE,pexit,wt,wh,wrho,wvin,wxin,wbetain,wkapain,wcpin)
rhol = wrho
drhodh = (rhor - rhol) / dh
drhohdh = (rhor * hr - rhol * hl) / dh

delrhoud =   (rhououtd - rhouind)
delrhohud = (rhohuoutd - rhohuind) 
hin = rhohuin / rhouin
alpha = - dzcetadt * drhodh
beta = rhouin - cetacr * delrhoud
tmpterm = 2._8 * hd - hin
gamma = alpha * tmpterm + 2 * beta
delta = beta * tmpterm

a = 2._8 * alpha
b = gamma + dzcetadt * drhohdh
c = cetacr * delrhohud + delta - rhohuin - hz(iz) * qc
sqterm = sqrt(b*b - 4._8 * a * c)
delh1 = half*(-b+sqterm)/a; delh2=half*(-b-sqterm)/a

!Check enthaly change
houtd = rhohuoutd / rhououtd
biga = dzcetadt*drhohdh+alpha*houtd
bigb= hz(iz) * qc - rhouin * houtd + rhohuin + cetacr*(delrhoud*houtd-delrhohud)
delh3 = bigb / biga

!Node Average Enhalpy Calculation
IF(abs(delh1-delh3) .lt. abs(delh2-delh3)) THEN
  delh = delh1
ELSE
  delh = delh2
ENDIF

h = hd + delh
wh = h
CALL steamtbl(FALSE,pexit,wt,wh,wrho,wvin,wxin,wbetain,wkapain,wcpin)
t = wt - CKELVIN;   rho = wrho

rhohd = rhod * hd; rhoh = rho * h
rhouout = rhouin - dzcetadt * (rho - rhod) - cetacr * delrhoud
rhohuout = rhohuin - dzcetadt * (rhoh - rhohd) - cetacr * delrhohud + hz(iz) * qc
hout = rhohuout / rhouout

wh = hout
CALL steamtbl(FALSE,pexit,wt,wh,wrho,wvin,wxin,wbetain,wkapain,wcpin)
tout = wt - CKELVIN; rhoout = wrho
uout = rhouout / rhoout; uavg = cetac * (uout + cetacr * uout0)

IF(((uoutd-uout0) * (uout-uout0) .LT. -1.e-8*uavg * uavg) .AND. cetac .LT. 1.0_8) THEN
  uout = uavg
  rhouout = rhoout * uout
  rhohuout = rhouout * hout
ENDIF

CoolantTH%Tcool(iz) = T; CoolantTH%DenCool(iz) = rho
CoolantTH%hcool(iz) = h; CoolantTH%u(iz) = uout
CoolantTH%TCoolInOut(2, iz) = Tout; CoolantTH%TCoolInOut(1, iz+1) = Tout
CoolantTH%Rhou(iz) = rhouout; CoolantTH%rhohu(iz) = rhohuout

END SUBROUTINE
  
SUBROUTINE Tfcaltr(tfvol, Tfuel, Tfueld, Tcool, htCoef, Qf, Qfd, nr, ThVar, ThOpt, hgap)
USE PARAM
USE TYPEDEF,           ONLY : ThVar_Type,     ThOpt_Type
USE SteamTBL_mod,      ONLY : steamtbl
USE FuelProperty_Mod,  ONLY : FTHCON,       CTHCON,       FRHOCP,      CRHOCP,         ftfavg ,   RSGAPCOND
USE LU_Mod,            ONLY : DirLu1Dsolver
IMPLICIT NONE
REAL :: tfvol(nr-3)
REAL :: Tfuel(nr), TfuelD(nr), tcool, htcoef, qf, Qfd
INTEGER :: nr
TYPE(ThVar_Type) :: ThVar
TYPE(ThOPT_Type) :: ThOPT
REAL, OPTIONAL :: hgap

REAL :: diag(100), L(100), U(100), x(100), xd(100), b(100), xvol(100)
REAL :: kf(100), kfb(100), kfm(100), kfmb(100), rhocp(100)
REAL :: kmr, kml, kmrb, kmlb, kgap, kgapb, kgap2, kgap4, kgap4b, kconv, kconv1, kconv4, kconv4b, tworm, kfi
REAL, POINTER :: R(:)
REAL :: delR, delR2,  DelR2oDt, tw2odt, ri,  q, alpha, alphab
REAL :: rs, rw, tw, rgap
REAL :: err
REAL :: cetaf, cetafb
INTEGER :: npr, npr1, npr2, npr3, npr4, npr5
INTEGER :: i, m

cetaf = 0.5_8; cetafb = 1._8 - cetaf
npr = ThVar%npr; npr1 =ThVar%npr1
npr2 = ThVar%npr2; npr3 = Thvar%npr3
npr4= ThVar%npr4; npr5 = ThVar%npr5
rs = ThVar%rs; rw = ThVar%rw
tw = ThVar%tw; rgap = ThVar%rgap; 
R => ThVar%R
DelR = R(2) - R(1); DelR2 = DelR * DelR
DelR2oDt = DelR2 / ThVar%DelT
tw2odt= tw * tw /ThVar%DelT

hgap = ThOpt%hgap
kgap = cetaf * ThOpt%hgap * delr; 
kgapb = cetafb * ThOpt%hgap * delr; 
kgap2 = ThOpt%hgap * tw * rs / rgap
kgap4 = ThOpt%hgap * tw * (4._8 - tw / rgap) * rs / rgap
kgap4 = kgap4 * cetaf
kgap4b = kgap4 * cetafb / cetaf
tworm = tw / (rgap + 0.5_8 * tw) 
q = (cetaf*qf + cetafb*qfd)* DelR2



!Convert C to kelvin
DO i = 1, npr4
  xd(i) = tfueld(i) + CKELVIN                 !Celsius Unit
  x(i) = tfuel(i) + CKELVIN          !Kelvin Unit for temporary
ENDDO

DO i = 1, npr1
  !kf(i) = CondFuel(x(i))
  !kfb(i) = CondFuel(xd(i))
  kfi = FTHCON(x(i), 0.)
  kf(i) = cetaf * kfi
  !kfi = CondFuel(xd(i))
  kfb(i) = cetafb * kfi
  rhocp(i) = FRHOCP(x(i)) *DelR2oDt
  !PRINT*, FuelCp(x(i)) 
ENDDO

DO i = 1, npr
  kfm(i) = 0.5_8 * (kf(i) + kf(i+1))
  kfmb(i) = 0.5_8 * (kfb(i) + kfb(i+1))
ENDDO
!
DO i = npr2, npr4
  kfi = CTHCON(x(i)) 
  kf(i) = cetaf * kfi
  !kfi = CondClad(xd(i)) 
  kfb(i) = cetafb * kfi
  rhocp(i) = CRHOCP(x(i)) * tw2odt
  !  kf(i) = CondClad(x(i))
ENDDO

kmr = 0.5_8 * (kf(npr3) + kf(npr4))
kml = 0.5_8 * (kf(npr3) + kf(npr2))

kmrb = 0.5_8 * (kfb(npr3) + kfb(npr4))
kmlb = 0.5_8 * (kfb(npr3) + kfb(npr2))

IF(PRESENT(hGap)) THEN
  IF(ThOpt%hGapModel .EQ. 3) THEN
    hgap = RSGAPCOND(kf(npr1), kf(npr2), x(npr1), x(npr2))
    kgap = cetaf * hgap * delr;
    kgapb = cetafb * hgap * delr;
    kgap2 = hgap * tw * rs / rgap
    kgap4 = hgap * tw * (4._8 - tw / rgap) * rs / rgap
    kgap4 = kgap4 * cetaf
    kgap4b = kgap4 * cetafb / cetaf
  END IF
END IF

xd(1:npr4) = tfueld(1:npr4)
x(1:npr4) = tfuel(1:npr4)

m=1;


!Center-line element
Diag(m)=rhocp(m)+4*kf(m)
U(m)=-4*kf(m)
b(m)=q+xd(m)*rhocp(m)-4*kfb(m)*xd(m)+4*kfb(m)*xd(m+1)

i=m  ! i = m - 1 point

DO m = 2, npr
  ri=1/dble(i)
  Diag(m) = rhocp(m) + (kfm(i)+kfm(m)+0.5*(kfm(m)-kfm(i))*ri)
  L(m) = -kfm(i)*(1-0.5*ri)
  U(m) = -kfm(m)*(1+0.5*ri)
  b(m) = q+xd(m)*rhocp(m)-(kfmb(i)+kfmb(m)+0.5*(kfmb(m)-kfmb(i))*ri)*xd(m)+kfmb(i)*(1-0.5*ri)*xd(m-1)+kfmb(m)*(1+0.5*ri)*xd(m+1)
  i = m
ENDDO

m = npr1
alpha = kgap*(1-kf(m-1)/kf(m))
alphab = alpha/cetaf*cetafb
Diag(m) = rhocp(m)+2._8*(kf(m)+kgap*(1+0.5/npr))+alpha
L(m) = -2._8*kf(m)
U(m) = -2._8*kgap*(1+0.5/npr)-alpha
b(m) = q + xd(m)*rhocp(m)-(2*(kfb(m)+kgapb*(1+0.5/npr))+alphab)*xd(m) +2*kfb(m)*xd(m-1)+(2*kgapb*(1+0.5/npr)+alphab)*xd(m+1)

m = npr2
alpha = 2._8*kgap2*(kf(m+1)/kf(m)-1)
alphab = alpha/cetaf*cetafb
Diag(m) = rhocp(m)+8._8*kf(m)+kgap4-alpha
L(m) = -kgap4+alpha
U(m) = -8._8*kf(m)
b(m)=xd(m)*rhocp(m)-(8._8*kfb(m)+kgap4b)*xd(m)+kgap4b*xd(m-1) +8._8*kfb(m)*xd(m+1)+xd(m)*alphab-xd(m-1)*alphab

m=npr3
Diag(m) = rhocp(m)+4._8*(kmr+kml)+tworm*(kmr-kml)
L(m) = -kml*(4._8-tworm)
U(m) = -kmr*(4._8+tworm)
b(m) = xd(m)*rhocp(m)-(4._8*(kmrb+kmlb)+tworm*(kmrb-kmlb))*xd(m)+kmlb*(4._8-tworm)*xd(m-1)+kmrb*(4._8+tworm)*xd(m+1)

m = npr4
kconv1 = cetaf*htcoef*tw
alpha = 2._8*kconv1*(1-kf(m-1)/kf(m))
alphab = alpha/cetaf*cetafb
kconv = htcoef*tw*(4+tw/rw)
kconv4 = cetaf*kconv
kconv4b = cetafb*kconv

Diag(m) = rhocp(m)+8*kf(m)+kconv4+alpha
L(m) = -8*kf(m)
b(m) = xd(m)*rhocp(m)-(8*kfb(m)+kconv4b+alphab)*xd(m)+8*kfb(m)*xd(m-1)+(kconv+alpha/cetaf)*Tcool

CALL DirLu1Dsolver(Diag(1:npr4), L(1:npr4), U(1:npr4), x(1:npr4), b(1:npr4), npr4)

DO m = 1, npr4
  tfuel(m) = x(m)  ! Celcious Unit
ENDDO

CALL ftfavg(x, xvol(1:npr5), r, npr)
tfvol(1:npr2) = xvol(1:npr2) + CKELVIN   !Kelvin Unit
tfuel(npr5) = xvol(npr5)

!NULLIFY(R)
!DO m = 1, 12
!PRINT'(5es14.6)', x(m), xd(m), q, rhocp(m), x(m)-xd(m)
!END DO 

END SUBROUTINE

SUBROUTINE Tfcaltr_pwshape(tfvol, Tfuel, Tfueld, Tcool, htCoef, Qf, Qfd, pwshape, pwshaped, nr, ThVar, ThOpt, hgap)
  USE PARAM
  USE TYPEDEF,           ONLY : ThVar_Type,     ThOpt_Type
  USE SteamTBL_mod,      ONLY : steamtbl
  USE FuelProperty_Mod,  ONLY : FTHCON,       CTHCON,       FRHOCP,      CRHOCP,         ftfavg ,   RSGAPCOND
  USE LU_Mod,            ONLY : DirLu1Dsolver
  IMPLICIT NONE
  REAL :: tfvol(nr-3)
  REAL :: Tfuel(nr), TfuelD(nr), tcool, htcoef, qf, Qfd, pwshape(nr), pwshaped(nr)
  INTEGER :: nr
  TYPE(ThVar_Type) :: ThVar
  TYPE(ThOPT_Type) :: ThOPT
  REAL, OPTIONAL :: hgap
  
  REAL :: diag(100), L(100), U(100), x(100), xd(100), b(100), xvol(100)
  REAL :: kf(100), kfb(100), kfm(100), kfmb(100), rhocp(100)
  REAL :: kmr, kml, kmrb, kmlb, kgap, kgapb, kgap2, kgap4, kgap4b, kconv, kconv1, kconv4, kconv4b, tworm, kfi
  REAL, POINTER :: R(:)
  REAL :: delR, delR2,  DelR2oDt, tw2odt, ri,  qdist(nr), alpha, alphab
  REAL :: rs, rw, tw, rgap
  REAL :: err
  REAL :: cetaf, cetafb
  INTEGER :: npr, npr1, npr2, npr3, npr4, npr5
  INTEGER :: i, m
  
  cetaf = 0.5_8; cetafb = 1._8 - cetaf
  npr = ThVar%npr; npr1 =ThVar%npr1
  npr2 = ThVar%npr2; npr3 = Thvar%npr3
  npr4= ThVar%npr4; npr5 = ThVar%npr5
  rs = ThVar%rs; rw = ThVar%rw
  tw = ThVar%tw; rgap = ThVar%rgap; 
  R => ThVar%R
  DelR = R(2) - R(1); DelR2 = DelR * DelR
  DelR2oDt = DelR2 / ThVar%DelT
  tw2odt= tw * tw /ThVar%DelT
  
  hgap = ThOpt%hgap
  kgap = cetaf * ThOpt%hgap * delr; 
  kgapb = cetafb * ThOpt%hgap * delr; 
  kgap2 = ThOpt%hgap * tw * rs / rgap
  kgap4 = ThOpt%hgap * tw * (4._8 - tw / rgap) * rs / rgap
  kgap4 = kgap4 * cetaf
  kgap4b = kgap4 * cetafb / cetaf
  tworm = tw / (rgap + 0.5_8 * tw) 

  qdist = 0.
  DO i = 1, nr
    qdist(i) = (cetaf*qf*pwshape(i) + cetafb*qfd*pwshaped(i))* DelR2
  END DO
  
  
  
  !Convert C to kelvin
  DO i = 1, npr4
    xd(i) = tfueld(i) + CKELVIN                 !Celsius Unit
    x(i) = tfuel(i) + CKELVIN          !Kelvin Unit for temporary
  ENDDO
  
  DO i = 1, npr1
    !kf(i) = CondFuel(x(i))
    !kfb(i) = CondFuel(xd(i))
    kfi = FTHCON(x(i), 0.)
    kf(i) = cetaf * kfi
    !kfi = CondFuel(xd(i))
    kfb(i) = cetafb * kfi
    rhocp(i) = FRHOCP(x(i)) *DelR2oDt
    !PRINT*, FuelCp(x(i)) 
  ENDDO
  
  DO i = 1, npr
    kfm(i) = 0.5_8 * (kf(i) + kf(i+1))
    kfmb(i) = 0.5_8 * (kfb(i) + kfb(i+1))
  ENDDO
  !
  DO i = npr2, npr4
    kfi = CTHCON(x(i)) 
    kf(i) = cetaf * kfi
    !kfi = CondClad(xd(i)) 
    kfb(i) = cetafb * kfi
    rhocp(i) = CRHOCP(x(i)) * tw2odt
    !  kf(i) = CondClad(x(i))
  ENDDO
  
  kmr = 0.5_8 * (kf(npr3) + kf(npr4))
  kml = 0.5_8 * (kf(npr3) + kf(npr2))
  
  kmrb = 0.5_8 * (kfb(npr3) + kfb(npr4))
  kmlb = 0.5_8 * (kfb(npr3) + kfb(npr2))
  
  IF(PRESENT(hGap)) THEN
    IF(ThOpt%hGapModel .EQ. 3) THEN
      hgap = RSGAPCOND(kf(npr1), kf(npr2), x(npr1), x(npr2))
      kgap = cetaf * hgap * delr;
      kgapb = cetafb * hgap * delr;
      kgap2 = hgap * tw * rs / rgap
      kgap4 = hgap * tw * (4._8 - tw / rgap) * rs / rgap
      kgap4 = kgap4 * cetaf
      kgap4b = kgap4 * cetafb / cetaf
    END IF
  END IF
  
  xd(1:npr4) = tfueld(1:npr4)
  x(1:npr4) = tfuel(1:npr4)
  
  m=1;
  
  
  !Center-line element
  Diag(m)=rhocp(m)+4*kf(m)
  U(m)=-4*kf(m)
  b(m)=qdist(m)+xd(m)*rhocp(m)-4*kfb(m)*xd(m)+4*kfb(m)*xd(m+1)
  
  i=m  ! i = m - 1 point
  
  DO m = 2, npr
    ri=1/dble(i)
    Diag(m) = rhocp(m) + (kfm(i)+kfm(m)+0.5*(kfm(m)-kfm(i))*ri)
    L(m) = -kfm(i)*(1-0.5*ri)
    U(m) = -kfm(m)*(1+0.5*ri)
    b(m) = qdist(m)+xd(m)*rhocp(m)-(kfmb(i)+kfmb(m)+0.5*(kfmb(m)-kfmb(i))*ri)*xd(m)+kfmb(i)*(1-0.5*ri)*xd(m-1)+kfmb(m)*(1+0.5*ri)*xd(m+1)
    i = m
  ENDDO
  
  m = npr1
  alpha = kgap*(1-kf(m-1)/kf(m))
  alphab = alpha/cetaf*cetafb
  Diag(m) = rhocp(m)+2._8*(kf(m)+kgap*(1+0.5/npr))+alpha
  L(m) = -2._8*kf(m)
  U(m) = -2._8*kgap*(1+0.5/npr)-alpha
  b(m) = qdist(m-1) + xd(m)*rhocp(m)-(2*(kfb(m)+kgapb*(1+0.5/npr))+alphab)*xd(m) +2*kfb(m)*xd(m-1)+(2*kgapb*(1+0.5/npr)+alphab)*xd(m+1)
  
  m = npr2
  alpha = 2._8*kgap2*(kf(m+1)/kf(m)-1)
  alphab = alpha/cetaf*cetafb
  Diag(m) = rhocp(m)+8._8*kf(m)+kgap4-alpha
  L(m) = -kgap4+alpha
  U(m) = -8._8*kf(m)
  b(m)=xd(m)*rhocp(m)-(8._8*kfb(m)+kgap4b)*xd(m)+kgap4b*xd(m-1) +8._8*kfb(m)*xd(m+1)+xd(m)*alphab-xd(m-1)*alphab
  
  m=npr3
  Diag(m) = rhocp(m)+4._8*(kmr+kml)+tworm*(kmr-kml)
  L(m) = -kml*(4._8-tworm)
  U(m) = -kmr*(4._8+tworm)
  b(m) = xd(m)*rhocp(m)-(4._8*(kmrb+kmlb)+tworm*(kmrb-kmlb))*xd(m)+kmlb*(4._8-tworm)*xd(m-1)+kmrb*(4._8+tworm)*xd(m+1)
  
  m = npr4
  kconv1 = cetaf*htcoef*tw
  alpha = 2._8*kconv1*(1-kf(m-1)/kf(m))
  alphab = alpha/cetaf*cetafb
  kconv = htcoef*tw*(4+tw/rw)
  kconv4 = cetaf*kconv
  kconv4b = cetafb*kconv
  
  Diag(m) = rhocp(m)+8*kf(m)+kconv4+alpha
  L(m) = -8*kf(m)
  b(m) = xd(m)*rhocp(m)-(8*kfb(m)+kconv4b+alphab)*xd(m)+8*kfb(m)*xd(m-1)+(kconv+alpha/cetaf)*Tcool
  
  CALL DirLu1Dsolver(Diag(1:npr4), L(1:npr4), U(1:npr4), x(1:npr4), b(1:npr4), npr4)
  
  DO m = 1, npr4
    tfuel(m) = x(m)  ! Celcious Unit
  ENDDO
  
  CALL ftfavg(x, xvol(1:npr5), r, npr)
  tfvol(1:npr2) = xvol(1:npr2) + CKELVIN   !Kelvin Unit
  tfuel(npr5) = xvol(npr5)
  
  !NULLIFY(R)
  !DO m = 1, 12
  !PRINT'(5es14.6)', x(m), xd(m), q, rhocp(m), x(m)-xd(m)
  !END DO 
  
  END SUBROUTINE



SUBROUTINE Tfcaltr_Qshape(tfvol, Tfuel, Tfueld, Tcool, htCoef, Qf, Qfd, Qshape, Qshaped, nr, ThVar, ThOpt)
USE PARAM
USE TYPEDEF,           ONLY : ThVar_Type,     ThOpt_Type
USE SteamTBL_mod,      ONLY : steamtbl
USE FuelProperty_Mod,  ONLY : CondFuel,       CondClad,       FuelCp,      CladCp,         ftfavg 
USE LU_Mod,            ONLY : DirLu1Dsolver
IMPLICIT NONE
REAL :: tfvol(nr-3)
REAL :: Tfuel(nr), TfuelD(nr), tcool, htcoef, qf, Qfd
REAL :: Qshaped(nr), Qshape(nr)
INTEGER :: nr
TYPE(ThVar_Type) :: ThVar
TYPE(ThOPT_Type) :: ThOPT

REAL :: diag(100), L(100), U(100), x(100), xd(100), b(100), xvol(100), Qdist(100)
REAL :: kf(100), kfb(100), kfm(100), kfmb(100), rhocp(100)
REAL :: kmr, kml, kmrb, kmlb, kgap, kgapb, kgap2, kgap4, kgap4b, kconv, kconv1, kconv4, kconv4b, tworm, kfi
REAL, POINTER :: R(:)
REAL :: delR, delR2,  DelR2oDt, tw2odt, ri,  alpha, alphab
REAL :: rs, rw, tw, rgap
REAL :: err
REAL :: cetaf, cetafb
INTEGER :: npr, npr1, npr2, npr3, npr4, npr5
INTEGER :: i, m

cetaf = 0.5_8; cetafb = 1._8 - cetaf
npr = ThVar%npr; npr1 =ThVar%npr1
npr2 = ThVar%npr2; npr3 = Thvar%npr3
npr4= ThVar%npr4; npr5 = ThVar%npr5
rs = ThVar%rs; rw = ThVar%rw
tw = ThVar%tw; rgap = ThVar%rgap; 
R => ThVar%R
DelR = R(2) - R(1); DelR2 = DelR * DelR
DelR2oDt = DelR2 / ThVar%DelT
tw2odt= tw * tw /ThVar%DelT
kgap = cetaf * ThOpt%hgap * delr; 
kgapb = cetafb * ThOpt%hgap * delr; 

kgap2 = ThOpt%hgap * tw * rs / rgap
kgap4 = ThOpt%hgap * tw * (4._8 - tw / rgap) * rs / rgap
kgap4 = kgap4 * cetaf
kgap4b = kgap4 * cetafb / cetaf
tworm = tw / (rgap + 0.5_8 * tw) 
!q = (cetaf*qf + cetafb*qfd)* DelR2
DO i = 1, nr
  Qdist(i) = (cetaf * qf * Qshape(i) + cetafb * qfd *Qshaped(i)) * DelR2
ENDDO

!Convert C to kelvin
DO i = 1, npr4
  xd(i) = tfueld(i) + CKELVIN                 !Celsius Unit
  x(i) = tfuel(i) + CKELVIN          !Kelvin Unit for temporary
ENDDO

DO i = 1, npr1
  !kf(i) = CondFuel(x(i))
  !kfb(i) = CondFuel(xd(i))
  kfi = CondFuel(x(i))
  kf(i) = cetaf * kfi
  !kfi = CondFuel(xd(i))
  kfb(i) = cetafb * kfi
  rhocp(i) = FuelCp(x(i)) *DelR2oDt
ENDDO

DO i = 1, npr
  kfm(i) = 0.5_8 * (kf(i) + kf(i+1))
  kfmb(i) = 0.5_8 * (kfb(i) + kfb(i+1))
ENDDO
!
DO i = npr2, npr4
  kfi = CondClad(x(i)) 
  kf(i) = cetaf * kfi
  !kfi = CondClad(xd(i)) 
  kfb(i) = cetafb * kfi
  rhocp(i) = CladCp(x(i)) * tw2odt
  !  kf(i) = CondClad(x(i))
ENDDO

kmr = 0.5_8 * (kf(npr3) + kf(npr4))
kml = 0.5_8 * (kf(npr3) + kf(npr2))

kmrb = 0.5_8 * (kfb(npr3) + kfb(npr4))
kmlb = 0.5_8 * (kfb(npr3) + kfb(npr2))

xd(1:npr4) = tfueld(1:npr4)
x(1:npr4) = tfuel(1:npr4)

m=1;

!Center-line element
Diag(m)=rhocp(m)+4*kf(m)
U(m)=-4*kf(m)
b(m)=qdist(m)+xd(m)*rhocp(m)-4*kfb(m)*xd(m)+4*kfb(m)*xd(m+1)

i=m  ! i = m - 1 point

DO m = 2, npr
  ri=1/dble(i)
  Diag(m) = rhocp(m) + (kfm(i)+kfm(m)+0.5*(kfm(m)-kfm(i))*ri)
  L(m) = -kfm(i)*(1-0.5*ri)
  U(m) = -kfm(m)*(1+0.5*ri)
  b(m) = qdist(m)+xd(m)*rhocp(m)-(kfmb(i)+kfmb(m)+0.5*(kfmb(m)-kfmb(i))*ri)*xd(m)+kfmb(i)*(1-0.5*ri)*xd(m-1)+kfmb(m)*(1+0.5*ri)*xd(m+1)
  i = m
ENDDO

m = npr1
alpha = kgap*(1-kf(m-1)/kf(m))
alphab = alpha/cetaf*cetafb
Diag(m) = rhocp(m)+2._8*(kf(m)+kgap*(1+0.5/npr))+alpha
L(m) = -2._8*kf(m)
U(m) = -2._8*kgap*(1+0.5/npr)-alpha
b(m) = qdist(m) + xd(m)*rhocp(m)-(2*(kfb(m)+kgapb*(1+0.5/npr))+alphab)*xd(m) +2*kfb(m)*xd(m-1)+(2*kgapb*(1+0.5/npr)+alphab)*xd(m+1)

m = npr2
alpha = 2._8*kgap2*(kf(m+1)/kf(m)-1)
alphab = alpha/cetaf*cetafb
Diag(m) = rhocp(m)+8._8*kf(m)+kgap4-alpha
L(m) = -kgap4+alpha
U(m) = -8._8*kf(m)
b(m)=xd(m)*rhocp(m)-(8._8*kfb(m)+kgap4b)*xd(m)+kgap4b*xd(m-1) +8._8*kfb(m)*xd(m+1)+xd(m)*alphab-xd(m-1)*alphab

m=npr3
Diag(m) = rhocp(m)+4._8*(kmr+kml)+tworm*(kmr-kml)
L(m) = -kml*(4._8-tworm)
U(m) = -kmr*(4._8+tworm)
b(m) = xd(m)*rhocp(m)-(4._8*(kmrb+kmlb)+tworm*(kmrb-kmlb))*xd(m)+kmlb*(4._8-tworm)*xd(m-1)+kmrb*(4._8+tworm)*xd(m+1)

m = npr4
kconv1 = cetaf*htcoef*tw
alpha = 2._8*kconv1*(1-kf(m-1)/kf(m))
alphab = alpha/cetaf*cetafb
kconv = htcoef*tw*(4+tw/rw)
kconv4 = cetaf*kconv
kconv4b = cetafb*kconv

Diag(m) = rhocp(m)+8*kf(m)+kconv4+alpha
L(m) = -8*kf(m)
b(m) = xd(m)*rhocp(m)-(8*kfb(m)+kconv4b+alphab)*xd(m)+8*kfb(m)*xd(m-1)+(kconv+alpha/cetaf)*Tcool

CALL DirLu1Dsolver(Diag(1:npr4), L(1:npr4), U(1:npr4), x(1:npr4), b(1:npr4), npr4)

DO m = 1, npr4
  tfuel(m) = x(m)  ! Celcious Unit
ENDDO

CALL ftfavg(x, xvol(1:npr5), r, npr)
tfvol(1:npr2) = xvol(1:npr2) + CKELVIN   !Kelvin Unit
tfuel(npr5) = xvol(npr5)

NULLIFY(R)

END SUBROUTINE


SUBROUTINE InitTransientTH(THInfo, ThVar, nxy, nz)
USE PARAM
USE TYPEDEF,             ONLY : ThInfo_Type,         ThVar_Type,       CoolantTH_Type,        &
                                FuelTh_Type
USE ALLOCS
IMPLICIT NONE
TYPE(ThInfo_Type) :: ThInfo
TYPE(ThVar_Type) :: ThVar
INTEGER :: nxy, nz

TYPE(CoolantTH_Type), POINTER :: CoolantTH(:)
TYPE(FuelTH_Type), POINTER :: FuelTH(:)

INTEGER :: ixy, iz

CoolantTH => ThInfo%CoolantTH
FuelTH => ThInfo%FuelTH
DO ixy = 1, nxy
  IF(.NOT. CoolantTH(ixy)%lfuel) CYCLE
  CoolantTH(ixy)%TCoold(1:nz) = CoolantTH(ixy)%TCool(1:nz)
  CoolantTH(ixy)%DenCoold(1:nz) = CoolantTH(ixy)%DenCool(1:nz)
  CoolantTH(ixy)%hCoold(1:nz) = CoolantTH(ixy)%hCool(1:nz)
  CoolantTH(ixy)%rhoud(0:nz) = CoolantTH(ixy)%rhou(0:nz)
  CoolantTH(ixy)%rhohud(0:nz) = CoolantTH(ixy)%rhohu(0:nz)
  CoolantTH(ixy)%qeffd(1:nz) = CoolantTH(ixy)%qeff(1:nz)
  FuelTH(ixy)%qvold(1:nz) = FuelTH(ixy)%qvol(1:nz)
  FuelTH(ixy)%tfueld(:, 1:nz) = FuelTH(ixy)%tfuel(:, 1:nz)
ENDDO

END SUBROUTINE

SUBROUTINE PrintHeatPower(Core, CmInfo, FmInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,                    ONLY : CoreInfo_Type,              FmInfo_Type,               CmInfo_Type,          &
                                       ThInfo_Type,                GroupInfo_Type,            TranCntl_Type,        &
                                       PE_TYPE,                                                                     &
                                       FuelTH_Type,                CoolantTH_Type,            Pin_Type,             &
                                       FxrInfo_Type,               PinXS_Type
USE CNTL,                       ONLY : nTracerCntl_Type
USE TH_Mod,                     ONLY : ThOpt,                      ThVar,                                           &
                                       CalRelPower,                SetPwShape,                tfcaltr, Tfcaltr_pwshape
USE FuelProperty_Mod,           ONLY : fhtcoef
USE timer,                      ONLY : nTracer_dclock,             TimeChk
USE FILES,                      ONLY : IO8
USE IOUTIL,                     ONLY : message
USE TIMER,                      ONLY : nTracer_dclock,             TimeChk
#ifdef MPI_ENV
USE MPIComm_Mod,                ONLY : MPI_SYNC
#endif
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(CmInfo_Type) :: CmInfo
TYPE(FmInfo_Type) :: FmInfo
TYPE(ThInfo_Type) :: ThInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Fxrinfo_type), POINTER :: Fxr(:,:)
TYPE(PinXS_Type), POINTER :: PinXS(:, :)
TYPE(FuelTh_Type), POINTER :: FuelTH(:)      
TYPE(CoolantTH_Type), POINTER :: CoolantTH(:)
REAL, POINTER :: RelPower(:, :)

REAL, POINTER :: PwShape(:, :), PwShaped(:, :)

INTEGER :: nowstep
INTEGER :: ng, nxy, nzth, nr, nsurf, nz, nchannel
INTEGER :: Comm
INTEGER :: ixy, iz, i
INTEGER :: ixya, ichtyp


REAL :: Tbeg, Tend
REAL :: dz, qprime, qf, qfd, qflux, qeff, qeffd
REAL :: Tdopld, TdoplMax, TdoplMax1
REAL :: TBulk, Tsurf, Tout, ToutAvg, Tfmax, htcoeff
REAL :: Deq, Fracdf, afp
REAL :: PEXIT, BTEMP

LOGICAL :: Master, Slave


IF(nTracerCntl%lBenchXS) RETURN

NowStep = TranCntl%NowStep
ThVar%DelT = TranCntl%DelT(NowStep)
Pin => Core%Pin;                  Fxr => FmInfo%FXR
PinXs => CmInfo%PinXS;            RelPower => ThInfo%RelPower
FuelTH => ThInfo%FuelTH;          CoolantTH => ThInfo%CoolantTH

Tbeg = nTracer_dclock(false, false)

ALLOCATE(PwShape(ThVar%npr5, ThVar%nzth), PwShaped(ThVar%npr5, ThVar%nzth))


Master = PE%CmfdMaster; Slave = .NOT. Master
Comm = PE%MPI_CMFD_COMM
WRITE(MESG, '(A)') 'Performing T/H Calculation...'
IF(Master) CALL MESSAGE(io8,TRUE,TRUE,'Performing T/H Calculation...')

PEXIT = ThInfo%PExit; Btemp = ThVar%BoilingTemp  
nxy = Core%nxy; nzth = ThVar%nzth; ng = GroupInfo%ng 

nr = ThVar%npr5;nsurf = Thvar%npr4; nz = ThVar%nzth

fracdf = THVar%fracdf;DEQ =ThVar%Deq; afp = ThVar%afp

CALL CalRelPower(Core, CmInfo, RelPower, GroupInfo%ng, nTracerCntl, PE, .TRUE.)

nchannel = 0; ToutAvg = 0; Tfmax = 0
DO ixy = 1, nxy
  IF(.NOT. CoolantTH(ixy)%lfuel) CYCLE
  Deq = ThVar%Deq
  IF(nTracerCntl%lThChConf) THEN 
    ixya = Pin(ixy)%iasy
    ichtyp = Core%ThChMap(ixya)
    IF(ichtyp .NE. 0) THEN
      Deq = ThVar%ThCh(ichtyp)%Deq
    END IF
  END IF
  nchannel = nchannel + 1
  !Set Power Shape
  CALL SetPwShape(Pwshape, Core, FmInfo%Fxr, FmInfo%Power, ixy, nzth, ThVar%npr5, Thvar, ThOpt, PE)
  !WRITE(88, '(100F10.5)') (PwShape(i, 1), i = 1, nr)
  WRITE(88, '(100F10.5)') (FuelTh(ixy)%tfuel(i, 1), i = 1, nr)
  WRITE(87, '(100F10.5)') CoolantTH(ixy)%TCool(1)
  !WRITE(88, '(100F10.5)')
ENDDO
DEALLOCATE(PwShape, PwShaped)
    
#ifdef MPI_ENV
CALL MPI_SYNC(Comm)
#endif

END SUBROUTINE
