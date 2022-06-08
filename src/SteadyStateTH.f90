#include <defines.h>
SUBROUTINE SteadyStateTH(Core, CmInfo, FmInfo, ThInfo, Eigv, ng, GroupInfo, nTRACERCntl, PE)

USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type, CMInfo_Type, FmInfo_Type, ThInfo_Type, GroupInfo_Type, PE_Type, FxrInfo_Type, PinXS_Type, Pin_Type, ThVar_Type, THOpt_Type, FuelTh_Type, CoolantTH_Type
USE TH_Mod,         ONLY : ThOpt, ThVar, CalRelPower, SteadyCoolantTH, SteadyFuelConduction, SubCellPowerProfile, SimpleTH, SetPwShape, SetBurnupDist, Grp_RelPower, SetUSErDefTH, recalculate_vars!, &
                           !hGapArray, SteadyCoolantTH_ThCh, SteadyFuelConduction_ThCh
USE SubChannel_mod, ONLY : SubChannelTH
USE MATRATH_mod,    ONLY : MATRA_TH
USE CTFTH_mod,      ONLY : CTF_TH
USE ESCOTTH_mod,    ONLY : ESCOT_TH
USE timer,          ONLY : nTracer_dclock, TimeChk
USE BasicOperation, ONLY : CP_CA, CP_VA
USE CNTL,           ONLY : nTracerCntl_Type
USE itrcntl_mod,    ONLY : ItrCntl
USE FILES,          ONLY : IO8
USE IOUTIL,         ONLY : message, Terminate
USE TIMER,          ONLY : nTracer_dclock, TimeChk

USE SubChCoupling_mod,              ONLY: is_coupled, CodeName
USE Anderson_Acceleration_SIMPLETH, ONLY: init_AA, Anderson_acc, m_AA

#ifdef MPI_ENV
USE MPIComm_Mod, ONLY : MPI_SYNC
#endif

IMPLICIT NONE

TYPE (CoreInfo_Type) :: Core
TYPE (CMInfo_Type) :: CmInfo
TYPE (FMInfo_Type) :: FmInfo
TYPE (ThInfo_Type) :: ThInfo
TYPE (GroupInfo_Type) :: GroupInfo
TYPE (nTracerCntl_Type) :: nTRACERCntl
TYPE (PE_Type) :: PE
REAL :: Eigv
INTEGER :: ng
! ----------------------------------------------------
TYPE (Pin_Type), POINTER :: Pin(:)
TYPE (Fxrinfo_type), POINTER :: Fxr(:,:)
TYPE (PinXS_Type), POINTER :: PinXS(:, :)
TYPE (FuelTh_Type), POINTER :: FuelTH(:)
TYPE (CoolantTH_Type), POINTER :: CoolantTH(:)

REAL, POINTER, DIMENSION(:)   :: RhoU
REAL, POINTER, DIMENSION(:,:) :: Profile, PwShape, BurnupDist, RelPower!, GrpRelPower 
INTEGER :: xyb, xye, nxy, nzth, nchannel, npr5   !--- CNJ Edit : Domain Decomposition + MPI
INTEGER :: ixy, iz, iter, n, Comm, ncall
REAL :: PEXIT, Tout, toutavg, tfmax1, tfmax, tdopd, tdoplmax, TdoplMax0, tbeg, tend, powlin, powlv, Btemp, tfavg_max, tmod_max, hgapavg, nhgap
LOGICAL :: Master, Slave, lsimpleTH, lMATRA
LOGICAL :: lftemp_ex = FALSE
LOGICAL :: lfirst    = TRUE
LOGICAL :: is_AA     = TRUE
LOGICAL :: is_reset  = FALSE
! ----------------------------------------------------

ncall     = 0
tfavg_max = 0.
tmod_max  = 0.

IF (nTracerCntl%libtyp.NE.11 .AND. nTRACERCntl%lBenchXs) RETURN

Tbeg      = nTracer_dclock(false, false)
lSimpleTh = nTracerCntl%lSimpleTH
Master    = PE%CmfdMaster
Slave     = .NOT. Master
Comm      = PE%MPI_CMFD_COMM
is_AA     = ThOpt%AA_STH

WRITE (MESG, '(A)') 'Performing T/H Calculation...'
IF (Master) CALL MESSAGE(io8,TRUE,TRUE,MESG)
IF (.NOT. is_coupled) THEN
  WRITE (MESG, '(A,i3,A)') 'T/H Calculation with Anderson Acceleration (Depth  =',m_AA,')'
  IF (lfirst .and. is_AA .and. MASTER) CALL MESSAGE(io8,TRUE,TRUE,MESG)
END IF

Pin       => Core%Pin
Fxr       => FmInfo%Fxr
PinXS     => CMInfo%PinXS
FuelTH    => ThInfo%FuelTH
CoolantTH => ThInfo%CoolantTH
RelPower  => ThInfo%RelPower

powlin = ThInfo%PowLin
PowLv  = ThInfo%PowLv
PEXIT  = ThInfo%PEXIT
nxy    = Core%nxy
nzth   = ThVar%nzth
npr5   = ThVar%npr5
Btemp  = ThVar%BoilingTemp

!--- CNJ Edit : Domain Decomposition + MPI
xyb = PE%myPinBeg
xye = PE%myPinEnd
IF (PE%RTMASTER) THEN
  xyb = 1; xye = nxy
END IF

IF (lSimpleTH) ALLOCATE (Profile (2*ThVar%npr, nzth))
ALLOCATE (Pwshape      (ThVar%npr5, nzth))
ALLOCATE (BurnUpDist (0:ThVar%npr5, nzth))
! ----------------------------------------------------
IF (ItrCntl%cmfdit .NE. 0) THEN
  ! CAL : Rel. Pw.
  CALL CalRelPower(Core, CmInfo, RelPower, ng, nTracerCntl, PE, TRUE)
    
  ! CAL : Grp. Rel. Pw.
  IF (nTracerCntl%ThCh_mod .GE. 1) THEN 
    !IF (nTRACERCntl%lthch_tf) THEN 
      !CALL Grp_RelPower(Core, CmInfo, RelPower,    RelPower, ng, nTracerCntl, PE, TRUE)
    CALL Grp_RelPower(Core, CmInfo, RelPower, ng, nTracerCntl, PE, TRUE)
    !ELSE
    !  ALLOCATE(GrpRelPower(0:nzth,0:nxy))
    !  CALL Grp_RelPower(Core, CmInfo, RelPower, GrpRelPower, ng, nTracerCntl, PE, TRUE)
    !END IF
  END IF
END IF
! ----------------------------------------------------
DO iter = 1, 1
  IF (lSimpleTH) THEN
    nchannel = 0
    toutavg  = 0
    tfmax    = 0
    tdoplmax = 0
    
    DO ixy = xyb, xye
      IF (.NOT. CoolantTH(ixy)%lfuel) CYCLE
      
      nchannel = nchannel + 1
      
      CALL SubCellPowerProfile(Core, FmInfo, Profile, ixy, THVar%npr, nzth, ng, GroupInfo, PE)
      CALL SimpleTH(PowLin, PowLv, Tout, Tfmax1, RelPower(1:nzth, ixy), Profile, FuelTh(ixy), CoolantTH(ixy), ThVar, ThOpt, PE)
      
      tfmax   = max(tfmax, tfmax1)
      ToutAvg = toutavg + Tout
    END DO
    
    CALL SetGTnGapCoolTemp(Core, ThInfo, PE)
    ToutAvg = ToutAvg / nchannel
	ELSE
    nchannel = 0
    toutavg = 0
    ! ------------------------------------------------
		IF (nTracerCntl%lMATRA) THEN
			CALL SubChannelTH(Core, FmInfo, CmInfo, ThInfo, GroupInfo, nTracerCntl, PE)
      
      !The improved subroutine for coupling nTRACER/MATRA is under development
			DO ixy = xyb, xye
				IF (.NOT. CoolantTH(ixy)%lfuel) CYCLE
        
				nchannel = nchannel + 1
				ToutAvg  = toutavg + ThInfo%Tcool(nzth, ixy)
      END DO
      
		ELSE IF (CodeName=='CTF') THEN
      tfmax = 0._8
      CALL CTF_TH(Core, ThInfo, nTracerCntl,ToutAvg, PE)
      nchannel = 1
    ELSE IF (CodeName=='ESCOT') THEN
      tfmax = 0._8
      CALL ESCOT_TH(Core, ThInfo, nTracerCntl, ToutAvg, lftemp_ex, PE, ng)
      nchannel = 1
    !elseif (CodeName=='WRAPPER' .and. is_coupled) then
    !  tfmax = 0._8
    !  CALL Wrapper_TH(Core, ThInfo, nTracerCntl, ToutAvg, PE, ng)
    !  lftemp_ex = .true.
    !  nchannel=1
    ELSE
      DO ixy = xyb, xye
        IF (.NOT. CoolantTH(ixy)%lfuel) CYCLE
        
        nchannel = nchannel + 1
        
        !IF (nTracerCntl%lthchConf) THEN
        !  IF (nTRACERCntl%ThCh_mod .GE. 1 .AND. .NOT. nTRACERCntl%lthch_tf) THEN 
        !    CALL SteadyCoolantTH_ThCh(Core, PowLin, PowLv, PEXIT, Tout, GrpRelPower(1:nzth, ixy), CoolantTH(ixy), ThVar, ThOpt, PE, ixy)
        !  ELSE 
        !    CALL SteadyCoolantTH_ThCh(Core, PowLin, PowLv, PEXIT, Tout,    RelPower(1:nzth, ixy), CoolantTH(ixy), ThVar, ThOpt, PE, ixy)
        !  END IF
        !ELSE
          !IF (nTRACERCntl%ThCh_mod .GE. 1 .AND. .NOT. nTRACERCntl%lthch_tf) THEN 
          !  CALL SteadyCoolantTH(PowLin, PowLv, PEXIT, Tout, GrpRelPower(1:nzth, ixy), CoolantTH(ixy), ThVar, ThOpt, PE)
          !ELSE
            CALL SteadyCoolantTH(PowLin, PowLv, PEXIT, Tout, RelPower(1:nzth, ixy), CoolantTH(ixy), ThVar, ThOpt, PE)
          !END IF
        !END IF
        
        IF (abs(CoolantTH(ixy)%TCoolInOut(2, nzth)-BTEMP) .LT. epsm4) THEN
          IF (Master) THEN
            WRITE(MESG, '(A12, I10)') 'Boiling :', ixy
            CALL MESSAGE(io8,FALSE, TRUE, mesg)
          END IF
        END IF
        
        CALL CP_VA(ThInfo%TCoolInOut(1:2, 1:nzth, ixy), CoolantTH(ixy)%TCoolInOut(1:2, 1:nzth), 2, nzth)
        
        ToutAvg = toutavg + Tout
      END DO
      
			CALL SetGTnGapCoolTemp(Core, ThInfo, PE)
    END IF
    ! ------------------------------------------------
    ToutAvg  = ToutAvg / nchannel
    tfmax    = 0.
    tdoplmax = 0.
    
    IF (.NOT. lftemp_ex) THEN
      DO ixy = xyb, xye
        IF (.NOT. CoolantTH(ixy)%lfuel) CYCLE
        
        !Condunction
        RhoU => CoolantTh(ixy)%rhou
        
        !It has problem for parallel excutition
        CALL SetPwShape(Pwshape, Core, FmInfo%Fxr, FmInfo%Power, ixy, nzth, ThVar%npr5, Thvar, ThOpt, PE)
        CALL SetBurnupDist(BurnUpDist, Core, FmInfo, ixy, nzth, ThVar%npr5, Thvar, ThOpt, PE)
        
        !IF (nTracerCntl%lThChConf) THEN
        !  CALL SteadyFuelConduction_ThCh(Core, powlin, PowLv, Tfmax1, RelPower(1:nzth, ixy), PwShape, BurnUpDist, RhoU, FuelTH(ixy), ThVar, ThOpt, nTracerCntl, PE, hGapArray(:,ixy), ixy)
        !ELSE
          !CALL SteadyFuelConduction(powlin, PowLv, Tfmax1, RelPower(1:nzth, ixy), PwShape, BurnUpDist, RhoU, FuelTH(ixy), ThVar, ThOpt, nTracerCntl, PE, hGapArray(:,ixy))
        CALL SteadyFuelConduction(powlin, PowLv, Tfmax1, RelPower(1:nzth, ixy), PwShape, BurnUpDist, RhoU, FuelTH(ixy), ThVar, ThOpt, nTracerCntl, PE)
        !END IF
        
        tfmax = max(tfmax, tfmax1)
        
        DO iz = 1, thvar%nzth
          Tfavg_max = MAX(FuelTh(ixy)%tfuel(THVar%npr5, iz), Tfavg_max)
          tmod_max  = max(coolantth(ixy)%tcool(iz), tmod_max)
        END DO
      END DO

      IF (lfirst) THEN
        CALL init_AA(is_AA,nxy,nzth,ThVar%npr+4)
        lfirst = FALSE
      END IF
      
      IF (is_AA) THEN
        CALL BurnupUpdate(is_reset)
        
        IF (is_reset) THEN
          ncall = 0
          WRITE(MESG, '(A,i3,A)') 'New BURNUP Step - AA has been RESET'
          IF (MASTER) CALL MESSAGE(io8,TRUE,TRUE,MESG)
        END IF
        
        CALL Anderson_acc(ncall,ThInfo,PE)
        
        ncall = ncall + 1
        
        IF (ncall .GT. 1) THEN
          tfmax = 0._8
          
          DO ixy = xyb, xye
            IF (.NOT. CoolantTH(ixy)%lfuel) CYCLE
            
            CALL recalculate_vars(FuelTH(ixy), BurnUpDist, Tfmax1, ThVar, ThOpt, nTracerCntl, PE)
            
            tfmax = max(tfmax, tfmax1)
          END DO
        END IF
      END IF
    ELSE
      DO ixy = xyb, xye
        IF (.NOT. CoolantTH(ixy)%lfuel) CYCLE
        
        DO iz = 1, nzth
          tfmax = max(tfmax, FuelTH(ixy)%tfuel(1, iz))
        END DO
      END DO
    END IF
  END IF
END DO
! ----------------------------------------------------
! Dopler Temp Update
tdoplmax = 0.
n = 0

DO ixy = xyb, xye
  IF (.NOT. CoolantTH(ixy)%lfuel) CYCLE
  
  n = n + 1
  TdoplMax0 = 0
  
  DO iz = 1, nzth
    tdopd = THInfo%Tdop(iz, ixy)

    THInfo%Tdop(iz, ixy) = SQRT(CKELVIN + FuelTH(ixy)%tfuel(npr5, iz))
    TdoplMax0 = max(TdoplMax0, ABS(1._8 - tdopd/THInfo%Tdop(iz, ixy)))
  END DO
  
  tdoplmax = tdoplmax + TdoplMax0 ** 2
END DO
tdoplmax = sqrt(tdoplmax) / n

ThInfo%Tfmax      = Tfmax
ThInfo%TdopChg    = Tdoplmax
ThInfo%TModoutAvg = ToutAvg
! ----------------------------------------------------
WRITE (mesg, 601) tfmax, toutavg, tdoplmax, FALSE
IF (master) CALL message(io8,false,TRUE,mesg)
601 format(2x, "Max Tf=", f7.1, "C,  Avg Outlet Temp", f7.2, "C", ", Dop. Change=", 1p, e9.2, l2, i3)

ThInfo%Tfavg_max  = Tfavg_max
ThInfo%TModoutAvg = tmod_max

! Assign User defined TH Condition
IF (nTracerCntl%lUSErDefTemp) CALL SetUserDefTH(Core, CmInfo, FmInfo, ThInfo, nTRACERCntl, PE)

#ifdef MPI_ENV
CALL MPI_SYNC(Comm)
#endif

Tend = nTracer_dclock(false, false)
TimeChk%ThTime = TimeChk%ThTime + (TEnd - Tbeg)
! ----------------------------------------------------
IF (lSimpleTH) DEALLOCATE (Profile)
!IF (nTRACERCntl%ThCh_mod .GE. 1 .AND. .NOT.nTRACERCntl%lthch_tf) DEALLOCATE (GrpRelPower)
DEALLOCATE (Pwshape)
DEALLOCATE (BurnUpDist)
NULLIFY (PIN)
NULLIFY (Fxr)
NULLIFY (Fxr)
NULLIFY (PinXS)
NULLIFY (FuelTH)
NULLIFY (CoolantTH)
NULLIFY (RelPower)
NULLIFY (RhoU)
! ----------------------------------------------------

END SUBROUTINE SteadyStateTH