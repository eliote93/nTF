module ESCOTTh_mod
use PARAM
use TYPEDEF,          only: CoreInfo_type, ThInfo_Type, PE_TYPE, AsyInfo_Type,    &
                            Asy_Type, Pin_Type, pininfo_type, FuelTh_Type
use CNTL,             only: nTracerCntl_Type
use SubChCoupling_mod,    only: SubChCode, last_TH, coupled_maxsteps, coupled_relconv
implicit none
private
public :: ESCOT_TH, Finalize_ESCOT

type(SubChCode) :: ESCOT
logical :: lfirst=.true., lfirsttemp = .TRUE.
integer :: fuel_cond = 0
contains

#ifdef __PGI

SUBROUTINE ESCOT_TH(Core, ThInfo, nTracerCntl, ToutAvg, is_Tfuel_cal, PE)
#ifdef HAVE_ESCOT
USE ESCOT_COUPLING_CUDA,    ONLY : ESCOT_GenMainDeck => GenMainDeck, &
                                   ESCOT_Initialize => Initialize, &
                                   ESCOT_SetRodPowers => SetRodPowers, &
                                   ESCOT_SolveStandalone => SolveStandAlone, &
                                   ESCOT_GetCoolantTemp => GetCoolantTemp, &
                                   ESCOT_GetCoolantDens => GetCoolantDens, &
                                   ESCOT_SetRelconvCrit => SetRelativeConvergenceCriteria, &
                                   ESCOT_SetCoupledMaxStpes => SetCoupledMaxSteps, &
                                   ESCOT_SetSomething => SetSomething
#endif
USE MPIComm_Mod,    ONLY : MPI_SYNC
IMPLICIT NONE

TYPE(CoreInfo_TYPE) :: Core
TYPE(ThInfo_TYPE) :: ThInfo
TYPE(nTracerCntl_TYPE) :: nTracerCntl
REAL(8), INTENT(INOUT) :: ToutAvg
LOGICAL, INTENT(INOUT) :: is_Tfuel_cal
TYPE(PE_TYPE) :: PE

INTEGER :: nProc

IF (lFirst) THEN
  fuel_cond = 0
  IF (fuel_cond .EQ. 0) THEN
    is_Tfuel_cal = .FALSE.
  ELSE
    is_Tfuel_cal = .TRUE.
  END IF

  CALL ESCOT%preproc(Core, PE)

#ifdef HAVE_ESCOT
  IF (PE%master) THEN
    CALL gen_preproc_input_ESCOT(Core, nTracerCntl)
    CALL ESCOT_GenMainDeck()
  END IF

#ifdef HAVE_MPI
  CALL MPI_SYNC(ESCOT%PE%MPI_COMM)
#endif
  IF (ESCOT%PE%lSubCh_proc) THEN
    CALL ESCOT_Initialize(ESCOT%PE%MPI_COMM)
    IF (coupled_maxsteps .GT. 0) CALL ESCOT_SetCoupledMaxStpes(coupled_maxsteps)
    CALL ESCOT_SetRelconvCrit(coupled_relconv)
  END IF
#endif
  lFirst = .FALSE.
END IF

#ifdef HAVE_ESCOT
IF (ESCOT%PE%lSubCh_proc) THEN
  CALL ESCOT%SetSubChData(Core, ThInfo, ESCOT_SetRodPowers)

  CALL ESCOT_SolveStandalone()

#ifdef HAVE_MPI
  CALL MPI_SYNC(ESCOT%PE%MPI_COMM)
#endif
  CALL ESCOT%GetSubChData(Core, ThInfo, ESCOT_GetCoolantTemp, ESCOT_GetCoolantDens)

  CALL ESCOT_GetOtherData(Core, ThInfo, ToutAvg)
END IF
#endif

CALL ESCOT%Bcast_SubChOut(ThInfo, .false., PE)

IF (Core%lGap) THEN
  CALL ESCOT%SetAsyGapCoolinfo(Core, ThInfo)
END IF

END SUBROUTINE

SUBROUTINE Finalize_ESCOT()
#ifdef HAVE_ESCOT
USE ESCOT_COUPLING_CUDA,      ONLY : ESCOT_Edit => Edit, &
                                     ESCOT_Cleanup => Cleanup
#endif
IMPLICIT NONE

IF (ESCOT%PE%lSubCh_proc) THEN
#ifdef HAVE_ESCOT
  CALL ESCOT_Edit()
  CALL ESCOT_Cleanup()
#endif
END IF

CALL ESCOT%Finalize()

END SUBROUTINE

SUBROUTINE ESCOT_getotherdata(Core, ThInfo, ToutAvg)
#ifdef HAVE_ESCOT
USE ESCOT_COUPLING_CUDA,      ONLY : ESCOT_GetExitCoolantTemp => GetExitCoolantTemp, &
                                     ESCOT_GetFuelTemp => GetFuelTemp
#endif
USE TH_Mod,   ONLY : ThVar
IMPLICIT NONE

TYPE(CoreInfo_TYPE) :: Core
TYPE(ThInfo_TYPE) :: ThInfo
REAL(8), INTENT(INOUT) :: ToutAvg

TYPE(FuelTh_TYPE), POINTER :: FuelTh(:) => NULL()

INTEGER :: ixy, iz, id_SubCh

#ifdef HAVE_ESCOT

CALL ESCOT_GetExitCoolantTemp(ToutAvg)

IF (fuel_cond .EQ. 0) RETURN


#endif

END SUBROUTINE

SUBROUTINE gen_preproc_input_ESCOT(Core, nTracerCntl)
USE geom,               ONLY : nCellX0, lCoreAng, lEdge, Albedo
USE ioutil,             ONLY : newunit
USE TH_Mod,             ONLY : ThVar, ThOpt
USE DEPL_MOD,           ONLY : DeplCntl
USE files,              ONLY : CaseId
USE SubChCoupling_mod,  ONLY : ActiveAsymIndex, Courant, sbch_outop
IMPLICIT NONE

TYPE(CoreInfo_TYPE) :: Core
TYPE(nTracerCntl_TYPE) :: nTracerCntl

TYPE(AsyInfo_TYPE), POINTER :: AsyInfo(:)
TYPE(Asy_TYPE), POINTER :: Asy(:)
TYPE(Pin_TYPE), POINTER :: Pin(:)
TYPE(PinInfo_TYPE), POINTER :: PinInfo(:)
INTEGER :: ESCOTinp
REAL(8) :: alphaP = 1.0, alphaU = 1.0, beta = 0.05
LOGICAL :: lWall = .TRUE.
LOGICAL :: lgt = .FALSE.
INTEGER :: iz, ixp, iyp, iya, ixa, ixya, m, ref_assem_id, ref_assem_type
INTEGER :: ixy_loc, ixy_glob, iasyType, ixe, ixb, iye, iyb, j, i, mm, centid
INTEGER :: ipintype
INTEGER :: nz_nonuni
INTEGER, ALLOCATABLE :: hz_no_node(:), gt_loc_temp(:,:), fuelrod_map(:,:)
REAL(8), ALLOCATABLE :: hz(:)
CHARACTER(20) :: fmt

TYPE gtinfo_type
  INTEGER :: x = 0, y = 0
END TYPE
TYPE(gtinfo_type), ALLOCATABLE :: gt_loc(:)

AsyInfo => Core%AsyInfo
Asy => Core%Asy
Pin => Core%Pin
PinInfo => Core%PinInfo

! Calculating Axial node information
nz_nonuni = 1
DO iz = 2, ESCOT%nz
  IF (Core%hz(iz) .NE. Core%hz(iz - 1)) THEN
    nz_nonuni = nz_nonuni + 1
  END IF
END DO

ALLOCATE(hz_no_node(nz_nonuni), hz(nz_nonuni))
hz_no_node = 1; hz = 0.0;
hz(1) = Core%hz(1)
nz_nonuni = 1
DO iz = 2, ESCOT%nz
  IF (Core%hz(iz) .NE. Core%hz(iz - 1)) THEN
    nz_nonuni = nz_nonuni + 1
    hz(nz_nonuni) = Core%hz(iz)
  ELSE
    hz_no_node(nz_nonuni) = hz_no_node(nz_nonuni) + 1
  END IF
END DO

hz = hz * 10.0 ! Convert cm to mm

! Guide Tube Information
IF (Core%nAsyGT .GT. 0) lgt = .TRUE.
IF (lgt) THEN
  ALLOCATE(gt_loc(Core%nAsyGT))
  sr_ref : &
  DO iya = 1, ESCOT%nya
    DO ixa = 1, ESCOT%nxa
      ixya = Core%CoreIdx(ixa,iya)
      IF (ixya .EQ. 0) CYCLE
      iasyType = Core%CoreMap(ixya)
      IF (.NOT. AsyInfo(iasyType)%lFuel) CYCLE
      ref_assem_id = ixya
      ref_assem_type = iasyType

      EXIT sr_ref
      IF (ixa .EQ. ESCOT%nxa .AND. iya .EQ. ESCOT%nya) &
        STOP "Cannot find Ref assembly for Guide tube modeling"
    END DO
  END DO sr_ref

  CALL ActiveAsymIndex(Core, Asy(ref_assem_id), ixb, ixe, iyb, iye)
  m = 0

  DO iyp = iyb, iye
    DO ixp = ixb, ixe
      ixy_loc = AsyInfo(ref_assem_type)%pin2dIdx(ixp,iyp)
      ixy_glob = Asy(ref_assem_id)%GlobalPinIdx(ixy_loc)
      IF (pin(ixy_glob)%lGT) THEN
        m = m + 1
        gt_loc(m)%x = ixp
        gt_loc(m)%y = iyp
      END IF
    END DO
  END DO

  IF (Asy(ref_assem_id)%lCentXY) THEN
    ALLOCATE(gt_loc_temp(Core%nAsyGT,2), source = 0)
    DO i = 1, m
      gt_loc_temp(i,1) = gt_loc(i)%x
      gt_loc_temp(i,2) = gt_loc(i)%y
    END DO

    centid = nCellX0 / 2 + mod(nCellX0, 2)
    DO i = 1, m
      gt_loc(i)%x = gt_loc_temp(i,1) + nCellX0 / 2
      gt_loc(i)%y = gt_loc_temp(i,2) + nCellX0 / 2
    END DO

    mm = m

    DO i = 1, m
      ixy_loc = AsyInfo(ref_assem_type)%pin2dIdx(gt_loc_temp(i,1),gt_loc_temp(i,2))
      ixy_glob = Asy(ref_assem_id)%GlobalPinIdx(ixy_loc)
      ipintype = pin(ixy_glob)%PinType
      IF (pinInfo(ipintype)%lCentXY .OR. pinInfo(ipintype)%lCentX) CYCLE
      mm = mm + 1
      gt_loc(mm)%x = centid - gt_loc_temp(i,1) + 1
      gt_loc(mm)%y = gt_loc(i)%y
    END DO

    DO i = 1, m
      ixy_loc=AsyInfo(ref_assem_type)%Pin2dIdx(gt_loc_temp(i,1),gt_loc_temp(i,2))
      ixy_glob=Asy(ref_assem_id)%GlobalPinIdx(ixy_loc)
      ipintype = pin(ixy_glob)%pintype
      IF (Pininfo(ipintype)%lCentX .OR. Pininfo(ipintype)%lCentXY) CYCLE
      mm = mm + 1
      gt_loc(mm)%x = gt_loc(i)%x
      gt_loc(mm)%y = centid - gt_loc_temp(i,2) + 1
    END DO

    DO i = 1, m
      ixy_loc=AsyInfo(ref_assem_type)%Pin2dIdx(gt_loc_temp(i,1),gt_loc_temp(i,2))
      ixy_glob=Asy(ref_assem_id)%GlobalPinIdx(ixy_loc)
      ipintype = pin(ixy_glob)%pintype
      IF (Pininfo(ipintype)%lCentX .OR. Pininfo(ipintype)%lCentY .OR. Pininfo(ipintype)%lCentXY) &
        CYCLE
      mm = mm + 1
      gt_loc(mm)%x = centid - gt_loc_temp(i,1) + 1
      gt_loc(mm)%y = centid - gt_loc_temp(i,2) + 1
    END DO

  ELSE IF (Core%lGap) THEN
    DO m = 1, Core%nAsyGT
      gt_loc(m)%x = gt_loc(m)%x - 1
      gt_loc(m)%y = gt_loc(m)%y - 1
    END DO
  END IF
END IF

OPEN(newunit(ESCOTinp), FILE = 'preproc.inp', STATUS = 'UNKNOWN')

! Write contrl variables
WRITE(ESCOTinp,'(A)') "CNTL"
WRITE(ESCOTinp,'(A)') "! Case ID"
WRITE(ESCOTinp,'(2X,A,2X,A)') "CaseID", TRIM(CaseId)
WRITE(ESCOTinp,'(A)') "! Total mass flow rate (kg/sec)"
WRITE(ESCOTinp,'(2X,A,X,F12.4)') "mdot", nTracerCntl%fMdotfa * ESCOT%actcore_wt
WRITE(ESCOTinp,'(A)') "! Outlet pressure (MPa)"
WRITE(ESCOTinp,'(2X,A,X,F12.4)') "pexit", nTracerCntl%Pexit * 1.E-6_8
WRITE(ESCOTinp,'(A)') "! Inlet temp (C)"
WRITE(ESCOTinp,'(2X,A,X,F12.4)') "Tin", nTracerCntl%TempInlet
WRITE(ESCOTinp,'(2X,A,X,F12.4)') "Courant", Courant
WRITE(ESCOTinp,'(2X,A,2(X,F11.4))') "alpha", alphaP, alphaU
WRITE(ESCOTinp,'(2X,A,X,F11.4)') "beta", beta
WRITE(ESCOTinp,'(2X,A,L3)') "lsteady", .TRUE.
IF (ESCOT%run_parallel) THEN
  WRITE(ESCOTinp,'(2X,A,2L3)') "parallel", .FALSE., .FALSE.
END IF
WRITE(ESCOTinp,'(2X,A,X,I4)') "fuel_temp", fuel_cond
WRITE(ESCOTinp,'(2X,A,I4)') "out_opt", sbch_outop(1)
IF (sbch_outop(2) == 0) THEN
  WRITE(ESCOTinp,'(2X,A,2X,A)') "3D_OUTPUT", "no"
ELSE
  WRITE(ESCOTinp,'(2X,A,2X,A)') "3D_OUTPUT", "ENSIGHT"
END IF

WRITE(ESCOTinp,*)

! Write core geometry information
WRITE(ESCOTinp,'(A)') "GEO"
WRITE(ESCOTinp,'(A)') "! Core Configuration"
WRITE(ESCOTinp,'(2X,A,2I4)') "core_dim", ESCOT%nxa_act%n, ESCOT%nya_act%n
DO j = 1, ESCOT%nya_act%n
  DO i = 1, ESCOT%nxa_act%n
    WRITE(ESCOTinp,'(I3)', ADVANCE = 'NO') ESCOT%core_conf(i,j)
  END DO
  WRITE(ESCOTinp,*)
END DO
WRITE(ESCOTinp,'(A)') "! Axial mesh (mm)"
WRITE(ESCOTinp,'(2X,A,X,2I4)') 'ax_mesh', ESCOT%nz, nz_nonuni
DO iz = 1, nz_nonuni
  WRITE(ESCOTinp,'(X,I3,F11.4)') hz_no_node(iz), hz(iz)
END DO

IF (ALL(ALBEDO(1:4) == 0.0)) THEN
  lWall = .FALSE.
ELSE
  lWall = .TRUE.
END IF

WRITE(ESCOTinp,'(A)') "! Shroud Modeling flag"
WRITE(ESCOTinp,'(2X,A,L3)') "wall", lWall
WRITE(ESCOTinp,'(A)') "! Symmetric option"
IF (lCoreAng == 360) THEN
  WRITE(ESCOTinp,'(2X,A,I5)') "sym_opt", lCoreAng
ELSE IF (lCoreAng == 90) THEN
  IF (lEdge) THEN
    WRITE(ESCOTinp,'(2X,A,I5,2X,A)') "sym_opt", lCoreAng, "edge"
  ELSE
    WRITE(ESCOTinp,'(2X,A,I5,2X,A)') "sym_opt", lCoreAng, "cent"
  END IF
END IF
WRITE(ESCOTinp,*)

! Write assembly geometry information
WRITE(ESCOTinp,'(A)') "ASSEM"
WRITE(ESCOTinp,'(A)') "! # of rods in a row / SA"
WRITE(ESCOTinp,'(2X,A,I4)') "npin", nCellX0
WRITE(ESCOTinp,'(A)') "! # of fuel rod in a SA"
WRITE(ESCOTinp,'(2X,A,I6)') "nfuel", Core%nAsyCell - Core%nAsyGT
WRITE(ESCOTinp,'(A)') "! Rod diam (mm)"
WRITE(ESCOTinp,'(2X,A,F11.4)') "pin_dim", ThVar%rw * 2._8 * 1000._8
WRITE(ESCOTinp,'(A)') "! Pin pitch (mm)"
WRITE(ESCOTinp,'(2X,A,F11.4)') "pin_pitch", ThVar%ChannelPitch * 1000._8
IF (lgt) THEN
  WRITE(ESCOTinp,'(A)') "! # of Guid tube/water rod in a SA"
  IF (ThVar%rgto .GE. SQRT(2.0) * 0.5 * ThVar%ChannelPitch .AND. MOD(nCellX0,2) == 0) THEN
    WRITE(ESCOTinp,'(2X,A,I4,3X,A)') "ngt", Core%nAsyGT, "CE_TYPE"
  ELSE
    WRITE(ESCOTinp,'(2X,A,I4)') "ngt", Core%nAsyGT
  END IF
  WRITE(ESCOTinp,'(A)') "! Guide tube diam (mm)"
  WRITE(ESCOTinp,'(2X,A,F11.4)') "gt_dim", ThVar%rgto * 2._8 * 1000._8
  WRITE(ESCOTinp,'(A)') "! Guide tube location"
  WRITE(ESCOTinp,'(2X,A)') "gt_loc"
  DO m = 1, Core%nAsyGT
    WRITE(ESCOTinp,'(2I4)') gt_loc(m)%x, gt_loc(m)%y
  END DO
END IF
IF (Core%lGap) THEN
  WRITE(ESCOTinp,'(A)') "! Assembly pitch (mm)"
  WRITE(ESCOTinp,'(2X,A,F11.4)') "asy_pitch", ThVar%AsyPitch * 1000._8
END IF
WRITE(ESCOTinp,*)

! Write power information
WRITE(ESCOTinp,'(A)') "POWER"
WRITE(ESCOTinp,'(A)') "! Nominal power of SA (MW)"
WRITE(ESCOTinp,'(2X,A,F11.4)') "asy_power", nTracerCntl%PowerFA * 1.E-6_8
WRITE(ESCOTinp,'(A)') "! Axial power profile option"
WRITE(ESCOTinp,'(2X,A,X,A)') "ax_power", "uniform"
WRITE(ESCOTinp,*)

! Write fuel conduction model
fuel_cond = 0
IF (fuel_cond .GT. 0) THEN
  WRITE(ESCOTinp,'(A)') "FUEL"
  WRITE(ESCOTinp,'(A)') "! Wall temperature correlation"
  WRITE(ESCOTinp,'(2X,A,I4,6X,A)') "nucl_boil", 1,"!1: Zuber-Forster  2:Thom  3: Jens-Lottes"
  WRITE(ESCOTinp,'(A)') "! Ring discretization scheme"
  WRITE(ESCOTinp,'(2X,A,I4,6X,A)') "ring_type", 1,"!1: Equi-space     2: Equi-volume"
  WRITE(ESCOTinp,'(A)') "! # of rings for fuel pellet"
  WRITE(ESCOTinp,'(2X,A,I4)') "nfcond", ThVar%npr
  WRITE(ESCOTinp,'(A)') "! Radius (mm)"
  WRITE(ESCOTinp,'(2X,A,f12.4,6X,A)') "rfo", ThVar%rs*1000._8, "! Fuel pellet radius"
  WRITE(ESCOTinp,'(2X,A,f12.4,6X,A)') "rci", (ThVar%rw-ThVar%tw)*1000., "! Cladding inner radius"
  WRITE(ESCOTinp,'(2X,A,f11.4,6X,A)') "rciw", (ThVar%rgti)*1000., "! Guide tube inner radius"
  WRITE(ESCOTinp,'(A)') "! Correlation model"
  WRITE(ESCOTinp,'(2X,A,I4,6X,A)') "mode_met", 1,"!FUEL and Clad, 1: FRAPCON/BISON correlations  2: User Input Data"
  WRITE(ESCOTinp,'(A)') "! Fuel composition/burnup/gadolinium/O-M ratio"
  WRITE(ESCOTinp,'(2X,A,f11.4,6X,A)') "UO2", DeplCntl%NowBurnUp(2), "! Material    Burnup (GWD/tHM)"
  WRITE(ESCOTinp,'(2X,A,I4,6X,A)') "mode_gap", 2,"!Gap            1: FRAPCON/BISON correlations  2: User Input Data"
  WRITE(ESCOTinp,'(A)') "! Gap thermal conductance (W/m^2/K)"
  WRITE(ESCOTinp,'(2X,A,f11.3)') "hgap", ThOpt%hgap

  ALLOCATE(fuelrod_map(NCELLX0,NCELLX0),source = 1)
  IF (lgt) THEN
    DO m = 1, Core%nAsyGT
      fuelrod_map(gt_loc(m)%x, gt_loc(m)%y) = 0
    END DO
  END IF
  WRITE(ESCOTinp,'(A)') "! Rod type map !0: not solve 1: UO2 2: UO2Gd 3: MOX 4: MOXGd"
  WRITE(ESCOTinp,'(2X,A)') "rodtype"
  WRITE(ESCOTinp,'(2X,A)') "{1}"

  WRITE(fmt,'(I0)') nCellX0

  DO j = 1, nCellX0
    WRITE(ESCOTinp,'('//TRIM(fmt)//'I3)') (fuelrod_map(i,j), i = 1, NCELLX0)
  END DO
  DEALLOCATE(fuelrod_map)
END IF

WRITE(ESCOTinp,'(A)') "."
CLOSE(ESCOTinp)

NULLIFY(Pin); NULLIFY(ASY); NULLIFY(AsyInfo)
END SUBROUTINE

#else

subroutine ESCOT_TH(Core, ThInfo, nTracerCntl,ToutAvg, is_Tfuel_cal,PE)
#ifdef HAVE_ESCOT
  use ESCOT_Coupling_Interface,       only: Gen_Main_Deck, ESCOT_Initialize,    &
                                            ESCOT_set_rodpowers_W_cm,           &
                                            ESCOT_Solve_Standalone,             &
                                            ESCOT_get_coolant_temp,             &
                                            ESCOT_get_coolant_dens,             &
                                            set_relconv_critera,                &
                                            set_coupled_maxsteps
#endif
  use MPIComm_Mod,      only: MPI_SYNC
  type(CoreInfo_Type) :: Core
  type(ThInfo_Type) :: ThInfo
  type(nTracerCntl_Type) :: nTracerCntl
  type(PE_Type) :: PE
  logical, intent(inout) :: is_Tfuel_cal
  real*8, intent(inout) :: ToutAvg
  integer :: nproc

  if (lfirst) then
    if (fuel_cond == 0) then
      is_Tfuel_cal = .false.
    else
      is_Tfuel_cal = .true.
    endif

    call ESCOT%preproc(Core,PE)

#ifdef HAVE_ESCOT
    if (PE%master) then
      call gen_preproc_input_ESCOT(Core,nTracerCntl)
      call Gen_Main_Deck
    endif
#ifdef HAVE_MPI
    call MPI_SYNC(ESCOT%PE%MPI_COMM)
#endif
    if (ESCOT%PE%lSubCh_proc) then
      call ESCOT_Initialize(ESCOT%PE%MPI_COMM)
      if (coupled_maxsteps>0) call set_coupled_maxsteps(coupled_maxsteps)
      call set_relconv_critera(coupled_relconv)
    endif
#endif
    lfirst=.false.
  endif
#ifdef HAVE_ESCOT
  if (ESCOT%PE%lSubCh_proc) then
    call ESCOT%SetSubChData(Core,ThInfo,ESCOT_set_rodpowers_W_cm)

    call ESCOT_Solve_Standalone
#ifdef HAVE_MPI
    call MPI_SYNC(ESCOT%PE%MPI_COMM)
#endif
    call ESCOT%GetSubChData(Core,ThInfo,ESCOT_get_coolant_temp,ESCOT_get_coolant_dens)

    call ESCOT_getotherdata(Core,ThInfo,ToutAvg)
  endif

#endif

  call ESCOT%Bcast_SubChOut(ThInfo,.true.,PE)

  if (Core%lgap) then
    call ESCOT%SetAsyGapCoolinfo(Core,ThInfo)
  endif

end subroutine

subroutine Finalize_ESCOT()
#ifdef HAVE_ESCOT
  use ESCOT_Coupling_Interface,     only: ESCOT_Edit,   ESCOT_cleanup
#endif
  if (ESCOT%PE%lSubCh_proc) then
#ifdef HAVE_ESCOT
    call ESCOT_Edit
    call ESCOT_cleanup
#endif
  endif
  call ESCOT%Finalize()

end subroutine Finalize_ESCOT

subroutine ESCOT_getotherdata(Core,ThInfo,ToutAvg)
#ifdef HAVE_ESCOT
  use ESCOT_Coupling_Interface,       only: ESCOT_get_coolant_temp_exit,    &
                                            ESCOT_get_fuel_temp
#endif
  use TH_mod,           only: Thvar
  type(CoreInfo_Type) :: Core
  type(ThInfo_Type) :: ThInfo
  type(FuelTh_Type), pointer :: FuelTH(:)=>null()
  real*8, intent(inout) :: ToutAvg
  real*8 :: Tfuel_temp(3)
  integer :: iz, iasy, iasytype, ixy, id_SubCh

#ifdef HAVE_ESCOT
  call ESCOT_get_coolant_temp_exit(ToutAvg)

  if (fuel_cond == 0) return

  FuelTH => ThInfo%FuelTH

  do iz=1,ESCOT%nz
    do ixy = 1, ESCOT%nxy
      id_SubCh=ESCOT%id_nT_to_SubCh(ixy)
      if(id_SubCh == 0) cycle
      call ESCOT_get_fuel_temp(id_SubCh,iz,ThInfo%Tfvol(1:Thvar%npr2,iz,ixy),Tfuel_temp)
      if (.not. FuelTH(ixy)%lFuel) cycle
      FuelTH(ixy)%Tfuel(1,iz) = Tfuel_temp(1) - CKELVIN
      FuelTH(ixy)%Tfuel(Thvar%npr1,iz) = Tfuel_temp(2) - CKELVIN
      FuelTH(ixy)%Tfuel(Thvar%npr5,iz) = Tfuel_temp(3) - CKELVIN
    enddo
  enddo

  nullify(FuelTH)
#endif
end subroutine

subroutine gen_preproc_input_ESCOT(Core,nTracerCntl)
  use geom,             only: NCELLX0, lCoreAng, lEdge, Albedo
  use ioutil,           only: newunit
  use TH_mod,           only: Thvar, ThOpt
  use DEPL_MOD,         only: DeplCntl
  use files,            only: CaseId
  use SubChCoupling_mod,  only: ActiveAsymIndex, Courant, sbch_outop
  type(CoreInfo_Type) :: Core
  type(nTracerCntl_Type) :: nTracerCntl
  type(AsyInfo_Type), pointer :: AsyInfo(:)
  type(Asy_Type), pointer :: Asy(:)
  type(Pin_Type), pointer :: Pin(:)
  type(pininfo_type), pointer :: PinInfo(:)
  integer :: ESCOTinp
  real*8 :: alphaP=1.0, alphaU=1.0, beta=0.05
  logical :: lwall=.true.
  logical :: lgt=.false.
  integer :: iz, ixp, iyp, iya, ixa, ixya, m, ref_assem_id, ref_assem_type,     &
             ixy_loc, ixy_glob, iasyType, ixe, ixb, iye, iyb, j, i, mm, centid, &
            ipintype
  integer :: nz_nonuni
  integer, allocatable :: hz_no_node(:), gt_loc_temp(:,:), fuelrod_map(:,:)
  real*8,allocatable :: hz(:)
  type gtinfo_type
    integer :: x=0, y=0
  endtype
  type(gtinfo_type), allocatable :: gt_loc(:)


  AsyInfo=>Core%AsyInfo
  Asy => Core%Asy
  Pin => Core%Pin
  Pininfo => Core%Pininfo

! Calculating Axial node information
  nz_nonuni=1
  do iz=2,ESCOT%nz
    if (Core%hz(iz)/=Core%hz(iz-1)) then
      nz_nonuni=nz_nonuni+1
    endif
  enddo

  allocate(hz_no_node(nz_nonuni),hz(nz_nonuni))
  hz_no_node=1; hz=0.
  hz(1)=Core%hz(1)
  nz_nonuni=1
  do iz=2,ESCOT%nz
    if (Core%hz(iz)/=Core%hz(iz-1)) then
      nz_nonuni=nz_nonuni+1
      hz(nz_nonuni)=Core%hz(iz)
    else
      hz_no_node(nz_nonuni)=hz_no_node(nz_nonuni)+1
    endif
  enddo

  hz=hz*10.   !Convert cm to mm

! Guide Tube Information
  if (Core%nAsyGT>0) lgt=.true.
  if (lgt) then
    allocate(gt_loc(Core%nAsyGT))
     sr_ref: do iya=1,ESCOT%nya
       do ixa=1,ESCOT%nxa
         ixya = Core%CoreIdx(ixa,iya)
         if (ixya == 0) cycle
         iasyType = Core%CoreMap(ixya)
         if (.not. AsyInfo(iasyType)%lFuel) cycle
         ref_assem_id=ixya
         ref_assem_type=iasyType
         exit sr_ref
 !Make Error Cannot find Ref assembly
         if (ixa == ESCOT%nxa .and. iya == ESCOT%nya) stop 'Cannot find Ref assembly for Guide tube modeling'
       enddo
     enddo sr_ref
    !ref_assem_id = 1
    call ActiveAsymIndex(Core, Asy(ref_assem_id),ixb, ixe, iyb, iye)
    m=0

    do iyp=iyb, iye
      do ixp=ixb, ixe
        ixy_loc=AsyInfo(ref_assem_type)%Pin2dIdx(ixp,iyp)
        ixy_glob=Asy(ref_assem_id)%GlobalPinIdx(ixy_loc)
        if (pin(ixy_glob)%lGT) then
          m=m+1
          gt_loc(m)%x = ixp
          gt_loc(m)%y = iyp
        endif
      enddo
    enddo

    if (Asy(ref_assem_id)%lCentXY) then       ! Quarter cut reference assembly case
      allocate(gt_loc_temp(Core%nAsyGT,2), source = 0)
      do i=1,m
        gt_loc_temp(i,1) = gt_loc(i)%x
        gt_loc_temp(i,2) = gt_loc(i)%y
      enddo

      centid = NCELLX0/2 + mod(NCELLX0,2)
      do i=1,m
        gt_loc(i)%x = gt_loc_temp(i,1) + NCELLX0/2
        gt_loc(i)%y = gt_loc_temp(i,2) + NCELLX0/2
      enddo

      mm = m

      do i=1,m
        ixy_loc=AsyInfo(ref_assem_type)%Pin2dIdx(gt_loc_temp(i,1),gt_loc_temp(i,2))
        ixy_glob=Asy(ref_assem_id)%GlobalPinIdx(ixy_loc)
        ipintype = pin(ixy_glob)%pintype
        if (Pininfo(ipintype)%lCentY .or. Pininfo(ipintype)%lCentXY) cycle
        mm = mm + 1
        gt_loc(mm)%x = centid - gt_loc_temp(i,1) + 1
        gt_loc(mm)%y = gt_loc(i)%y
      enddo

      do i=1,m
        ixy_loc=AsyInfo(ref_assem_type)%Pin2dIdx(gt_loc_temp(i,1),gt_loc_temp(i,2))
        ixy_glob=Asy(ref_assem_id)%GlobalPinIdx(ixy_loc)
        ipintype = pin(ixy_glob)%pintype
        if (Pininfo(ipintype)%lCentX .or. Pininfo(ipintype)%lCentXY) cycle
        mm = mm + 1
        gt_loc(mm)%x = gt_loc(i)%x
        gt_loc(mm)%y = centid - gt_loc_temp(i,2) + 1
      enddo

      do i=1,m
        ixy_loc=AsyInfo(ref_assem_type)%Pin2dIdx(gt_loc_temp(i,1),gt_loc_temp(i,2))
        ixy_glob=Asy(ref_assem_id)%GlobalPinIdx(ixy_loc)
        ipintype = pin(ixy_glob)%pintype
        if (Pininfo(ipintype)%lCentX .or. Pininfo(ipintype)%lCentY .or. Pininfo(ipintype)%lCentXY) cycle
        mm = mm + 1
        gt_loc(mm)%x = centid - gt_loc_temp(i,1) + 1
        gt_loc(mm)%y = centid - gt_loc_temp(i,2) + 1
      enddo

    elseif (Core%lgap) then
      do m = 1, Core%nAsyGT
        gt_loc(m)%x = gt_loc(m)%x-1
        gt_loc(m)%y = gt_loc(m)%y-1
      enddo
    endif
  endif


  open(newunit(ESCOTinp),file='preproc.inp', status='unknown')

! Write Control Variables
  write(ESCOTinp,'(A)') 'CNTL'
  write(ESCOTinp,'(A)') '!Case ID'
  write(ESCOTinp,'(2X,A,2X,A)') 'CaseID', trim(CaseId)
  write(ESCOTinp,'(A)') '! Total mass flow rate (kg/sec)'
  write(ESCOTinp,'(2X,A,X,f12.4)') 'mdot', nTracerCntl%fmdotfa*ESCOT%actcore_wt
  write(ESCOTinp,'(A)') '! Outlet pressure (MPa)'
  write(ESCOTinp,'(2X,A,X,f12.4)') 'pexit', nTracerCntl%Pexit*1.E-6_8
  write(ESCOTinp,'(A)') '! Inlet temp (C)'
  write(ESCOTinp,'(2X,A,X,f12.4)') 'Tin', nTracerCntl%TempInlet
  write(ESCOTinp,'(2X,A,X,f12.4)') 'Courant', Courant
  write(ESCOTinp,'(2X,A,2(X,f11.4))') 'alpha', alphaP, alphaU
  write(ESCOTinp,'(2X,A,X,f11.4)') 'beta', beta
  write(ESCOTinp,'(2X,A,L3)') 'lsteady', .true.
  if (ESCOT%run_parallel) then
    write(ESCOTinp,'(2X,A,2L3)') 'parallel', .false., .true.
  endif
  write(ESCOTinp,'(2X,A,X,I4)') 'fuel_temp', fuel_cond
  write(ESCOTinp,'(2X,A,I4)') 'out_opt', sbch_outop(1)
  if (sbch_outop(2) == 0 ) then
    write(ESCOTinp,'(2X,A,2X,A)') '3D_OUTPUT', 'no'
  else
    write(ESCOTinp,'(2X,A,2X,A)') '3D_OUTPUT', 'ENSIGHT'
  endif

  write(ESCOTinp,*)

! Write Core geometry information
  write(ESCOTinp,'(A)') 'GEO'
  write(ESCOTinp,'(A)') '! Core Configuration'
  write(ESCOTinp,'(2X,A,2I4)') 'core_dim', ESCOT%nxa_act%n, ESCOT%nya_act%n
  do j=1,ESCOT%nya_act%n
    do i=1,ESCOT%nxa_act%n
      write(ESCOTinp,'(I3)',advance='no') ESCOT%core_conf(i,j)
    enddo
    write(ESCOTinp,*)
  enddo
  write(ESCOTinp,'(A)') '! Axial mesh (mm)'
  write(ESCOTinp,'(2X,A,X,2I4)') 'ax_mesh', ESCOT%nz, nz_nonuni
  do iz=1,nz_nonuni
    write(ESCOTinp,'(X,I3,f11.4)') hz_no_node(iz), hz(iz)
  enddo

  if (all(ALBEDO(1:4)==0.)) then
    lwall =.false.
  else
    lwall =.true.
  endif

  write(ESCOTinp,'(A)') '! Shroud Modeling flag'
  write(ESCOTinp,'(2X,A,L3)') 'wall', lwall
  write(ESCOTinp,'(A)') '! Symmetric option'
  if (lCoreAng == 360) then
    write(ESCOTinp,'(2X,A,I5)') 'sym_opt', lCoreAng
  elseif (lCoreAng == 90) then
    if (lEdge) then
      write(ESCOTinp,'(2X,A,I5,2X,A)') 'sym_opt', lCoreAng, 'edge'
    else
      write(ESCOTinp,'(2X,A,I5,2X,A)') 'sym_opt', lCoreAng, 'cent'
    endif
  endif
  write(ESCOTinp,*)

! Write assembly geometry information
  write(ESCOTinp,'(A)') 'ASSEM'
  write(ESCOTinp,'(A)') '! # of rods in a row / SA'
  write(ESCOTinp,'(2X,A,I4)') 'npin', NCELLX0
  write(ESCOTinp,'(A)') '! # of fuel rod in a SA'
  write(ESCOTinp,'(2X,A,I6)') 'nfuel', Core%nAsyCell-Core%nAsyGT
  write(ESCOTinp,'(A)') '! Rod diam (mm)'
  write(ESCOTinp,'(2X,A,f11.4)') 'pin_dim', ThVar%rw*2._8*1000._8
  write(ESCOTinp,'(A)') '! Pin pitch (mm)'
  write(ESCOTinp,'(2X,A,f11.4)') 'pin_pitch', ThVar%ChannelPitch*1000._8
  if (lgt) then
    write(ESCOTinp,'(A)') '! # of Guide tube/water rod in a SA'
    if (ThVar%rgto >= sqrt(2.)*0.5*ThVar%ChannelPitch .and. mod(NCELLX0,2) ==0 ) then
      write(ESCOTinp,'(2X,A,I4,3X,A)') 'ngt', Core%nAsyGT, "CE_TYPE"
    else
      write(ESCOTinp,'(2X,A,I4)') 'ngt', Core%nAsyGT
    endif
    write(ESCOTinp,'(A)') '! Guide tube diam (mm)'
    write(ESCOTinp,'(2X,A,f11.4)') 'gt_dim', ThVar%rgto*2._8*1000._8
    write(ESCOTinp,'(A)') '! Guide tube location'
    write(ESCOTinp,'(2X,A)') 'gt_loc'
    do m=1,Core%nAsyGT
      write(ESCOTinp,'(2I4)') gt_loc(m)%x, gt_loc(m)%y
    enddo
  endif
  if (Core%lgap) then
    write(ESCOTinp,'(A)') '! Assembly pitch (mm)'
    write(ESCOTinp,'(2X,A,f11.4)') 'asy_pitch', ThVar%AsyPitch*1000._8
  endif
  write(ESCOTinp,*)

! Write Power information
  write(ESCOTinp,'(A)') 'POWER'
  write(ESCOTinp,'(A)') '! Nominal power of SA (MW)'
  write(ESCOTinp,'(2X,A,f11.4)') 'asy_power', nTracerCntl%PowerFA*1.E-6_8
  write(ESCOTinp,'(A)') '! Axial power profile option'
  write(ESCOTinp,'(2X,A,X, A)') 'ax_power', 'uniform'
  write(ESCOTinp,*)

! Write Fuel conduction model
  if (fuel_cond>0) then
    write(ESCOTinp,'(A)') 'FUEL'
    write(ESCOTinp,'(A)') '! Wall temperature correlation'
    write(ESCOTinp,'(2X,A,I4,6X,A)') 'nucl_boil', 1,'!1: Zuber-Forster  2:Thom  3: Jens-Lottes'
    write(ESCOTinp,'(A)') '! Ring discretization scheme'
    write(ESCOTinp,'(2X,A,I4,6X,A)') 'ring_type', 1,'!1: Equi-space     2: Equi-volume'
    write(ESCOTinp,'(A)') '! # of rings for fuel pellet'
    write(ESCOTinp,'(2X,A,I4)') 'nfcond', ThVar%npr
    write(ESCOTinp,'(A)') '! Radius (mm)'
    write(ESCOTinp,'(2X,A,f12.4,6X,A)') 'rfo', ThVar%rs*1000._8, '! Fuel pellet radius'
    write(ESCOTinp,'(2X,A,f12.4,6X,A)') 'rci', (ThVar%rw-ThVar%tw)*1000., '! Cladding inner radius'
    write(ESCOTinp,'(2X,A,f11.4,6X,A)') 'rciw', (ThVar%rgti)*1000., '! Guide tube inner radius'
    write(ESCOTinp,'(A)') '! Correlation model'
    write(ESCOTinp,'(2X,A,I4,6X,A)') 'mode_met', 1,'!FUEL and Clad, 1: FRAPCON/BISON correlations  2: User Input Data'
    write(ESCOTinp,'(A)') '! Fuel composition/burnup/gadolinium/O-M ratio'
    write(ESCOTinp,'(2X,A,f11.4,6X,A)') 'UO2', DeplCntl%NowBurnUp(2), '! Material    Burnup (GWD/tHM)'
    write(ESCOTinp,'(2X,A,I4,6X,A)') 'mode_gap', 2,'!Gap            1: FRAPCON/BISON correlations  2: User Input Data'
    write(ESCOTinp,'(A)') '! Gap thermal conductance (W/m^2/K)'
    write(ESCOTinp,'(2X,A,f11.3)') 'hgap', ThOpt%hgap

    allocate(fuelrod_map(NCELLX0,NCELLX0),source = 1)
    if (lgt) then
      do m = 1,Core%nAsyGT
        fuelrod_map(gt_loc(m)%x, gt_loc(m)%y) = 0
      enddo
    endif
    write(ESCOTinp,'(A)') '! Rod type map !0: not solve 1: UO2 2: UO2Gd 3: MOX 4: MOXGd'
    write(ESCOTinp,'(2X,A)') 'rodtype'
    write(ESCOTinp,'(2X,A)') '{1}'
    do j=1,NCELLX0
      write(ESCOTinp,'(<NCELLX0>I3)') (fuelrod_map(i,j), i=1,NCELLX0)
    enddo
    deallocate(fuelrod_map)
  endif


  write(ESCOTinp,'(A)') '.'
  close(ESCOTinp)

  nullify(Pin); nullify(Asy); nullify(AsyInfo)

end subroutine

#endif

end module
