module ESCOTTh_mod
use PARAM
use TYPEDEF,              only: CoreInfo_type, ThInfo_Type, PE_TYPE, AsyInfo_Type,    &
                                Asy_Type, Pin_Type, pininfo_type, FuelTh_Type
use CNTL,                 only: nTracerCntl_Type
use SubChCoupling_mod,    only: SubChCode, coupled_maxsteps, coupled_relconv, LBUUPDATED, is_wall, is_asht
USE IOUTIL,               ONLY: message
USE FILES,                ONLY: IO8
USE TH_Mod,               ONLY: ThOpt
implicit none
private
public :: ESCOT_TH, Finalize_ESCOT

type(SubChCode) :: ESCOT
logical :: lfirst=.true., lpin_resolved
integer :: fuel_cond = 1
contains

subroutine ESCOT_TH(Core, ThInfo, nTracerCntl,ToutAvg, is_Tfuel_cal,PE,ng)
#ifdef HAVE_ESCOT  
  use ESCOT_Coupling_Interface,       only: Gen_Main_Deck, ESCOT_Initialize,    &
                                            ESCOT_set_Gd_Frac,                  &
                                            ESCOT_set_rodpowers_W_cm,           &
                                            ESCOT_set_pinresolved_power,        &
                                            ESCOT_set_BU_GWD_tHM,               &
                                            ESCOT_set_pinresolved_burnup,       &
                                            ESCOT_Solve_Standalone,             &
                                            ESCOT_get_coolant_temp,             &
                                            ESCOT_get_coolant_dens,             &
                                            set_relconv_critera,                &
                                            set_coupled_maxsteps,               &
                                            set_AA_info
#endif                                            
  use MPIComm_Mod,      only: MPI_SYNC
  type(CoreInfo_Type) :: Core
  type(ThInfo_Type) :: ThInfo
  type(nTracerCntl_Type) :: nTracerCntl
  type(PE_Type) :: PE
  logical, intent(inout) :: is_Tfuel_cal
  logical :: is_reset
  real*8, intent(inout) :: ToutAvg
  integer, intent(in) :: ng
  integer :: nproc, nrf, MASTER, bpin=0, epin=0
  integer, parameter :: m_depth = 1
  real*8, parameter :: beta_aa = 1._8
  
  Master = PE%CmfdMaster

  if (lfirst) then
    if (fuel_cond == 0) then
      is_Tfuel_cal = .false.
      lpin_resolved = .false.
    else
      is_Tfuel_cal = .true.
      lpin_resolved = THOpt%LPINRESOLV
    endif
    
    if (nTracerCntl%lHex) call Copy_hCoretoCore(Core,bpin,epin)
    
    call ESCOT%preproc(Core,PE,nTracerCntl)
    
    if (is_Tfuel_cal) call ESCOT%init_subpin(Core,PE,lpin_resolved)
    
#ifdef HAVE_ESCOT
    if (nTracerCntl%lHex) then 
      if (PE%master) then     
        call gen_preproc_input_ESCOT_Hex(Core,nTracerCntl,bpin,epin)
        call Gen_Main_Deck 
      endif
    else
      if (PE%master) then     
        call gen_preproc_input_ESCOT(Core,nTracerCntl)
        call Gen_Main_Deck  
      endif
    endif
#ifdef HAVE_MPI 
    call MPI_SYNC(ESCOT%PE%MPI_COMM)
#endif
    if (ESCOT%PE%lSubCh_proc) then
      WRITE(MESG, '(A,i3)') 'nTRACER/ESCOT running with AA -', m_depth
      IF(Master) CALL MESSAGE(io8,FALSE,TRUE,MESG)
      call ESCOT%setAA(set_AA_info,m_depth,beta_aa)
      call ESCOT_Initialize(ESCOT%PE%MPI_COMM)
      if (coupled_maxsteps>0) call set_coupled_maxsteps(coupled_maxsteps)
      call set_relconv_critera(coupled_relconv)
      if (is_Tfuel_cal) call ESCOT%SetGDFrac(ESCOT_set_Gd_Frac)
    endif
#endif    
    if (nTracerCntl%lHex) deallocate(Core%CoreMap,Core%CoreIdx)
    lfirst=.false.
  endif
#ifdef HAVE_ESCOT    
  if (is_Tfuel_cal) then    
    
    call ESCOT%UpdateRodInfo(Core,PE,lpin_resolved,ng,is_reset)

    if (lpin_resolved) then
      call ESCOT%SetPinResolvedPower(ESCOT_set_pinresolved_power)
  
      WRITE(MESG, '(A)') 'nTRACER is setting the Pin-Resolved Power...'
      IF(Master) CALL MESSAGE(io8,FALSE,TRUE,MESG)
    endif
    
    call ESCOT%SetBurnup(lpin_resolved, ESCOT_set_BU_GWD_tHM, ESCOT_set_pinresolved_burnup)
    if (lpin_resolved .and. LBUUPDATED) then
      WRITE(MESG, '(A)') 'nTRACER is Setting the Pin-Resolved Burnup...'
      IF(Master) CALL MESSAGE(io8,FALSE,TRUE,MESG)
    end if
    
  end if
  
  call ESCOT%SetSubChData(Core,ThInfo,ESCOT_set_rodpowers_W_cm)
  
  WRITE(MESG, '(A)') 'nTRACER is setting the Linear Power...'
  IF(Master) CALL MESSAGE(io8,FALSE,TRUE,MESG)
  
  if (.not. lpin_resolved) then
    WRITE(MESG, '(A)') 'nTRACER is setting the Pin-Average Power...'
    IF(Master) CALL MESSAGE(io8,FALSE,TRUE,MESG)
    if (LBUUPDATED) then
      WRITE(MESG, '(A)') 'nTRACER is Setting the Pin-Average Burnup...'
      IF(Master) CALL MESSAGE(io8,FALSE,TRUE,MESG)
    end if
  end if
#ifdef HAVE_MPI    
    call MPI_SYNC(ESCOT%PE%MPI_COMM)  
#endif  
  if (ESCOT%PE%lSubCh_proc) then       

    call ESCOT_Solve_Standalone(is_reset)
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
  use CNTL,             only: nTracerCntl  
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
  
subroutine gen_preproc_input_ESCOT_Hex(Core,nTracerCntl,bpin,epin)
  use geom,             only: NCELLX0, lCoreAng, lEdge, Albedo
  use ioutil,           only: newunit
  use TH_mod,           only: Thvar, ThOpt
  use DEPL_MOD,         only: DeplCntl
  use files,            only: CaseId
  use SubChCoupling_mod,  only: ActiveAsymIndex, Courant, sbch_outop
  use HexData,          only: hAsyTypInfo, hPinInfo, gCel, aoF2F, hCel
  type(CoreInfo_Type) :: Core
  type(nTracerCntl_Type) :: nTracerCntl
  type(AsyInfo_Type), pointer :: AsyInfo(:)
  type(Asy_Type), pointer :: Asy(:)
  type(Pin_Type), pointer :: Pin(:)
  type(pininfo_type), pointer :: PinInfo(:)
  integer :: ESCOTinp, nxpin, nypin, bpin, epin, offset=0
  real*8 :: alphaP=1.0, alphaU=1.0, beta=0.03
  logical :: lwall=.false., lwht=.false.
  logical :: lgt=.false.
  integer :: iz, ixp, iyp, iya, ixa, ixya, m, ref_assem_id, ref_assem_type,     &
             ixy_loc, ixy_glob, iasyType, ixe, ixb, iye, iyb, j, i, mm, centid, &
            ipintype, ixyp
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

    do ixyp= bpin, epin
      if (.not. pin(ixyp)%lfuel) then
        m=m+1  
        gt_loc(m)%x = hPinInfo(ixyp)%ix
        gt_loc(m)%y = hPinInfo(ixyp)%iy
      endif
    enddo
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
  write(ESCOTinp,'(2X,A,2I4)') 'core_dim', max(ESCOT%nxa_act%n,ESCOT%nya_act%n), max(ESCOT%nxa_act%n,ESCOT%nya_act%n)
  do j=1,max(ESCOT%nxa_act%n,ESCOT%nya_act%n)
    do i=1,max(ESCOT%nxa_act%n,ESCOT%nya_act%n)
      write(ESCOTinp,'(I3)',advance='no') ESCOT%core_conf(i,j)
    enddo
    write(ESCOTinp,*)
  enddo
  write(ESCOTinp,'(A)') '! Axial mesh (mm)'
  write(ESCOTinp,'(2X,A,X,2I4)') 'ax_mesh', ESCOT%nz, nz_nonuni
  do iz=1,nz_nonuni
    write(ESCOTinp,'(X,I3,f11.4)') hz_no_node(iz), hz(iz)
  enddo
  
  if(is_wall) lwall = .true.
  if(is_asht) lwht = .true.
  
  write(ESCOTinp,'(A)') '! Shroud Modeling flag'
  write(ESCOTinp,'(2X,A,L3)') 'wall', lwall     
  write(ESCOTinp,'(A)') '! Symmetric option'
  if (lCoreAng == 360) then
    write(ESCOTinp,'(2X,A,I5)') 'sym_opt', lCoreAng
  else
    if (lEdge) then
      write(ESCOTinp,'(2X,A,I5,2X,A)') 'sym_opt', lCoreAng, 'edge'
    else
      write(ESCOTinp,'(2X,A,I5,2X,A)') 'sym_opt', lCoreAng, 'cent'
    endif
  endif
  write(ESCOTinp,*)  

! Write assembly geometry information  
  write(ESCOTinp,'(A)') 'ASSEM'  
  write(ESCOTinp,'(A)') '! Hexagonal geometry'
  write(ESCOTinp,'(2X,A,A3)') 'hex', 'T'
  write(ESCOTinp,'(A)') '! # of rods in a row / SA'
  write(ESCOTinp,'(2X,A,I4)') 'npin', hCel(1)%nPin
  write(ESCOTinp,'(A)') '! # of fuel rod in a SA'
  write(ESCOTinp,'(2X,A,I6)') 'nfuel', Core%nAsyCell-Core%nAsyGT
  write(ESCOTinp,'(A)') '! Rod diam (mm)'
  write(ESCOTinp,'(2X,A,f11.4)') 'pin_dim', ThVar%rw*2._8*1000._8
  write(ESCOTinp,'(A)') '! Pin pitch (mm)'
  write(ESCOTinp,'(2X,A,f11.4)') 'pin_pitch', hCel(1)%pF2F*10._8  
  if (lgt) then
    write(ESCOTinp,'(A)') '! # of Guide tube/water rod in a SA'
    write(ESCOTinp,'(2X,A,I4)') 'ngt', Core%nAsyGT
    write(ESCOTinp,'(A)') '! Guide tube diam (mm)'
    write(ESCOTinp,'(2X,A,f11.4)') 'gt_dim', ThVar%rgto*2._8*1000._8
    write(ESCOTinp,'(A)') '! Guide tube location'
    write(ESCOTinp,'(2X,A)') 'gt_loc'
    do m=1,Core%nAsyGT
      write(ESCOTinp,'(2I4)') gt_loc(m)%x, gt_loc(m)%y
    enddo
  endif
  write(ESCOTinp,'(A)') '! Assembly pitch (mm)'
  !if (Gap_mod .ne. 0.)then
  !  write(ESCOTinp,'(2X,A,3f9.4)') 'asy_pitch', hCel(1)%aiF2F*10._8, aoF2F*10._8,  Gap_mod*10._8   
  !else
    write(ESCOTinp,'(2X,A,2f9.4)') 'asy_pitch', hCel(1)%aiF2F*10._8, aoF2F*10._8
  !endif 
  if (lwht) then
    write(ESCOTinp,'(A)') '! Inter-assembly heat transfer'  
	write(ESCOTinp,'(2X,A,A3)') 'asy_heat', 'T'
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
    if (ThOpt%PinDim(5) /= 0.) then
      write(ESCOTinp,'(2X,A,I4)') 'nfcond', ThVar%npr-1
    else
      write(ESCOTinp,'(2X,A,I4)') 'nfcond', ThVar%npr
    endif
    write(ESCOTinp,'(A)') '! Radius (mm)'
    write(ESCOTinp,'(2X,A,f12.4,6X,A)') 'rfo', ThVar%rs*1000._8, '! Fuel pellet radius'
    if (ThOpt%PinDim(5) /= 0.) write(ESCOTinp,'(2X,A,f12.4,6X,A)') 'rfi', ThOpt%PinDim(5), '! Fuel pellet inner radius'
    write(ESCOTinp,'(2X,A,f12.4,6X,A)') 'rci', (ThVar%rw-ThVar%tw)*1000., '! Cladding inner radius'      
    write(ESCOTinp,'(2X,A,f11.4,6X,A)') 'rciw', (ThVar%rgti)*1000., '! Guide tube inner radius'           
    write(ESCOTinp,'(A)') '! Correlation model'
    write(ESCOTinp,'(2X,A,I4,6X,A)') 'mode_met', 1,'!FUEL and Clad, 1: FRAPCON/BISON correlations  2: User Input Data'
    write(ESCOTinp,'(A)') '! Fuel composition/burnup/gadolinium/O-M ratio'
    write(ESCOTinp,'(2X,A,f11.4,6X,A)') 'UO2', DeplCntl%NowBurnUp(2), '! Material    Burnup (GWD/tHM)'
    write(ESCOTinp,'(2X,A,I4,6X,A)') 'mode_gap', 2,'!Gap            1: FRAPCON/BISON correlations  2: User Input Data'
    write(ESCOTinp,'(A)') '! Gap thermal conductance (W/m^2/K)'
    write(ESCOTinp,'(2X,A,f11.3)') 'hgap', ThOpt%hgap
    
    allocate(fuelrod_map(2*hAsyTypInfo(1)%nPin-1,2*hAsyTypInfo(1)%nPin-1),source = 1)
    if (lgt) then
      do m = 1,Core%nAsyGT
        fuelrod_map(gt_loc(m)%x, gt_loc(m)%y) = 0
      enddo
    endif
    write(ESCOTinp,'(A)') '! Rod type map !0: not solve 1: UO2 2: UO2Gd 3: MOX 4: MOXGd'
    write(ESCOTinp,'(2X,A)') 'rodtype'
    write(ESCOTinp,'(2X,A)') '{1}'
    nxpin=hAsyTypInfo(1)%nPin; nypin=2*hAsyTypInfo(1)%nPin-1
    do j=1,nypin
      write(ESCOTinp,'(<2*hAsyTypInfo(1)%nPin-1>I3)') (fuelrod_map(i+offset,j), i=1,nxpin)
      if(j < hAsyTypInfo(1)%nPin)then
        nxpin=nxpin+1
      else
        nxpin=nxpin-1
        offset=offset+1
      endif
    enddo
    deallocate(fuelrod_map)
  endif


  write(ESCOTinp,'(A)') '.'
  close(ESCOTinp) 
  
  nullify(Pin); nullify(Asy); nullify(AsyInfo)
  
end subroutine

subroutine Copy_hCoretoCore(Core,bpin,epin)
  use HexData,                 only: hCore, hAsyTypInfo, hAsy
  USE geom,                    ONLY : nZ
  use SubChCoupling_mod,  only: ActiveAsymIndex
  type(CoreInfo_Type) :: Core
  integer :: iya,ixa, ixya, ixy, ixyp, iasy, bpin, epin
  allocate(Core%CoreMap(Core%nxya),Core%CoreIdx(Core%nxa, Core%nya))
  
  ixya=0; iasy=0
  Core%CoreIdx=0
  do iya=1, Core%nya
    do ixa=1, Core%nxa
      if(hCore(ixa,iya) == 0) cycle
      ixya=ixya+1 	  
      Core%CoreMap(ixya) = hCore(ixa,iya)
      Core%CoreIdx(ixa, iya) = ixya
    enddo
  enddo
  
  Core%nAsyCell = hAsyTypInfo(1)%nRodPin(1)
  
  do iasy=1,ixya
    if (hAsy(iasy)%nRodPin ==  Core%nAsyCell) then
      bpin=bpin+1
      epin=epin+hAsy(iasy)%nRodPin
      exit
    endif 
    bpin=bpin+hAsy(iasy)%nTotPin
    epin=epin+hAsy(iasy)%nTotPin
  enddo    
   
  do ixyp=bpin, epin
    if (.not. Core%pin(ixyp)%lfuel) then
      Core%nAsyGT = Core%nAsyGT+1
    endif
  enddo
  
end subroutine

end module