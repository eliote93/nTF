module ESCOTTh_mod
use PARAM
use TYPEDEF,          only: CoreInfo_type, ThInfo_Type, PE_TYPE, AsyInfo_Type,    &
                            Asy_Type, Pin_Type
use CNTL,             only: nTracerCntl_Type
use SubChCoupling_mod,    only: SubChCode, last_TH, coupled_maxsteps, coupled_relconv
implicit none
private
public :: ESCOT_TH

type(SubChCode) :: ESCOT
logical :: lfirst=.true.
contains

subroutine ESCOT_TH(Core, ThInfo, nTracerCntl,ToutAvg, PE)
#ifdef HAVE_ESCOT  
  use ESCOT_Coupling_Interface,       only: Gen_Main_Deck, ESCOT_Initialize,    &
                                            ESCOT_set_rodpowers_W_cm,           &
                                            ESCOT_Solve_Standalone,             &
                                            ESCOT_get_coolant_temp,             &
                                            ESCOT_get_coolant_dens,             &
                                            ESCOT_Edit,   ESCOT_cleanup,        &
                                            ESCOT_get_coolant_temp_exit,        &
                                            set_relconv_critera,                &
                                            set_coupled_maxsteps
#endif                                            
  use MPIComm_Mod,      only: MPI_SYNC
  type(CoreInfo_Type) :: Core
  type(ThInfo_Type) :: ThInfo
  type(nTracerCntl_Type) :: nTracerCntl
  type(PE_Type) :: PE
  real*8 :: ToutAvg
  integer :: nproc
  
  if (lfirst) then
    call ESCOT%preproc(Core,PE)
#ifdef HAVE_ESCOT        
    if (PE%master) then
      call gen_preproc_input_ESCOT(Core,nTracerCntl)
      call Gen_Main_Deck  
    endif
    call MPI_SYNC(ESCOT%PE%MPI_COMM)
    if (ESCOT%PE%lSubCh_proc) then
      call ESCOT_Initialize(ESCOT%PE%MPI_COMM)
      if (coupled_maxsteps>0) call set_coupled_maxsteps(coupled_maxsteps)
      if (coupled_relconv>0.) call set_relconv_critera(coupled_relconv)
    endif
#endif    
    lfirst=.false.
  endif
#ifdef HAVE_ESCOT    
  if (ESCOT%PE%lSubCh_proc) then
    call ESCOT%SetSubChData(Core,ThInfo,ESCOT_set_rodpowers_W_cm)

    call ESCOT_Solve_Standalone

    call MPI_SYNC(ESCOT%PE%MPI_COMM)  

    call ESCOT%GetSubChData(Core,ThInfo,ESCOT_get_coolant_temp,ESCOT_get_coolant_dens)

    call ESCOT_get_coolant_temp_exit(ToutAvg)
  endif  

#endif  

  call ESCOT%Bcast_SubChOut(ThInfo,PE)
  
  if (Core%lgap) then
    call ESCOT%SetAsyGapCoolinfo(Core,ThInfo)
  endif
  
  if (last_TH) then
    if (ESCOT%PE%lSubCh_proc) then
#ifdef HAVE_ESCOT    
    call ESCOT_Edit
    call ESCOT_cleanup
#endif
    endif
    call ESCOT%Finalize()
  endif
  
end subroutine

subroutine gen_preproc_input_ESCOT(Core,nTracerCntl)
  use geom,             only: NCELLX0
  use ioutil,           only: newunit
  use TH_mod,           only: Thvar
  use SubChCoupling_mod,  only: ActiveAsymIndex
  type(CoreInfo_Type) :: Core
  type(nTracerCntl_Type) :: nTracerCntl
  type(AsyInfo_Type), pointer :: AsyInfo(:)
  type(Asy_Type), pointer :: Asy(:)
  type(Pin_Type), pointer :: Pin(:)
  integer :: ESCOTinp
  real*8 :: Courant_number=0.6, alphaP=1.0, alphaU=1.0, beta=0.05
  logical :: lwall=.false. 
  logical :: lgt=.false.
  integer :: iz, ixp, iyp, iya, ixa, ixya, m, ref_assem_id, ref_assem_type,     &
             ixy_loc, ixy_glob, iasyType, ixe, ixb, iye, iyb, j, i
  integer :: nz_nonuni
  integer, allocatable :: hz_no_node(:)
  real*8,allocatable :: hz(:)
  type gtinfo_type
    integer :: x=0, y=0
  endtype
  type(gtinfo_type), allocatable :: gt_loc(:)


  AsyInfo=>Core%AsyInfo
  Asy => Core%Asy
  Pin => Core%Pin

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
        if (.not. (Asy(ixya)%lCentX .or. Asy(ixya)%lCentY .or. Asy(ixya)%lCentXY)) then
          ref_assem_id=ixya
          ref_assem_type=iasyType
          exit sr_ref
        else
!Make Error Cannot find Ref assembly          
        stop 'Cannot find Ref assembly for Guide tube modeling'
        endif
      enddo
    enddo sr_ref

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

    if (Core%lgap) then
      do m = 1, Core%nAsyGT
        gt_loc(m)%x = gt_loc(m)%x-1
        gt_loc(m)%y = gt_loc(m)%y-1
      enddo
    endif
  endif
  
  
  open(newunit(ESCOTinp),file='preproc.inp', status='unknown')

! Write Control Variables  
  write(ESCOTinp,'(A)') 'CNTL'
  write(ESCOTinp,'(A)') '! Total mass flow rate (kg/sec)'
  write(ESCOTinp,'(2X,A,X,f12.4)') 'mdot', nTracerCntl%fmdotfa*ESCOT%nasy_act
  write(ESCOTinp,'(A)') '! Outlet pressure (MPa)'
  write(ESCOTinp,'(2X,A,X,f12.4)') 'pexit', nTracerCntl%Pexit*1.E-6_8
  write(ESCOTinp,'(A)') '! Inlet temp (C)'
  write(ESCOTinp,'(2X,A,X,f12.4)') 'Tin', nTracerCntl%TempInlet
  write(ESCOTinp,'(2X,A,X,f12.4)') 'Courant', Courant_number
  write(ESCOTinp,'(2X,A,2(X,f11.4))') 'alpha', alphaP, alphaU
  write(ESCOTinp,'(2X,A,X,f11.4)') 'beta', beta
  write(ESCOTinp,'(2X,A,L3)') 'lsteady', .true.
  if (ESCOT%run_parallel) then
    write(ESCOTinp,'(2X,A,L3)') 'parallel', .true.
  endif
  write(ESCOTinp,*)

! Write Core geometry information  
  write(ESCOTinp,'(A)') 'GEO'
  write(ESCOTinp,'(A)') '! # of Fuel Assembly Type'
  write(ESCOTinp,'(2X,A,I3)') 'nasy_type', 1
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
  write(ESCOTinp,'(A)') '! Shroud Modeling flag'
  write(ESCOTinp,'(2X,A,L3)') 'wall', lwall
  write(ESCOTinp,*)  

! Write assembly geometry information  
  write(ESCOTinp,'(A)') 'ASSEM'  
  write(ESCOTinp,'(A)') '! # of rods in a row / SA'
  write(ESCOTinp,'(2X,A,I4)') 'npin', NCELLX0
  write(ESCOTinp,'(A)') '! # of fuel rod in a SA'
  write(ESCOTinp,'(2X,A,I6)') 'nfuel', Core%nAsyCell-Core%nAsyGT
  write(ESCOTinp,'(A)') '! # of Guide tube/water rod in a SA'
  write(ESCOTinp,'(2X,A,I4)') 'ngt', Core%nAsyGT  
  write(ESCOTinp,'(A)') '! Rod diam (mm)'
  write(ESCOTinp,'(2X,A,f11.4)') 'pin_dim', ThVar%rw*2._8*1000._8
  write(ESCOTinp,'(A)') '! Pin pitch (mm)'
  write(ESCOTinp,'(2X,A,f11.4)') 'pin_pitch', ThVar%ChannelPitch*1000._8  
  if (lgt) then
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
  write(ESCOTinp,'(A)') '.'
  close(ESCOTinp) 
  
  nullify(Pin); nullify(Asy); nullify(AsyInfo)
  
end subroutine


  
end module