module SubChCoupling_mod
use HexData, only : hAsy    
use PARAM
use TYPEDEF,      only: CoreInfo_type,    ThInfo_Type,    Asy_type,       &
                        AsyInfo_type, PE_type 
use CNTL,         only: nTracerCntl_Type
implicit none
private
public :: is_coupled,             &
          CodeName,               &
          SubChCode,              &
          ActiveAsymIndex,        &
          coupled_maxsteps,       &
          coupled_relconv,        &
          Courant,                &
          sbch_outop,             &
          Init_subpin,            &
          SetGDFrac,              &
          UpdateRodInfo,          &
          SetPinResolvedPower,    &
          SetBurnup,              &
          LBUUPDATED,             &
          LGAD_TH,                &
          burnupn,                &
          qvoln,                  &
          gad_fracn,              &
          burnup_avg,             &
          setAA,                  &
          is_wall,                &
          is_asht

logical :: is_coupled=.false., LBUUPDATED, is_wall=.false., is_asht=.false.   
logical, save :: LGAD_TH = .false., LMOX_TH = .false.
character (len=10) :: CodeName
integer :: coupled_maxsteps=-1, MXNXSR = 1
real*8 :: coupled_relconv=0.02, Courant=0.6
integer :: nxy, nz,   &
          sbch_outop(2) = [0 ,0]      !1: Txt output 2: 3D output
real*8 :: PowLin, PowLv
real*8, pointer :: qvoln(:,:,:), ihmmassn(:,:,:), burnupn(:,:,:), gad_fracn(:,:), burnup_avg(:,:)
REAL, pointer :: Power_WATT(:,:)

type (Asy_type), pointer :: Asy(:)
type (AsyInfo_type), pointer :: AsyInfo(:)

type bd_idx
  integer :: ibeg, iend, n
end type bd_idx

type SubChCode
  integer :: nz, nxy, nxya, nxa, nya, nz_act, nasy_act
  integer, allocatable :: id_SubCh_to_nT(:), id_nT_to_SubCh(:), &
                          core_conf(:,:)
  logical :: run_parallel=.false.
  real*8 :: actcore_wt
  character(len=10) :: CodeName
  type(PE_type) :: PE
  type(bd_idx) :: nxa_act, nya_act
  
  contains
  
  procedure :: preproc
  procedure :: SetSubChData
  procedure :: GetSubChData
  procedure :: SetAsyGapCoolinfo
  procedure :: Bcast_SubChOut
  procedure :: Finalize
  procedure :: Init_subpin
  procedure :: SetGDFrac
  procedure :: UpdateRodInfo
  procedure :: SetPinResolvedPower
  procedure :: SetBurnup
  procedure :: SetAA
  
end type SubChCode

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!preprocessing: Index converter, MPIGrouping
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine preproc(self, Core, PE_nT,nTracerCntl)
  class (SubChCode) :: self
  type(CoreInfo_type), intent (in) :: Core
  type(PE_type), intent (in) :: PE_nT
  integer :: iz, iasy, iasytype, id
  integer :: ierr
  integer, allocatable :: proc_list(:)
  type(nTracerCntl_Type) :: nTracerCntl
#ifdef HAVE_MPI
  include 'mpif.h'
#endif

  Asy => Core%Asy
  AsyInfo => Core%AsyInfo
  self%nz=Core%nz
  self%nxy=Core%nxy
  self%nxya=Core%nxya
  self%nxa=Core%nxa
  self%nya=Core%nya
  self%CodeName=CodeName
  self%run_parallel=PE_nT%AxialParallel

  call SetActCoreReg(self, Core)
  if (nTracerCntl%lHex) then
    call IndexConversionHex(self, Core,PE_nT)  
  else    
    call IndexConversion(self, Core)
  endif
  
  if (self%run_parallel) then
    self%PE%myrank = PE_nT%myrank
    select case(self%CodeName)
    case('ESCOT')
      self%PE%nproc=PE_nT%nproc
    case('CTF')
      self%PE%nproc=self%nasy_act
      if (self%PE%nproc<=1) self%run_parallel=.false.
    end select
#ifdef HAVE_MPI
    call MPI_Comm_Group(PE_nT%MPI_CMFD_COMM, self%PE%WorldGroup, ierr)
#endif
    if (PE_nT%myrank < self%PE%nproc) then
      allocate(proc_list(self%PE%nproc))
      do id = 0, self%PE%nproc-1
        proc_list(id+1)=id
      enddo
#ifdef HAVE_MPI
      call MPI_Group_incl(self%PE%WorldGroup, self%PE%nproc, proc_list,    &
                          self%PE%LocalGroup, ierr)
      call MPI_Comm_Create_Group(PE_nT%MPI_CMFD_COMM, self%PE%LocalGroup, &
                                  0, self%PE%MPI_COMM, ierr)
#endif
      self%PE%lSubCh_proc=.true.
    else
      allocate(proc_list(PE_nT%nproc-self%PE%nproc))
      do id=0,size(proc_list)-1
        proc_list(id+1)=self%PE%nproc+id
      enddo
#ifdef HAVE_MPI
      call MPI_Group_incl(self%PE%WorldGroup, PE_nT%nproc-self%PE%nproc,   &
                          proc_list, self%PE%LocalGroup, ierr)
      call MPI_Comm_Create_Group(PE_nT%MPI_CMFD_COMM, self%PE%LocalGroup, &
                                  1, self%PE%MPI_COMM, ierr)
#endif
      self%PE%lSubCh_proc=.false.
    endif
    deallocate(proc_list)
  else
    self%PE%nproc=1
    self%PE%lSubCh_proc=.true.
  endif

end subroutine

subroutine SetSubChData(self,Core,ThInfo, set_rodpowers_W_cm)
  class (SubChCode) :: self
  type(CoreInfo_type), intent (in) :: Core
  type(ThInfo_Type), intent (in) :: ThInfo
  interface 
    subroutine set_rodpowers_W_cm(pin_number,axial_level,rodpower)
      integer,intent(in) :: pin_number, axial_level
      real*8,intent(in) :: rodpower
    end subroutine
  end interface
  integer :: iz, iasy, iasytype, ix, iy, id_Subch, ixy
  real*8 :: qprime

  do iz=1,self%nz
    do ixy = 1, self%nxy
      id_SubCh=self%id_nT_to_SubCh(ixy)
      if(id_SubCh == 0) cycle
      qprime=ThInfo%RelPower(iz, ixy)*ThInfo%PowLin*ThInfo%Powlv/100._8   !W/m -> W/cm
      call set_rodpowers_W_cm(id_SubCh,iz,qprime)
    enddo    
  enddo

end subroutine SetSubChData

subroutine SetGDFrac(self,ESCOT_set_Gd_Frac)
  implicit none
  class (SubChCode) :: self
  interface
    subroutine ESCOT_set_Gd_Frac(pin_number,axial_level,gdfrac)
      integer,intent(in) :: pin_number, axial_level
      real*8,intent(in) :: gdfrac
    end subroutine
  end interface
  integer :: iz, ixy, id_SubCh
  
  if (LGAD_TH) then
    do iz=1,self%nz
      do ixy = 1, self%nxy
        id_SubCh=self%id_nT_to_SubCh(ixy)
        if(id_SubCh == 0) cycle
        call ESCOT_set_Gd_Frac(id_SubCh,iz,gad_fracn(ixy,iz))
      enddo    
    enddo
  else
    return
  end if
  
end subroutine

subroutine setAA(self,set_AA_info,m_depth,beta_aa)
  implicit none
  class (SubChCode) :: self
  interface
    subroutine set_AA_info(m_depth,beta)
      integer, intent(in) :: m_depth
      REAL*8, intent(in) :: beta
    end subroutine
  end interface
  integer, intent(in) :: m_depth
  real*8, intent(in) :: beta_aa
  call set_AA_info(m_depth, beta_aa)
end subroutine
  
subroutine SetPinResolvedPower(self, ESCOT_set_pinresolved_power)
  use TH_mod,           only: Thvar
  implicit none
  class (SubChCode) :: self
  interface
    subroutine ESCOT_set_pinresolved_power(pin_number,axial_level,PinPower)
      integer,intent(in) :: pin_number, axial_level
      real*8 :: PinPower(:)
    end subroutine
  end interface
  integer :: nrf, iz, ixy, id_SubCh

  nrf = ThVar%npr 
  
  do iz=1,self%nz
    do ixy = 1, self%nxy
      id_SubCh=self%id_nT_to_SubCh(ixy)
      if(id_SubCh == 0) cycle
      call ESCOT_set_pinresolved_power(id_SubCh,iz,qvoln(1:nrf,ixy,iz))
    enddo
  enddo
  
end subroutine

subroutine SetBurnup(self, lpin_resolved, ESCOT_set_BU_GWD_tHM, ESCOT_set_pinresolved_burnup)
  use TH_mod,           only: Thvar
  implicit none
  class (SubChCode) :: self
  interface
    subroutine ESCOT_set_BU_GWD_tHM(pin_number,axial_level,bu)
      integer,intent(in) :: pin_number, axial_level
      real*8,intent(in) ::bu
    end subroutine
    subroutine ESCOT_set_pinresolved_burnup(pin_number,axial_level,burnup)
      integer,intent(in) :: pin_number, axial_level
      real*8 :: burnup(:)
    end subroutine
  end interface
  logical, intent(in) :: lpin_resolved
  integer :: nrf, iz, ixy, id_SubCh
  
  if (.not. LBUUPDATED) return
  
  nrf = ThVar%npr
  
  do iz=1,self%nz
    do ixy = 1, self%nxy
      id_SubCh=self%id_nT_to_SubCh(ixy)
      if(id_SubCh == 0) cycle
      if (lpin_resolved) then
        call ESCOT_set_pinresolved_burnup(id_SubCh,iz,burnupn(1:nrf,ixy,iz))
      else
        call ESCOT_set_BU_GWD_tHM(id_SubCh,iz,burnup_avg(ixy,iz))
      end if
    end do
  end do    
  
  end subroutine

subroutine GetSubChData(self,Core,ThInfo,get_coolant_temp,get_coolant_dens, &
                        axial_idx_calib)
  class (SubChCode) :: self
  type(CoreInfo_type), intent (in) :: Core
  type(ThInfo_Type), intent (inout) :: ThInfo
  interface 
  subroutine get_coolant_temp(pin_number,axial_level,temp)
    integer, intent (in) :: pin_number, axial_level
    real*8, intent (inout) :: temp
  end subroutine
  subroutine get_coolant_dens(pin_number,axial_level,density)
    integer, intent (in) :: pin_number, axial_level
    real*8, intent (inout) :: density
  end subroutine
  end interface
  integer, optional, intent(in) :: axial_idx_calib
  integer :: iz, iasy, iasytype, ix, iy, id_Subch, ixy, idum, iz_subch
  real*8 :: qprime

  if(present(axial_idx_calib)) then
    idum=axial_idx_calib
  else
    idum=0
  endif

  do iz=1,self%nz
    do ixy = 1, self%nxy
      id_SubCh=self%id_nT_to_SubCh(ixy)
      if(id_SubCh == 0) cycle
      iz_subch=iz+idum
      call get_coolant_temp(id_SubCh,iz_subch,ThInfo%Tcool(iz,ixy))
      call get_coolant_dens(id_SubCh,iz_subch,ThInfo%DenCool(iz,ixy))
    enddo    
  enddo

end subroutine


subroutine SetAsyGapCoolinfo(self,Core,ThInfo)
  class (SubChCode) :: self
  type(CoreInfo_type), intent (in) :: Core
  type(ThInfo_Type), intent (inout) :: ThInfo
  integer :: iasy, iasytype
  integer :: nx_loc, ny_loc
  integer, allocatable, dimension(:) :: actedge_col_loc, asygap_col_loc, actedge_col_glob,  &
                                    asygap_col_glob, actedge_row_loc, asygap_row_loc,   &
                                    actedge_row_glob, asygap_row_glob
                                    
!Copy T/H info from boundary active cells to gap cells    
  do iasy = 1, self%nxya
    iasytype = Core%CoreMap(iasy)
    if(.not. AsyInfo(iasyType)%lFuel) cycle
    nx_loc = AsyInfo(iasyType)%nx
    ny_loc = AsyInfo(iasyType)%ny
!Copy bottem row edge ells
    allocate(actedge_row_loc(nx_loc), asygap_row_loc(nx_loc),     &
              actedge_row_glob(nx_loc), asygap_row_glob(nx_loc))
    actedge_row_loc = AsyInfo(iasyType)%Pin2dIdx(1:nx_loc,ny_loc-1)
    asygap_row_loc  = AsyInfo(iasyType)%Pin2dIdx(1:nx_loc,ny_loc)
    actedge_row_glob = Asy(iasy)%GlobalPinIdx(actedge_row_loc)
    asygap_row_glob  = Asy(iasy)%GlobalPinIdx(asygap_row_loc)
	
    ThInfo%Tcool(:,asygap_row_glob)=ThInfo%Tcool(:,actedge_row_glob)
    ThInfo%Dencool(:,asygap_row_glob)=ThInfo%Dencool(:,actedge_row_glob)
!Copy left column edge cells
    allocate(actedge_col_loc(ny_loc),asygap_col_loc(ny_loc),      &
              actedge_col_glob(ny_loc), asygap_col_glob(ny_loc))
    actedge_col_loc = AsyInfo(iasyType)%Pin2dIdx(nx_loc-1,1:ny_loc)
    asygap_col_loc = AsyInfo(iasyType)%Pin2dIdx(nx_loc,1:ny_loc)
    actedge_col_glob = Asy(iasy)%GlobalPinIdx(actedge_col_loc)
    asygap_col_glob = Asy(iasy)%GlobalPinIdx(asygap_col_loc)

    ThInfo%Tcool(:,asygap_col_glob)=ThInfo%Tcool(:,actedge_col_glob)
    ThInfo%Dencool(:,asygap_col_glob)=ThInfo%Dencool(:,actedge_col_glob)
    
!Treatment of half or whole assemly region
    if (.not. asy(iasy)%LCENTXY) then
	    if (.not. asy(iasy)%LCENTX) then
!Copy top row edge cells
	      actedge_row_loc = AsyInfo(iasyType)%Pin2dIdx(1:nx_loc,2)
	      asygap_row_loc = AsyInfo(iasyType)%Pin2dIdx(1:nx_loc,1)
	      actedge_row_glob = Asy(iasy)%GlobalPinIdx(actedge_row_loc)
	      asygap_row_glob = Asy(iasy)%GlobalPinIdx(asygap_row_loc)

	      ThInfo%Tcool(:,asygap_row_glob)=ThInfo%Tcool(:,actedge_row_glob)
	      ThInfo%Dencool(:,asygap_row_glob)=ThInfo%Dencool(:,actedge_row_glob)
	    endif
		
	    if(.not. asy(iasy)%LCENTY) then
!Copy right column edge cells      
	      actedge_col_loc = AsyInfo(iasyType)%Pin2dIdx(2,1:ny_loc)
	      asygap_col_loc = AsyInfo(iasyType)%Pin2dIdx(1,1:ny_loc)
	      actedge_col_glob = Asy(iasy)%GlobalPinIdx(actedge_col_loc)
	      asygap_col_glob = Asy(iasy)%GlobalPinIdx(asygap_col_loc)
		
	      ThInfo%Tcool(:,asygap_col_glob)=ThInfo%Tcool(:,actedge_col_glob)
	      ThInfo%Dencool(:,asygap_col_glob)=ThInfo%Dencool(:,actedge_col_glob)
	    endif
    endif
    
		deallocate(actedge_col_loc,asygap_col_loc,actedge_col_glob,         &
                asygap_col_glob, actedge_row_loc, asygap_row_loc,       &
                actedge_row_glob, asygap_row_glob)
  enddo
  
end subroutine

subroutine Bcast_SubChOut(self,ThInfo,is_ftemp_comm,PE)
  use MPIComm_Mod, only : BCAST
  class (SubChCode) :: self
  type(PE_TYPE) :: PE
  type(ThInfo_Type) :: ThInfo
  logical, intent (in) :: is_ftemp_comm
  integer :: comm, i, j, npr2

  comm = PE%MPI_CMFD_COMM

  call BCAST(ThInfo%Tcool, self%nz, self%nxy+2, COMM)
  call BCAST(ThInfo%Dencool, self%nz, self%nxy+2, COMM)
  
  if (.not. is_ftemp_comm) return
  
  npr2 = size(ThInfo%Tfvol,1)
  call BCAST(ThInfo%Tfvol, npr2, self%nz, self%nxy, COMM)
  do i = 1, self%nxy
    if (.not. ThInfo%FuelTH(i)%lFuel) cycle
    call BCAST(ThInfo%FuelTH(i)%Tfuel, npr2+3, self%nz, COMM)
  enddo

end subroutine

subroutine Finalize(self)
  class (SubChCode) :: self
#ifdef HAVE_MPI
  include 'mpif.h'  
#endif  
  integer :: ierr 
  
  nullify(Asy)
  nullify(AsyInfo)
  
  if (.not. self%run_parallel) return
#ifdef HAVE_MPI
  call MPI_Group_Free(self%PE%LocalGroup,ierr)
  call MPI_Group_Free(self%PE%WorldGroup,ierr)
  call MPI_Comm_Free(self%PE%MPI_COMM,ierr)
#endif
end subroutine
!Rule for indexing
!
!  
subroutine IndexConversion(self, Core)
  class (SubChCode) :: self
  type(CoreInfo_type) :: Core
  integer, allocatable :: max_nrod_row(:)
  integer :: total_nrod, iya, ixa, iassem, nrod_row, ixya, iasyType,        &
              nxr, nyr, nxy, ystart, xstart, iy, ix, ixy0, ibx, iex, iby,        &
              iey, id_Subch, id_nT


  allocate(max_nrod_row(self%nya))
  max_nrod_row=0
  total_nrod=0
  do iya = 1, self%nya
    iassem=0
    nrod_row=0
    do ixa = 1, self%nxa
      ixya = Core%CoreIdx(ixa, iya)
      if (ixya == 0) cycle
      iasyType = Core%CoreMap(ixya)
      if (.not. AsyInfo(iasyType)%lFuel) cycle

      call ActiveAsymIndex(Core, Asy(ixya), ibx, iex, iby, iey)

      nxr=iex-ibx+1
      nyr=iey-iby+1
      total_nrod=total_nrod+nxr*nyr
      nrod_row=nrod_row+nxr
      iassem=iassem+1
    enddo
    max_nrod_row(iya)=nrod_row    
  enddo

  allocate(self%id_SubCh_to_nT(total_nrod),self%id_nT_to_SubCh(self%nxy))
  self%id_SubCh_to_nT=0
  self%id_nT_to_SubCh=0

  ystart=0
  do iya = 1,self%nya
    xstart=0
    do ixa = 1,self%nxa
      ixya = Core%CoreIdx(ixa, iya)
      if (ixya == 0) cycle
      iasyType = Core%CoreMap(ixya)
      if (.not. AsyInfo(iasyType)%lFuel) cycle

      call ActiveAsymIndex(Core, Asy(ixya), ibx, iex, iby, iey)
      id_Subch=xstart+ystart
      do iy=iby, iey
        do ix=ibx, iex
          ixy0 = AsyInfo(iasyType)%Pin2dIdx(ix,iy)
          id_nT= Asy(ixya)%GlobalPinIdx(ixy0)
          id_Subch=id_Subch+1
          self%id_SubCh_to_nT(id_Subch)=id_nT
          self%id_nT_to_SubCh(id_nT)=id_Subch
        enddo
        id_Subch=id_Subch+max_nrod_row(iya)-(iex-ibx+1)
      enddo
      xstart=xstart+(iex-ibx+1)
    enddo
    ystart=ystart+max_nrod_row(iya)*(iey-iby+1)
  enddo

end subroutine

subroutine IndexConversionHex(self, Core,PE)
  use HexData, only : hAsy, hPinInfo, hAsyTypInfo
  use CORE_mod,only : THInfo
  class (SubChCode) :: self
  type(CoreInfo_type) :: Core
  integer, allocatable :: max_nrod_row(:)
  integer :: total_nrod, iya, ixa, nrod_row, ixya, iasyType,        &
              nxr, nyr, nxy, ystart, xstart, iy, ix, ixy0, ibx, iex, iby,        &
              iey, id_Subch, id_nT, nCx, id_pin, ixy
  type(PE_type) :: PE
  
  allocate(max_nrod_row(self%nya))
  total_nrod=0
  do iya = 1, self%nya
    nrod_row=0
    do ixa = 1, self%nxa
      ixya = Core%CoreIdx(ixa, iya)
      if (ixya == 0) cycle
      iasyType = Core%CoreMap(ixya)
      if (.not. AsyInfo(iasyType)%lFuel) cycle
  
      total_nrod=total_nrod+hAsy(ixya)%nRodPin
    enddo
  enddo
  
  allocate(self%id_SubCh_to_nT(total_nrod),self%id_nT_to_SubCh(self%nxy))
  self%id_SubCh_to_nT=0
  self%id_nT_to_SubCh=0
  
  id_Subch=0; id_pin=0
  do iya = 1,self%nya
    do ixa = 1,self%nxa
      ixya = Core%CoreIdx(ixa, iya)
      if (ixya == 0) cycle
      iasyType = Core%CoreMap(ixya)
      
      if (.not. AsyInfo(iasyType)%lFuel) then
        id_pin=id_pin+hAsy(ixya)%nTotPin
        cycle
      endif
      
      do ixy=1, hAsy(ixya)%ntotpin
	    id_pin=id_pin+1
        if (hPinInfo(id_pin)%lRod) id_Subch=id_Subch+1
		
        if (hPinInfo(id_pin)%lRod) self%id_SubCh_to_nT(id_Subch)=id_pin
        if (hPinInfo(id_pin)%lRod) self%id_nT_to_SubCh(id_pin)=id_Subch
      enddo
    enddo
  enddo
  
end subroutine

subroutine SetActCoreReg(self,Core)
  class(SubChCode) :: self
  type(CoreInfo_type) :: Core
  integer :: iz, iasy, ixa, iya, iasytype
  
  self%nz_act=0
  self%nasy_act=0
  self%actcore_wt = 0.
  do iz=1, self%nz
    if (.not. Core%lFuelPlane(iz)) cycle
    self%nz_act = self%nz_act +1
  enddo

  do iasy = 1, self%nxya
    iasytype = Core%CoreMap(iasy)
    if (.not. AsyInfo(iasytype)%lFuel) cycle
    self%nasy_act=self%nasy_act+1
    self%actcore_wt = self%actcore_wt + Asy(iasy)%wt 
  enddo

  self%nxa_act%ibeg=self%nxa
  self%nxa_act%iend=1
  
  do iya = 1, self%nya
    left_sweep:do ixa = 1, self%nxa
      iasy = Core%CoreIdx(ixa, iya)
      if (iasy==0) cycle
      iasytype = Core%CoreMap(iasy)
      if (.not. AsyInfo(iasytype)%lFuel) cycle
      self%nxa_act%ibeg=min(self%nxa_act%ibeg,ixa)
      exit left_sweep
    enddo left_sweep

    right_sweep: do ixa = self%nxa, 1, -1
      iasy = Core%CoreIdx(ixa, iya)
      if (iasy==0) cycle
      iasytype = Core%CoreMap(iasy)
      if (.not. AsyInfo(iasytype)%lFuel) cycle
      self%nxa_act%iend=max(self%nxa_act%iend,ixa)
      exit right_sweep
    enddo right_sweep
  enddo
  self%nxa_act%n=self%nxa_act%iend-self%nxa_act%ibeg+1

  self%nya_act%ibeg=self%nya
  self%nya_act%iend=1
  do ixa = 1, self%nxa
    top_sweep: do iya = 1, self%nya
      iasy = Core%CoreIdx(ixa, iya)
      if (iasy==0) cycle
      iasytype = Core%CoreMap(iasy)
      if (.not. AsyInfo(iasytype)%lFuel) cycle
      self%nya_act%ibeg=min(self%nya_act%ibeg,iya)
      exit top_sweep
    enddo top_sweep

    botm_sweep: do iya = self%nya, 1, -1
      iasy = Core%CoreIdx(ixa, iya)
      if (iasy==0) cycle
      iasytype = Core%CoreMap(iasy)
      if (.not. AsyInfo(iasytype)%lFuel) cycle
      self%nya_act%iend=max(self%nya_act%iend,iya)
      exit botm_sweep
    enddo botm_sweep
  enddo  
  self%nya_act%n=self%nya_act%iend-self%nya_act%ibeg+1

  allocate(self%core_conf(max(self%nxa_act%n,self%nya_act%n),max(self%nxa_act%n,self%nya_act%n)))

  self%core_conf=0                          
  do iya=self%nya_act%ibeg, self%nya_act%iend
    do ixa=self%nxa_act%ibeg, self%nxa_act%iend
      iasy = Core%CoreIdx(ixa, iya)
      if (iasy==0) cycle
      iasytype = Core%CoreMap(iasy)
      if (AsyInfo(iasytype)%lFuel) then
        self%core_conf(ixa-self%nxa_act%ibeg+1,iya-self%nya_act%ibeg+1)=1
      endif
    enddo
  enddo
end subroutine

subroutine ActiveAsymIndex(Core, Asy_local, ibx, iex, iby, iey)
  type (CoreInfo_type) :: Core
  type (Asy_type) :: Asy_local
  integer, intent (inout) :: ibx, iex, iby, iey
  integer :: iasyType, ncell_x, ncell_y

  iasyType = Asy_local%asytype
  ncell_x=AsyInfo(iasyType)%nx
  ncell_y=AsyInfo(iasyType)%ny

  if (Core%lgap) then
    if (Asy_local%LCENTXY) then
      ibx=1; iex=ncell_x-1; iby=1; iey=ncell_y-1; 
    elseif (Asy_local%LCENTX) then
      ibx=2; iex=ncell_x-1; iby=1; iey=ncell_y-1
    elseif (Asy_local%LCENTY) then
      ibx=1; iex=ncell_x-1; iby=2; iey=ncell_y-1
    else
      ibx=2; iex=ncell_x-1; iby=2; iey=ncell_y-1
    ENDIF
  else
    ibx=1; iex=ncell_x; iby =1; iey=ncell_y
  endif

end subroutine ActiveAsymIndex

subroutine Init_subpin(self,Core,PE,lpin_res)
  use TH_mod,           only: Thvar
  implicit none
  class(SubChCode) :: self
  type(PE_Type) :: PE
  TYPE(CoreInfo_Type) :: Core
  logical, intent(in) :: lpin_res
  integer :: nrf, nxypin, nz, myzb, myze, nCellType, i, nCoreFsr
  
  nCoreFsr = Core%nCoreFsr
  nrf = ThVar%npr
  nxypin = self%nxy
  nz = self%nz
  myzb = PE%myzb; myze = PE%myze
  MXNXSR = 0
  nCellType = Core%nCellType 
  do i = 1, nCellType
    MXNXSR = MAX(MXNXSR, Core%CellInfo(i)%nFxr)
  end do   
  
  if (lpin_res) then
    if(PE%Master) print*, 'Pin-Resolved Analysis'
    allocate(qvoln(1:nrf, 1:nxypin, 1:nz), source = 0._8)
    allocate(ihmmassn(1:nrf, 1:nxypin, myzb:myze), source = 0._8)
    allocate(burnupn(1:nrf, 1:nxypin, 1:nz), source = 0._8)
    allocate(Power_WATT(1:nCoreFsr,1:nz), source = 0._8)
    call InitIHM(Core,nxypin,myzb,myze,nrf)
  else
    if(PE%Master) print*, 'Pin-Average Analysis'
    allocate(burnup_AVG(1:nxypin, 1:nz), source = 0._8)
  end if
  call Gadolinium_presence(Core,nxypin,myzb,myze,nz)
  !Insert MOX presence
  if (LGAD_TH) call CalcGadFrac(Core,nxypin,myzb,myze,PE%MPI_CMFD_COMM)
  
end subroutine Init_subpin

subroutine InitIHM(Core,npin,myzb,myze,nrf)
  USE Typedef, only : CoreInfo_Type, Cell_Type, Pin_Type,FxrInfo_Type
  USE Core_mod, only : FmInfo
  implicit none
  TYPE(CoreInfo_Type) :: Core
  TYPE(Cell_Type), POINTER :: CellInfo(:)
  TYPE(Pin_Type), POINTER :: Pin(:)
  TYPE(FxrInfo_Type), POINTER :: FXR(:,:)
  integer, intent(in) :: npin,myzb,myze,nrf
  real*8 :: IHM(MXNXSR), vol_tot=0., HZ
  integer :: k, ixy, icel, nLocalFxr, j, FxrIdxSt, LX, IX2, IFUEL
  
  CellInfo => Core%CellInfo
  Pin => Core%Pin
  FXR => FmInfo%fxr
  IHM = 0._8
  
  do k = myzb,myze
    if (.NOT. Core%lFuelPlane(k)) cycle
    HZ = Core%hz(k)
    do ixy = 1, npin
      icel = Pin(ixy)%Cell(k)
      IF(.NOT. CellInfo(icel)%lfuel) CYCLE
      IF(CellInfo(icel)%lRect) CYCLE
      IF(CellInfo(icel)%lCad) CYCLE
      nLocalFxr = CellInfo(icel)%nFxr
      FxrIdxSt = Pin(ixy)%FxrIdxSt
      IFUEL = -99
      vol_tot=0.
        do j = 1,nLocalFxr
          LX = FxrIdxSt + j - 1
          IX2 = nLocalFxr - j + 1
          if (FXR(LX,k)%LFUEL) then
            IFUEL = max(IX2, IFUEL)
            vol_tot = vol_tot + FXR(LX,k)%area*HZ
            IHM(IX2) = FXR(LX,k)%Hmkg0/(FXR(LX,k)%area*HZ)
          end if
        end do
        call transfer_info(icel, nLocalFxr, IHM(1:IFUEL), IFUEL, ihmmassn(1:nrf,IXY,K), nrf, 2,&
                            CellInfo(icel))
        ihmmassn(1:nrf,IXY,K) = ihmmassn(1:nrf,IXY,K)*vol_tot
      end do
  end do
  
end subroutine InitIHM

subroutine Gadolinium_presence(Core,npin,myzb, myze,nz)
  USE Typedef,         only: CoreInfo_Type, Cell_Type, Pin_Type,FxrInfo_Type
  use Material_Mod,    only: Mixture
  use Core_mod,        only: FmInfo
  implicit none
  TYPE(CoreInfo_Type) :: Core
  TYPE(Cell_Type), POINTER :: CellInfo(:)
  TYPE(Pin_Type), POINTER :: Pin(:)
  TYPE(FxrInfo_Type), POINTER :: FXR(:,:)
  integer, intent(in) :: npin, myzb, myze, nz
  integer:: k, ixy, icel, nLocalFxr, FxrIdxSt, j, LX, IM
  
  CellInfo => Core%CellInfo
  Pin => Core%Pin
  FXR => FmInfo%fxr
  
  do k = myzb,myze
    if (.NOT. Core%lFuelPlane(k)) cycle
    do ixy = 1, npin
      icel = Pin(ixy)%Cell(k)
      if (CellInfo(icel)%lfuel) then
        nLocalFxr = CellInfo(icel)%nFxr
        FxrIdxSt = Pin(ixy)%FxrIdxSt
        do j = 1,nLocalFxr
          LX = FxrIdxSt + j - 1
          IM = FXR(LX,K)%IMIX
          if (Mixture(IM)%LGD) then
            LGAD_TH = .TRUE.
            allocate(gad_fracn(1:npin, 1:nz), source = 0._8)
            return
          else
            continue
          end if
        end do
      else
        continue
      end if
    end do
  end do
  
end subroutine

subroutine UpdateRodInfo(Self,Core,PE,lpin_resolved,ng,is_reset)
  use TYPEDEF,          only: FxrInfo_Type
  use TH_mod,           only: Thvar
  USE CNTL,             ONLY: nTracerCntl
  USE CORE_MOD,         ONLY: GroupInfo, FmInfo
  USE MOC_MOD,          ONLY: PowerUpdate_Watt
  implicit none
  class(SubChCode) :: self
  type(PE_Type) :: PE
  TYPE(CoreInfo_Type) :: Core
  TYPE(FxrInfo_Type), POINTER :: FXR(:,:)
  REAL, POINTER :: PHIS(:,:,:)
  logical, intent(in) :: lpin_resolved
  integer, intent(in) :: ng
  logical, intent (inout) :: is_reset
  logical :: lxslib
  integer :: nrf, nxypin, nz, myzb, myze

  is_reset = .false.
  
  nrf = ThVar%npr
  nxypin = self%nxy
  nz = self%nz
  myzb = PE%myzb; myze = PE%myze
  lxslib = nTracerCntl%lXslib
  FXR => FmInfo%fxr
  PHIS => FmInfo%PHIS
  call BurnupUpdate
  is_reset = LBUUPDATED
  if (lpin_resolved) then
    call PowerUpdate_Watt(Core, Fxr, phis, Power_WATT, myzb, myze, ng, lxslib, GroupInfo, PE, .false.)
    call UpdatePinResolvedPower(Core,nxypin,myzb, myze, nrf, nz, PE%MPI_CMFD_COMM)
    if (LBUUPDATED) call UpdatePinResolvedBurnup(Core,nxypin,myzb, myze, nrf, nz, PE%MPI_CMFD_COMM)
  else
    if (LBUUPDATED) call UpdateAvgBurnup(Core,nxypin,myzb, myze, nz, PE%MPI_CMFD_COMM)
  end if

  NULLIFY(FXR, PHIS)
  
end subroutine

subroutine BurnupUpdate
  use DEPL_MOD,         only: DeplCntl
  implicit none
  integer, save :: bu_step_old = -1
  integer :: bu_step_new, PC_OPT
  logical :: LCORRECT
  logical :: LCORRECT_SAVE = .false.
  
  bu_step_new = DeplCntl%NowStep
  PC_OPT = DeplCntl%PC_OPT
  
  if (bu_step_new > bu_step_old) then
    if (PC_OPT == 1) then
      LBUUPDATED = .true.
    else
      LBUUPDATED = .true.
      LCORRECT_SAVE = .false.
    end if
  else
    if (PC_OPT == 1) then
      LBUUPDATED = .false.
    else
      LCORRECT = .not. DeplCntl%LPREDICT
      if (LCORRECT_save /= LCORRECT) then
        LBUUPDATED = .true.
      else
        LBUUPDATED = .false.
      end if 
      LCORRECT_SAVE = LCORRECT
    endif
  end if
  
  bu_step_old = bu_step_new

end subroutine

subroutine CalcGadFrac(Core,npin,myzb, myze, COMM)
  USE Typedef,         only: CoreInfo_Type, Cell_Type, Pin_Type,FxrInfo_Type
  USE Core_mod,        only: FmInfo
  use Material_Mod,    only: Mixture
  USE NuclidMap_mod,   only: AtomicWeight
#ifdef HAVE_MPI
  use MPIComm_Mod,     only: ALL_REDUCE_REAL
#endif
  implicit none
  TYPE(CoreInfo_Type) :: Core
  TYPE(Cell_Type), POINTER :: CellInfo(:)
  TYPE(Pin_Type), POINTER :: Pin(:)
  TYPE(FxrInfo_Type), POINTER :: FXR(:,:)
  integer, intent(in) :: npin,myzb,myze, COMM
  real*8 :: GD(MXNXSR), AW, gad_frac, VXSR(MXNXSR), HZ
  integer :: k, ixy, icel, nLocalFxr, j, FxrIdxSt, MM, LX, IX2, IM, NISO, i, ID
  logical :: lgadpin
#ifdef HAVE_MPI
  include 'mpif.h'
#endif  
  CellInfo => Core%CellInfo
  Pin => Core%Pin
  FXR => FmInfo%fxr
  GD = 0._8
  VXSR = 0._8

  do j=1,nz
    do i=1,npin
      gad_fracn(i,j) = 0._8
    end do
  end do
  
  do k = myzb,myze
    if (.NOT. Core%lFuelPlane(k)) cycle
    HZ = Core%hz(k)
    do ixy = 1, npin
      icel = Pin(ixy)%Cell(k)
      nLocalFxr = CellInfo(icel)%nFxr
      FxrIdxSt = Pin(ixy)%FxrIdxSt
      lgadpin = .false.
      gad_frac = 0._8
      MM = 0
      if(CellInfo(icel)%lfuel) then
        do j = 1,nLocalFxr
          LX = FxrIdxSt + j - 1
          IM = FXR(LX,K)%IMIX
          if (.not. Mixture(IM)%LFUEL) cycle
          MM = MM + 1
          IX2 = nLocalFxr - j + 1
          if (Mixture(IM)%LGD) then
            lgadpin = .TRUE.
            gad_frac = 0._8
            NISO = FXR(LX,k)%ndim
            do i = 1, NISO
              ID = FXR(LX,k)%idiso(i)
              if (ID<=64000 .or. ID>=64163) cycle
              AW = AtomicWeight(ID)
              gad_frac = gad_frac+Mixture(IM)%pnum(i)*AW
            end do
            VXSR(IX2) = FXR(LX,k)%area*HZ
            GD(IX2) = gad_frac/AVOGADRO/Mixture(IM)%dens
          end if
        end do
      end if
      if (MM > 0 .and. lgadpin) then
        gad_frac = 0._8
        do i = 1, MM
          gad_frac = gad_frac + GD(i)*VXSR(i)
        end do
        gad_fracn(IXY,K) = gad_frac/sum(VXSR(1:MM))
        VXSR(:) = 0. 
        GD(:) = 0.
      end if
    end do
  end do
#ifdef HAVE_MPI    
  call ALL_REDUCE_REAL(gad_fracn, npin*nz, COMM)
#endif
end subroutine

subroutine UpdatePinResolvedPower(Core,npin,myzb, myze, nrf, nz, COMM)
  USE Typedef,         only: CoreInfo_Type, Cell_Type, Pin_Type,FxrInfo_Type
  USE Core_mod,        only: FmInfo
  use Material_Mod,    only: Mixture
#ifdef HAVE_MPI
  use MPIComm_Mod,     only: ALL_REDUCE_REAL
#endif
  implicit none
  TYPE(CoreInfo_Type) :: Core
  TYPE(Cell_Type), POINTER :: CellInfo(:)
  TYPE(Pin_Type), POINTER :: Pin(:)
  TYPE(FxrInfo_Type), POINTER :: FXR(:,:)
  integer, intent(in) :: npin,myzb,myze,nrf,nz,COMM
  real*8 :: PXSR(MXNXSR), VOL(MXNXSR), HZ
  real*8, pointer :: FSR_PW(:,:)
  integer :: j, i, k, ixy, icel, nLocalFxr, FxrIdxSt, IFUEL, LX, IM, IX2
  integer :: IREG1, IREGFINAL, IR, IX, IFSR
#ifdef HAVE_MPI
  include 'mpif.h'
#endif  
  CellInfo => Core%CellInfo
  Pin => Core%Pin
  FXR => FmInfo%fxr
  
  PXSR(1:MXNXSR) = 0._8
  VOL(1:MXNXSR) = 0._8
  do j=1,nz
    do i=1,npin
      do k=1,nrf
        qvoln(k,i,j) = 0._8
      end do
    end do
  end do
  
  do k = myzb, myze
    if (.NOT. Core%lFuelPlane(k)) cycle
    HZ = Core%hz(k)
    do ixy = 1, npin
      icel = Pin(ixy)%Cell(k)
      nLocalFxr = CellInfo(icel)%nFxr
      FxrIdxSt = Pin(ixy)%FxrIdxSt
      if(CellInfo(icel)%lfuel) then
        PXSR(1:MXNXSR) = 0.
        VOL(1:MXNXSR) = 0.
        IFUEL = -99
        do j = 1,nLocalFxr
          LX = FxrIdxSt + j - 1
          IM = FXR(LX,K)%IMIX
          if (.not. Mixture(IM)%LFUEL) cycle
          IX2 = nLocalFxr - j + 1
          if (Mixture(IM)%LFUEL) then
            IFUEL = MAX(IX2,IFUEL)
            VOL(IX2) = VOL(IX2) + FXR(LX,k)%area*HZ
            IREG1 = FXR(LX,k)%FsrIdxSt
            IREGFINAL = IREG1 + FXR(LX,k)%nFsrInFxr - 1
            do IR = IREG1, IREGFINAL
              PXSR(IX2) = PXSR(IX2) + Power_WATT(IR,k)
            end do
          end if
        end do
        do IX = 1, IFUEL
          if(VOL(IX) /=0.) PXSR(IX) = PXSR(IX)/VOL(IX)*1.D+06
        end do
        call transfer_info(icel, nLocalFxr, PXSR(1:IFUEL), IFUEL, qvoln(1:nrf,IXY,K), nrf, 1,&
                            CellInfo(icel))
      end if
    end do
  end do
#ifdef HAVE_MPI    
  call ALL_REDUCE_REAL(qvoln, nrf*npin*nz, COMM)
#endif
  NULLIFY(FXR, PIN, CELLINFO)
end subroutine

subroutine UpdatePinResolvedBurnup(Core,npin,myzb, myze, nrf, nz, COMM)
  USE Typedef,         only: CoreInfo_Type, Cell_Type, Pin_Type,FxrInfo_Type
  USE Core_mod,        only: FmInfo
  use Material_Mod,    only: Mixture
#ifdef HAVE_MPI
  use MPIComm_Mod,     only: ALL_REDUCE_REAL
#endif
  implicit none
  TYPE(CoreInfo_Type) :: Core
  TYPE(Cell_Type), POINTER :: CellInfo(:)
  TYPE(Pin_Type), POINTER :: Pin(:)
  TYPE(FxrInfo_Type), POINTER :: FXR(:,:)
  integer, intent(in) :: npin,myzb,myze,nrf,nz, COMM
  real*8 :: BUXSR(MXNXSR), VOL(MXNXSR), HZ
  integer :: j, i, k, ixy, icel, nLocalFxr, FxrIdxSt, IFUEL, LX, IM, IX2, irf
  real*8 :: vol_tot, HMDEN0
#ifdef HAVE_MPI
  include 'mpif.h'
#endif  
  CellInfo => Core%CellInfo
  Pin => Core%Pin
  FXR => FmInfo%fxr
  
  BUXSR(1:MXNXSR) = 0._8
  VOL(1:MXNXSR) = 0._8
  vol_tot = 0._8
  HMDEN0 = 0._8
  do j=1,nz
    do i=1,npin
      do k=1,nrf
        burnupn(k,i,j) = 0._8
      end do
    end do
  end do
  
  do k = myzb, myze
    if (.NOT. Core%lFuelPlane(k)) cycle
    HZ = Core%hz(k)
    do ixy = 1, npin
      icel = Pin(ixy)%Cell(k)
      nLocalFxr = CellInfo(icel)%nFxr
      FxrIdxSt = Pin(ixy)%FxrIdxSt
      if(CellInfo(icel)%lfuel) then
        BUXSR(1:nLocalFxr) = 0.
        VOL(1:nLocalFxr) = 0.
        vol_tot = 0._8
        IFUEL = -99
        do j = 1,nLocalFxr
          LX = FxrIdxSt + j - 1
          IM = FXR(LX,K)%IMIX
          if (.not. Mixture(IM)%LFUEL) cycle
          IX2 = nLocalFxr - j + 1
          if (Mixture(IM)%LFUEL) then
            IFUEL = MAX(IX2,IFUEL)
            VOL(IX2) = FXR(LX,k)%area*HZ
            if (VOL(IX2) /= 0.) HMDEN0 = FXR(LX,K)%Hmkg0/VOL(IX2)
            BUXSR(IX2) = FXR(LX,k)%Burnup*HMDEN0
            vol_tot = vol_tot + VOL(IX2)
          end if
        end do
        call transfer_info(icel, nLocalFxr, BUXSR(1:IFUEL), IFUEL, burnupn(1:nrf,IXY,K), nrf, 2,&
                            CellInfo(icel))
		do irf=1,nrf					
          if (ihmmassn(irf,IXY,K) /= 0.) burnupn(irf,IXY,K)=burnupn(irf,IXY,K)/ihmmassn(irf,IXY,K)*vol_tot
		end do
      end if
    end do
  end do  
#ifdef HAVE_MPI    
  call ALL_REDUCE_REAL(burnupn, nrf*npin*nz, COMM)
#endif  
  NULLIFY(FXR, PIN, CELLINFO)
end subroutine

subroutine UpdateAvgBurnup(Core,npin,myzb, myze, nz, COMM)
  USE Typedef,         only: CoreInfo_Type, Cell_Type, Pin_Type,FxrInfo_Type
  USE Core_mod,        only: FmInfo
  use Material_Mod,    only: Mixture
#ifdef HAVE_MPI
  use MPIComm_Mod,     only: ALL_REDUCE_REAL
#endif
  implicit none
  TYPE(CoreInfo_Type) :: Core
  TYPE(Cell_Type), POINTER :: CellInfo(:)
  TYPE(Pin_Type), POINTER :: Pin(:)
  TYPE(FxrInfo_Type), POINTER :: FXR(:,:)
  integer, intent(in) :: npin,myzb,myze,nz, COMM
  integer :: j, i, k, ixy, icel, nLocalFxr, FxrIdxSt, IFUEL, LX, IM, IX2
  real*8 :: HMTOT0, BU_PIN
#ifdef HAVE_MPI
  include 'mpif.h'
#endif  
  CellInfo => Core%CellInfo
  Pin => Core%Pin
  FXR => FmInfo%fxr
  
  BU_PIN = 0._8
  HMTOT0 = 0._8
  do j=1,nz
    do i=1,npin
      burnup_avg(i,j) = 0._8
    end do
  end do
  
  do k = myzb, myze
    if (.NOT. Core%lFuelPlane(k)) cycle
    do ixy = 1, npin
      icel = Pin(ixy)%Cell(k)
      nLocalFxr = CellInfo(icel)%nFxr
      FxrIdxSt = Pin(ixy)%FxrIdxSt
      if(CellInfo(icel)%lfuel) then
        BU_PIN = 0.
        HMTOT0 = 0.
        IFUEL = -99
        do j = 1,nLocalFxr
          LX = FxrIdxSt + j - 1
          IM = FXR(LX,K)%IMIX
          if (.not. Mixture(IM)%LFUEL) cycle
          IX2 = nLocalFxr - j + 1
          if (Mixture(IM)%LFUEL) then
            IFUEL = MAX(IX2,IFUEL)
            HMTOT0 = HMTOT0 + FXR(LX,K)%Hmkg0
            BU_PIN = BU_PIN + FXR(LX,k)%Burnup*FXR(LX,K)%Hmkg0
          end if
        end do
        BU_PIN = BU_PIN/HMTOT0
        burnup_avg(ixy,k) = BU_PIN
      end if
    end do
  end do  
#ifdef HAVE_MPI    
  call ALL_REDUCE_REAL(burnup_avg, npin*nz, COMM)
#endif
  NULLIFY(FXR, PIN, CELLINFO)
end subroutine

subroutine transfer_info(icell, NXSR_loc, NEUT_data,ifuel,TH_data,iring_th,imode,CellInfo)
  use TYPEDEF,     only: Cell_Type, THCell_Type
  implicit none
  TYPE(Cell_Type) :: CellInfo
  TYPE(THCell_Type), POINTER :: THCell
  integer, intent (in) :: icell, ifuel, iring_th, NXSR_loc, imode
  real*8, intent(in) :: NEUT_data(1:ifuel)     !quantity / volume
  real*8, intent(inout) :: TH_data(1:iring_th) !imode = 1, quantity/volume, imode = 2, quantity/voltot
  integer :: ix, j, ix2
  real*8 :: sumfrac
  
  ThCell => CellInfo%ThCell
  
  do j = 1,iring_th
    TH_data(j) = 0.
    sumfrac = 0.
    do ix = 1,IFUEL
      IX2 = NXSR_loc-IX+1
      TH_data(j) = TH_data(j) + Thcell%Invfrac(ix2,j)*NEUT_data(ix)
      sumfrac = sumfrac + Thcell%Invfrac(ix2,j)
    enddo
    if (imode == 1) then
      if(sumfrac /= 0.) TH_data(j) = TH_data(j)/sumfrac
    else
      continue
    endif
  enddo
end subroutine transfer_info
  
  
end module