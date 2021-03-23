module SubChCoupling_mod
use PARAM
use TYPEDEF,      only: CoreInfo_type,    ThInfo_Type,    Asy_type,       &
                        AsyInfo_type, PE_type
implicit none
private
public :: is_coupled,             &
          CodeName,               &
          SubChCode,              &
          ActiveAsymIndex,        &
          last_TH,                &
          coupled_maxsteps,       &
          coupled_relconv,        &
          Courant,                &
          sbch_outop

logical :: is_coupled=.false., last_TH=.false.
character (len=10) :: CodeName
integer :: coupled_maxsteps=-1
real*8 :: coupled_relconv= 0.02, Courant = 0.6
integer :: nxy, nz, &
           sbch_outop(2) = [0, 0]      !1: Txt output 2: 3D output
real*8 :: PowLin, PowLv

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
  
end type SubChCode

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!preprocessing: Index converter, MPIGrouping
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef __PGI

SUBROUTINE Preproc(self, Core, PE_nT)
IMPLICIT NONE

CLASS(SubChCode) :: self
TYPE(CoreInfo_TYPE), INTENT(IN) :: Core
TYPE(PE_TYPE), INTENT(IN) :: PE_nT

INTEGER :: iz, iasy, iasytype, id
INTEGER :: ierr
INTEGER, POINTER :: proc_list(:)
#ifdef HAVE_MPI
INCLUDe 'mpif.h'
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

CALL SetActCoreReg(self, Core)

CALL IndexConversion(self, Core)

IF (self%run_parallel) THEN
  SELECT CASE(self%CodeName)
  CASE("ESCOT")
    self%PE%nproc = PE_nT%nproc
  CASE("CTF")
    self%PE%nproc = self%nasy_act
    IF (self%PE%nproc .LE. 1) self%run_parallel = .FALSE.
  END SELECT
#ifdef HAVE_ESCOT
  CALL MPI_Comm_Group(PE_nT%MPI_CMFD_COMM, self%PE%WorldGroup, ierr)
#endif
  IF (PE_nT%myrank .LE. self%PE%nProc) THEN
    ALLOCATE(proc_list(self%PE%nproc))
    DO id = 0, self%PE%nproc - 1
      proc_list(id + 1) = id
    END DO
#ifdef HVAE_MPI
    CALL MPI_Group_incl(self%PE%WorldGroup, self%PE%nproc, proc_list, &
                        self%PE%LocalGroup, ierr)
    self%PE%MPI_COMM = PE_nT%MPI_CMFD_COMM
#endif
    self%PE%lSubCh_proc = .TRUE.
    DEALLOCATE(proc_list)
  END IF
ELSE
  self%PE%nproc = 1
  self%PE%lSubCh_proc = .TRUE.
END IF

END SUBROUTINE

#else
subroutine preproc(self, Core, PE_nT)
  class (SubChCode) :: self
  type(CoreInfo_type), intent (in) :: Core
  type(PE_type), intent (in) :: PE_nT
  integer :: iz, iasy, iasytype, id
  integer :: ierr
  integer, allocatable :: proc_list(:)
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
  
  call IndexConversion(self, Core)

  if (self%run_parallel) then
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

#endif
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

!  do iz=1,self%nz
!    do ixy = 1, self%nxy
!      id_SubCh=self%id_nT_to_SubCh(ixy)
!      if(id_SubCh == 0) cycle
!      qprime=ThInfo%RelPower(iz, ixy)*ThInfo%PowLin*ThInfo%Powlv/100._8   !W/m -> W/cm
!      call set_rodpowers_W_cm(id_SubCh,iz,qprime)
!    enddo
!  enddo

  !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(iz, ixy, id_SubCh, qprime)
      do ixy = 1, self%nxy
        id_SubCh=self%id_nT_to_SubCh(ixy)
    IF (id_SubCh .EQ. 0) CYCLE
    DO iz = 1, self%nz
        qprime=ThInfo%RelPower(iz, ixy)*ThInfo%PowLin*ThInfo%Powlv/100._8   !W/m -> W/cm
        call set_rodpowers_W_cm(id_SubCh,iz,qprime)
      enddo    
    enddo
  !$OMP END PARALLEL DO

end subroutine SetSubChData

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

!  do iz=1,self%nz
!    do ixy = 1, self%nxy
!      id_SubCh=self%id_nT_to_SubCh(ixy)
!      if(id_SubCh == 0) cycle
!      iz_subch=iz+idum
!      call get_coolant_temp(id_SubCh,iz_subch,ThInfo%Tcool(iz,ixy))
!      call get_coolant_dens(id_SubCh,iz_subch,ThInfo%DenCool(iz,ixy))
!    enddo
!  enddo

  !$OMP PARALLEL
  !$OMP DO SCHEDULE(GUIDED) PRIVATE(iz, ixy, id_SubCh)
  DO ixy = 1, self%nxy
    id_SubCh = self%id_nT_to_SubCh(ixy)
    IF (id_SubCh .EQ. 0) CYCLE
  do iz=1,self%nz
      CALL get_coolant_temp(id_SubCh, iz, ThInfo%Tcool(iz,ixy))
    END DO
  END DO
  !$OMP END DO

  !$OMP DO SCHEDULE(GUIDED) PRIVATE(iz, ixy, id_SubCh)
      do ixy = 1, self%nxy
        id_SubCh=self%id_nT_to_SubCh(ixy)
    IF (id_SubCh .EQ. 0) CYCLE
    DO iz = 1, self%nz
      CALL get_coolant_dens(id_SubCh, iz, ThInfo%DenCool(iz,ixy))
      enddo    
    enddo
  !$OMP END DO
  !$OMP END PARALLEL

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
	logical, intent(in) :: is_ftemp_comm
  integer :: comm, i, j, npr2

  comm = PE%MPI_CMFD_COMM

  call BCAST(ThInfo%Tcool, self%nz, self%nxy+2, COMM)
  call BCAST(ThInfo%Dencool, self%nz, self%nxy+2, COMM)
  
	if (.not. is_ftemp_comm) return

	npr2 = size(Thinfo%Tfvol,1)
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
#ifndef __PGI
  call MPI_Group_Free(self%PE%LocalGroup,ierr)
  call MPI_Group_Free(self%PE%WorldGroup,ierr)
  call MPI_Comm_Free(self%PE%MPI_COMM,ierr)
#endif
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
  
  allocate(self%core_conf(self%nxa_act%n, self%nya_act%n))
  
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


end module