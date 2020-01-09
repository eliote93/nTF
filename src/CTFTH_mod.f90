#include <defines.h>

module CTFTH_mod
use PARAM
use typedef, ONLY : CoreInfo_Type,      Pin_Type,       PinInfo_TYPE,   &
                    ThInfo_Type,        Asy_Type,       AsyInfo_Type,   &   
                    PE_TYPE,            FuelTh_Type     
use CNTL,    ONLY : nTRACERCntl_Type
use TH_MOD,       ONLY : ThVar
use SubChCoupling_mod,      only: last_TH, SubChCode, coupled_maxsteps

implicit none
private
public :: CTF_TH


logical :: lfirst = .true.

type, extends(SubChCode) :: SubChCode_modified
  contains
  procedure :: SetSubChData_CTF

end type SubChCode_modified

type(SubChCode_modified) :: CTF

contains

subroutine SetSubChData_CTF(self,Core,ThInfo,set_rodpowers_W_cm)
  class(SubChCode_modified) :: self
  type(CoreInfo_type), intent (in) :: Core
  type(ThInfo_Type), intent (in) :: ThInfo
  interface 
  subroutine set_rodpowers_W_cm(pin_number,axial_level,rodpower)
    integer,intent(in) :: pin_number, axial_level
    real*8,intent(in) :: rodpower
  end subroutine
  end interface
  integer :: iz, iasy, iasytype, ix, iy, id_Subch, ixy, iz_subCh
  real*8 :: qprime
  do iz=1,self%nz
    do iasy=1,self%nxya
      iasytype = Core%CoreMap(iasy)
      if (.not.core%AsyInfo(iasytype)%lFuel) cycle
      do ixy = 1, self%nxy
        id_SubCh=self%id_nT_to_SubCh(ixy)
        if(id_SubCh == 0) cycle
        qprime=ThInfo%RelPower(iz, ixy)*ThInfo%PowLin*ThInfo%Powlv/100._8   !W/m -> W/cm
        iz_subch=iz+1
        call set_rodpowers_W_cm(id_SubCh,iz_subch,qprime)
      enddo    
    enddo
  enddo

  do iasy=1,self%nxya
    iasytype = Core%CoreMap(iasy)
    if (.not.Core%AsyInfo(iasytype)%lFuel) cycle
    do ixy = 1, self%nxy
      id_SubCh=self%id_nT_to_SubCh(ixy)
      if(id_SubCh == 0) cycle
      qprime=ThInfo%RelPower(1, ixy)*ThInfo%PowLin*ThInfo%Powlv/100._8   !W/m -> W/cm
      call set_rodpowers_W_cm(id_SubCh,1,qprime)
      qprime=ThInfo%RelPower(self%nz, ixy)*ThInfo%PowLin*ThInfo%Powlv/100._8   !W/m -> W/cm
      call set_rodpowers_W_cm(id_SubCh,self%nz+2,qprime)
    enddo    
  enddo
end subroutine SetSubChData_CTF




subroutine CTF_TH(Core, ThInfo, nTracerCntl,ToutAvg, PE)
  use MPIComm_Mod,      only: MPI_SYNC
  type(CoreInfo_Type) :: Core
  type(ThInfo_Type) :: ThInfo
  type(nTracerCntl_Type) :: nTracerCntl
  type(PE_Type) :: PE
  real*8 :: ToutAvg  
  interface
  subroutine CTF_set_rodpowers_W_cm(pin_number,axial_level,rodpower)
    integer, intent(in) :: pin_number, axial_level
    real*8, intent(in) :: rodpower
  end subroutine
  subroutine CTF_get_coolant_temp(pin_number,axial_level,temp)
    integer, intent (in) :: pin_number, axial_level
    real*8, intent (inout) :: temp
  end subroutine
  subroutine CTF_get_coolant_dens(pin_number,axial_level,density)
    integer, intent (in) :: pin_number, axial_level
    real*8, intent (inout) :: density
  end subroutine
  subroutine CTF_set_basefn_F(input_basefn)
    character(len=*) :: input_basefn
  end subroutine CTF_set_basefn_F
  subroutine CTF_Initialize(cobra_subcomm,hdf5_has_been_initialized,&
                          model_CRUD)
   integer, intent(in) :: cobra_subcomm
   logical, intent(in) :: hdf5_has_been_initialized
   logical, optional, intent(in) :: model_CRUD
  end subroutine
  subroutine set_coupled_maxsteps(maxsteps)
   integer, intent(in) :: maxsteps
  end subroutine set_coupled_maxsteps
  end interface

  if (lfirst) then
    call CTF%preproc(Core,PE)
#ifdef HAVE_CTF
    if (PE%master) then
      call gen_preproc_input_CTF(Core,nTracerCntl)
      call CTF_preproc()      
    endif
    call MPI_SYNC(CTF%PE%MPI_COMM)

    if(CTF%PE%lSubCh_proc) then

      call CTF_set_basefn_F('deck')
      CALL CTF_Initialize(CTF%PE%MPI_COMM,.FALSE.)
      if (coupled_maxsteps>0) call set_coupled_maxsteps(coupled_maxsteps)
    endif
#endif    
    lfirst=.false.
  endif
#ifdef HAVE_CTF
  if (CTF%PE%lSubCh_proc) then
    call CTF%SetSubChData_CTF(Core,ThInfo,CTF_set_rodpowers_W_cm)

    call CTF_Solve_Standalone

    call MPI_SYNC(CTF%PE%MPI_COMM)  

    call CTF%GetSubChData(Core,ThInfo,CTF_get_coolant_temp,CTF_get_coolant_dens,    &
                          axial_idx_calib=1)

    call CTF_get_coolant_temp_exit(ToutAvg)
  endif  

#endif  

  call CTF%Bcast_SubChOut(ThInfo,PE)
  
  if (Core%lgap) then
    call CTF%SetAsyGapCoolinfo(Core,ThInfo)
  endif
  
  if (last_TH) then
#ifdef HAVE_CTF
    if (CTF%PE%lSubCh_proc) then
      CALL CTF_Edits(1)
      CALL CTF_Cleanup()
    endif
#endif      
    call CTF%Finalize()
  endif

  call MPI_SYNC(PE%MPI_CMFD_COMM)  

end subroutine



subroutine gen_preproc_input_CTF(Core,nTracerCntl)
  use files,              only: caseid
  use geom,               only: NCELLX0, lCoreAng, lEdge, lRot
  use ioutil,             only: newunit, terminate, toupper, nfields
  use TH_mod,             only: Thvar
  use SubChCoupling_mod,  only: ActiveAsymIndex
  use files,              only: io5
  type(CoreInfo_Type) :: Core
  type(nTracerCntl_Type) :: nTracerCntl
  type(AsyInfo_Type), pointer :: AsyInfo(:)
  type(Asy_Type), pointer :: Asy(:)
  type(Pin_Type), pointer :: Pin(:)

  INTEGER :: iz, ixa, iya, iasytype, nxa_act_WholeCore, nya_act_WholeCore, nasy_act_WholeCore,    &
              nactasy_col_centline, nactasy_row_centline
  INTEGER :: indev, ctfinp1, ctfinp2, ctfinp3, ctfinp4
  INTEGER, allocatable :: WholeCore_Conf(:,:)
  logical :: lgt, uniformz, flagspc=.false.
  INTEGER :: ib, ixb,ixe,iyb,iye, k, i, j,    &
              ixp, iyp, ixy_loc, ixy_glob, ixya, m, ref_assem_id, ref_assem_type
  character :: oneline*512, card*10, probe
  character(len=30) :: casename
  integer :: nfield, dum_int
  integer,dimension(50) :: idum=0
  real*8 :: hzact_start, hzact_tot, hzact_end
  !=================================================================================!
  !========CONTROL CARD=========!
  !Main Control data
  INTEGER ::  MAPS, ICOBRA, OITMAX, IITMAX, IRFC, EDMOD
  REAL*8 :: EPSO, COURANT
  CHARACTER*3 :: opt_par
  !OTHER CONTROL
  INTEGER :: ISOL
  !Mixing and void drift model
  INTEGER :: IMIXV
  REAL*8 :: AAAK, BETA, THETM
  !BOUNDARY CONTITION & INITIAL CONDITION
  REAL*8 :: GINIT, TINIT, HGIN, DHFRAC, IBCVAL1, IBCVAL2, IBCVAL3, OBCVAL2, OBCVAL3, PREF
  INTEGER :: ISPEC, OSPEC
  !TIME DOMAIN DATA
  REAL*8 :: DTMIN, DTMAX, TEND, RTWFP
  INTEGER :: MAXITS
  !Boron control
  INTEGER :: ntypeboron, nasyrow, nasycol, iBTM
  REAL*8 :: BRIN, RDIF
  integer, allocatable :: iasy_br(:,:)
  REAL*8, allocatable :: asy_br_con(:)
  !CONVERGENCE DATA
  REAL*8 :: GEB, GMB,FES, SES, MSO
  !OUTPUT CONTROL
  INTEGER :: OCCHAN, OCGAPS, OCRODS, OCDNB, OCVTK, OChdf
  REAL*8 :: boron_con
  !========GEOMETRY CARD==========!
  INTEGER :: OPT_SYMM
  !AXIAL MESH INFORAMTION
  INTEGER :: NONODE, IVARDX, nspc
  INTEGER, allocatable :: JLEV(:)
  REAL*8, allocatable :: VARDX(:), spcz(:), CDL(:)
  REAL*8 :: dxs
  !=========ASSEMBLY CARD=========
  INTEGER :: Casing,HTYPE, cmf
  REAL*8 :: FTDENS,HGAP
  type gtinfo_type
    integer :: x=0, y=0
  endtype
  type(gtinfo_type), allocatable :: gt_loc(:)


  indev=io5

  AsyInfo=>Core%AsyInfo
  Asy => Core%Asy
  Pin => Core%Pin


if(lCoreAng.eq.90 .and. .not.lEdge) THEN
  if (lRot) THEN
    OPT_Symm=5
  else
    OPT_Symm=4
  endif
  nactasy_col_centline=0
  do i=1,CTF%nxa_act%n
    if(CTF%core_conf(i,1)==0) cycle
    nactasy_col_centline=nactasy_col_centline+1
  enddo
  nactasy_row_centline=0
  do j=1,CTF%nya_act%n
    if(CTF%core_conf(1,j)==0) cycle
    nactasy_row_centline=nactasy_row_centline+1
  enddo
  if (CTF%Core_conf(1,1)==0) then
    dum_int=0
  else
    dum_int=1
  endif

  nasy_act_WholeCore= CTF%nasy_act*4-nactasy_col_centline*2       &
                        -nactasy_row_centline*2+dum_int
  nxa_act_WholeCore=2*CTF%nxa_act%n-1
  nya_act_WholeCore=2*CTF%nya_act%n-1
elseif (lCoreang.eq.360) THEN
  OPT_Symm=1
  nasy_act_WholeCore=CTF%nasy_act
  nxa_act_WholeCore=CTF%nxa_act%n
  nya_act_WholeCore=CTF%nya_act%n
elseif (lCoreang.eq.90 .and. lEdge) THEN
  OPT_Symm=1
  nasy_act_WholeCore=CTF%nasy_act
  nxa_act_WholeCore=CTF%nxa_act%n
  nya_act_WholeCore=CTF%nya_act%n
else
  call terminate("nTRACER/CTF Coupling cannot handle this type of Core Radial Config. Check Input Card ''RAD_CONF'' ")
ENDIF




lgt=.false.
! Guide Tube Information
  if (Core%nAsyGT>0) lgt=.true.
  if (lgt) then
    allocate(gt_loc(Core%nAsyGT))
    sr_ref: do iya=1,CTF%nya
      do ixa=1,CTF%nxa
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

  hzact_tot=0.
  hzact_start=0.

  axial_sweep: do iz=1,CTF%nz
    if(.not. Core%lFuelPlane(iz)) then
      hzact_start=hzact_start+Core%hz(iz)*0.01  ![cm]-<[m]
      cycle
    endif
    exit axial_sweep
  enddo axial_sweep


  do iz=1,CTF%nz
    if(.not. Core%lFuelPlane(iz)) cycle
    hzact_tot=hzact_tot+Core%hz(iz)*0.01      ![cm]->[m]
  enddo
  hzact_end=hzact_start+hzact_tot
  
  ivardx=1
  do iz=2,CTF%nz
    if(Core%hz(iz)/=Core%hz(iz-1)) then
      IVARDX=IVARDX+1
    endif
  enddo
  allocate(JLEV(IVARDX),VARDX(IVARDX))
  ivardx=1
  VARDX=0.;
  JLEV(1)=1
  VARDX(1)=Core%hz(1)
  do iz=2,Core%nz
    if(Core%hz(iz)/=Core%hz(iz-1)) then
      ivardx=ivardx+1
      JLEV(ivardx)=1
      VARDX(ivardx)=VARDX(ivardx-1)+Core%hz(iz)
    else
      JLEV(ivardx)=JLEV(ivardx)+1
      VARDX(ivardx)=VARDX(ivardx)+Core%hz(iz)
    endif
  enddo
  VARDX=VARDX*10.       !Convert Unit [cm] -> [mm]




allocate(WholeCore_Conf(nxa_act_WholeCore,nya_act_WholeCore))
WholeCore_Conf=0
if (lCoreAng==90 .and. .not.lEdge) THEN
  ixb=CTF%nxa_act%n
  iyb=CTF%nya_act%n
  do iya=iyb,nya_act_WholeCore
    do ixa=ixb, nxa_act_WholeCore
      if(CTF%Core_conf(ixa-ixb+1,iya-iyb+1)/=0) WholeCore_Conf(ixa,iya)=1
    enddo
  enddo
  do iya=iyb,nya_act_WholeCore
    do ixa=ixb-1,1,-1
      WholeCore_Conf(ixa,iya)=WholeCore_Conf(nxa_act_WholeCore-ixa+1,iya)
    enddo
  enddo
  do iya=1,iyb-1
    do ixa=1,nxa_act_WholeCore
      WholeCore_Conf(ixa,iya)=WholeCore_Conf(ixa,nya_act_WholeCore-iya+1)
    enddo
  enddo
else
  do iya=1,CTF%nya_act%n
    do ixa=1,CTF%nxa_act%n
      if(CTF%Core_conf(ixa,iya)/=0) WholeCore_Conf(ixa,iya)=1
    enddo
  enddo
endif


!========INITIALIZATION AND SET DEFAULT VALUE========    
!======control.inp
  if (CTF%run_parallel) THEN
    opt_par='yes'
  ELSE
    opt_par='no '
  endif
  MAPS=1
  ICOBRA=1
  EPSO=0.001
  OITMAX=5
  IITMAX=40
  COURANT=0.800
  IRFC=1
  EDMOD=0
  IMIXV=1
  AAAK=1.400
  BETA=0.05000
  if (CTF%run_parallel) THEN
    ISOL=5
  ELSE
    ISOL=3
  endif
  ibtm=0
!INITIAL CONDITIONS
  GINIT=nTracerCntl%fmdotfa*nasy_act_WholeCore
  PREF=nTracerCntl%Pexit*1.E-5
  TINIT=-nTracerCntl%TempInlet
  HGIN=288.42000
  DHFRAC=0.000
  !GLOBAL BOUNDARY CONDITIONS
  ispec=2
  Ospec=1
  IBCVAL1=nTracerCntl%fmdotfa*nasy_act_WholeCore
  IBCVAL2=-nTracerCntl%TempInlet
  OBCVAL2=nTracerCntl%TempInlet
  IBCVAL3=0.000
  OBCVAL3=nTracerCntl%Pexit*1.E-5
  DTMIN=0.0000001
  DTMAX=0.1
  TEND=0.1
  RTWFP=1000.0
  MAXITS=10000
  GEB=0.01
  GMB=0.01
  FES=0.5
  SES=0.5
  MSO=0.5
  !OUTPUT CONTROL
  OCCHAN=1
  OCGAPS=1
  OCRODS=1
  OCDNB=1
  OCVTK=0
  OChdf=0


  !=============assembly.inp=========
  Casing=0
  HTYPE=0
  CMF=1
  FTDENS=95.0
  HGAP=5678.3   !W/(m^2*K)


! read nTRACER's input deck to take additional users' input ...
rewind(indev)
do while(.TRUE.)
  read(indev,'(a256)',end=1300) oneline
  read(oneline,'(a1)') probe
       
  if(probe .eq. '.') exit
  if(probe.eq.'!') cycle
  if(oneline.eq.'') cycle
  if(probe.eq.'/') exit
       
  read(oneline,*) card; call toupper(card)
  nfield=nfields(oneline)-1
  select case(card)
    case('CONTROL')
      do while(.TRUE.)
        read(indev,'(a256)',end=1300) oneline
        read(oneline,'(a1)') probe
        if(probe .eq. '.') exit
        if(probe.eq.'!') cycle
        if(oneline.eq.'') cycle
        if(probe.eq.'/') exit
        read(oneline,*) card; call toupper(card)
        if((card.NE.'CMAXOUTER') .and. (card.NE.'BORONTRACK')) exit
        select case(card)
          case('CMAXOUTER')
            read(oneline,*) card, coupled_maxsteps
          case('BORONTRACK')
            read(oneline,*) card, ibtm
            if (ibtm.EQ.1 ) THEN
              brin=nTracerCntl%BoronPPM; Rdif=0.0; nTracerCntl%lborontrack= .TRUE.
            elseif (ibtm.EQ.2) THEN
              brin=nTracerCntl%BoronPPM; Rdif=1.0; nTracerCntl%lborontrack= .TRUE.
            ENDIF
            read(indev,*) ntypeboron
            allocate(asy_br_con(ntypeboron))
            asy_br_con=0.0
            read(indev,*) (asy_br_con(i), i=1,ntypeboron)
            read(indev,*) nasyrow, nasycol
            allocate(iasy_br(nasyrow,nasycol))
            iasy_br=0
            do i=1, nasyrow
              read(indev,*) (iasy_br(i,j), j=1,nasycol)
            ENDDO
        END SELECT
      ENDDO
    case('CTFGEOM')
      do while(.TRUE.)
        read(indev,'(a256)',end=1300) oneline
        read(oneline,'(a1)') probe
        if(probe .eq. '.') exit
        if(probe.eq.'!') cycle
        if(oneline.eq.'') cycle
        if(probe.eq.'/') exit
        read(oneline,*) card; call toupper(card)
        if(card.NE.'SPACEGRID') exit
        select case(card)
          case('SPACEGRID')
            read(oneline,*) card, nspc
            if (nspc > 0) THEN
              flagspc=.TRUE.
              allocate(spcz(nspc),CDL(nspc))
              spcz=0.;CDL=0.9026
              do i=1,nspc
                read(indev,'(a256)') oneline
                nfield=nfields(oneline)
                if(nfield.eq.2) then
                  read(oneline,*) spcz(i), CDL(i)
                elseif (nfield.eq.1) then
                  read(oneline,*) spcz(i)
                else
                  call terminate("Wrong elements are input at SPACERGRID. Check input card") 
                endif
              ENDDO
            endif
        END SELECT
      ENDDO 
  end select
enddo
!     
1300 continue


  open(newunit(ctfinp1),file='control.inp', status='unknown')
  open(newunit(ctfinp2),file='geo.inp', status='unknown')
  open(newunit(ctfinp3),file='assembly.inp',status='unknown')
  open(newunit(ctfinp4),file='power.inp',status='unknown')
  write(ctfinp1,'(a15)') caseid
  write(ctfinp1,'(a10)') '{parallel}'
  write(ctfinp1,'(a3)') opt_par
  write(ctfinp1,'(i3)') MAPS
  if (MAPS.EQ.1) THEN
    write(ctfinp1,'(a15)') caseid
    write(ctfinp1,'(a15)') caseid
  ENDIF
  write(ctfinp1,'(i3)') ICOBRA
  write(ctfinp1,'(f10.4)') EPSO
  write(ctfinp1,'(i3)') OITMAX
  write(ctfinp1,'(i3)') IITMAX
  write(ctfinp1,'(f10.4)') COURANT
  write(ctfinp1,'(i3)') IRFC
  write(ctfinp1,'(i3)') EDMOD
  write(ctfinp1,'(i3)') IMIXV
  IF (ANY(IMIXV.EQ.(/1,2,3/))) write(ctfinp1,'(f10.4)') AAAK
  IF(ANY(IMIXV.EQ.(/1,3/))) write(ctfinp1,'(f10.4)')BETA
  IF(ANY(IMIXV.EQ.(/2,3/))) write(ctfinp1,'(f10.4)')THETM
  write(ctfinp1,'(i3)')ISOL
  write(ctfinp1,'(f10.4)')GINIT
  write(ctfinp1,'(f10.4)')-TINIT
  write(ctfinp1,'(f10.4)')PREF
  write(ctfinp1,'(f10.4)')TINIT
  write(ctfinp1,'(f10.4)')HGIN
  write(ctfinp1,'(f10.4)')DHFRAC
  write(ctfinp1,'(i3)') ISPEC
  write(ctfinp1,'(i3)') OSPEC
  write(ctfinp1,'(f10.4)') IBCVAL1
  write(ctfinp1,'(f10.4)') IBCVAL2
  write(ctfinp1,'(f10.4)') OBCVAL2
  write(ctfinp1,'(f10.4)') IBCVAL3
  write(ctfinp1,'(f10.4)') OBCVAL3
  write(ctfinp1,'(i3)') ibtm
  if(ibtm>0) THEN
    write(ctfinp1,'(2f10.4)') brin, Rdif
    write(ctfinp1,'(i3)') ntypeboron
    write(ctfinp1,'(f10.2)') (asy_br_con(i), i=1,ntypeboron)
    write(ctfinp1,'(2i5)') nasyrow, nasycol
    do i=1, nasyrow
      write(ctfinp1,'(i3)') (iasy_br(i,j), j=1,nasycol)
    ENDDO
  ENDIF
  write(ctfinp1,'(ES11.4)') DTMIN
  write(ctfinp1,'(f10.4)') DTMAX
  write(ctfinp1,'(f10.4)') TEND
  write(ctfinp1,'(f10.4)') RTWFP
  write(ctfinp1,'(i10)') MAXITS
  write(ctfinp1,'(a22)') '{convergence criteria}'
  write(ctfinp1,'(f10.4,/,f10.4,/f10.4,/f10.4,/f10.4)') GEB, GMB,FES, SES, MSO
  write(ctfinp1,'(a15)') '{edit channels}'
  write(ctfinp1,'(i3)') OCCHAN
  write(ctfinp1,'(a11)') '{edit gaps}'
  write(ctfinp1,'(i3)') OCGAPS
  write(ctfinp1,'(a11)') '{edit rods}'
  write(ctfinp1,'(i3)') OCRODS
  write(ctfinp1,'(a10)') '{edit dnb}'
  write(ctfinp1,'(i3)') OCDNB
  write(ctfinp1,'(a10)') '{rods vtk}'
  write(ctfinp1,'(i3)') OCVTK
  write(ctfinp1,'(a11)') '{edit hdf5}'
  write(ctfinp1,'(i3)') OChdf

  !Goemotry input
  write(ctfinp2,'(i5)') nasy_act_WholeCore
  write(ctfinp2,'(i5)') 1
  write(ctfinp2,'(2i5)') nxa_act_WholeCore, nya_act_WholeCore
  write(ctfinp2,'(a17)') '{symmetry option}'
  write(ctfinp2,'(i5)') OPT_SYMM
  write(ctfinp2,'(a2)',ADVANCE='no') '** '
  DO i=1,nxa_act_WholeCore
    write(ctfinp2,'(i3)',ADVANCE='no') i
  ENDDO
  write(ctfinp2,'(/)',ADVANCE='no')
  do iya=1,nya_act_WholeCore
    write(ctfinp2,'(i3)',ADVANCE='no') iya
    do ixa=1,nxa_act_WholeCore
      write (ctfinp2,'(i3)',ADVANCE='no') WholeCore_Conf(ixa,iya)
    enddo
    write(ctfinp2,*)
  enddo
  write(ctfinp2,'(i3)') IVARDX
  do i=1,IVARDX
    write(ctfinp2,'(f10.4,i5)') VARDX(i), JLEV(i)
  enddo
  write(ctfinp2,'(a12)') 'assembly.inp'

  !Assemlby card
  write(ctfinp3,'(i10)') Core%nAsyCell-Core%nAsyGT
  write(ctfinp3,'(i10)') nCellX0
  write(ctfinp3,'(i10)') Core%nAsyGT
  write(ctfinp3,'(f10.4)') hzact_tot*1000._8
  write(ctfinp3,'(A)') '{active region start}'
  write(ctfinp3,'(f10.4)') hzact_start*1000._8
  write(ctfinp3,'(f10.4)') ThVar%AsyPitch*1000._8
  write(ctfinp3,'(i10)') Casing
  write(ctfinp3,'(i10)')HTYPE
  write(ctfinp3,'(i10)')CMF
  write(ctfinp3,'(f10.4)') ThVar%rs*2._8*1000._8
  write(ctfinp3,'(i10)') ThVar%npr
  write(ctfinp3,'(f10.4)') ThVar%rgap*2._8*1000._8
  write(ctfinp3,'(f10.4)') ThVar%rw*2._8*1000._8
  write(ctfinp3,'(f10.4)') ThVar%ChannelPitch*1000._8
  write(ctfinp3,'(f10.4)') FTDENS
  write(ctfinp3,'(f10.4)') HGAP
  write(ctfinp3,'(a8)') 'Zircaloy'
  if (lgt) THEN
    write(ctfinp3,'(f10.4)') Thvar%rgti*2._8*1000._8
    write(ctfinp3,'(f10.4)') Thvar%rgto*2._8*1000._8
    write(ctfinp3,'(a8)') 'Zircaloy'
    do i=1,Core%nAsyGT
      write(ctfinp3,'(2i6)') gt_loc(i)%x, gt_loc(i)%y
    enddo
  endif
  if (flagspc) THEN        !!!
    write(ctfinp3,'(i3)') nspc
    do i=1,nspc
      write(ctfinp3,'(i5,3f10.4)') i, spcz(i)*1000, 0.0 , CDL(i)
    enddo
  ELSE
    write(ctfinp3,'(i3)') 0
  endif

  write(ctfinp4,'(f10.4)') nTracerCntl%PowerLevel*nTracerCntl%PowerFA*1.E-6_8*nasy_act_WholeCore
  write(ctfinp4,'(i5)')2
  write(ctfinp4,'(2f10.4)') 0.,  1.0
  write(ctfinp4,'(2f10.4)') hzact_tot*1000._8, 1.0
  do iya=1,nya_act_WholeCore
    do ixa=1,nxa_act_WholeCore
      if(WholeCore_Conf(ixa,iya) == 0) then
        write(ctfinp4,'(f10.4)',ADVANCE='no') 0.0 
      else
        write(ctfinp4,'(f10.4)',ADVANCE='no') 1.0 
      endif
    enddo
    write(ctfinp4,*)
  enddo
  close(ctfinp1)
  close(ctfinp2)
  close(ctfinp3)
  close(ctfinp4)
end subroutine


END MODULE














