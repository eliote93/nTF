module MATRATH_mod
  use param
  USE PARAM
  USE typedef, ONLY : CoreInfo_Type,       FmInfo_Type,        Pin_Type,        &
                      CmInfo_Type,         ThInfo_Type,        GroupInfo_Type,  &
                      Asy_Type,            AsyInfo_Type,       CoolantTH_Type,  &
                      PE_TYPE, FuelTh_Type, PinInfo_TYPE
  USE CNTL,    ONLY : nTRACERCntl_Type
  USE FILES,   ONLY : FILENAME, InputFileIdx
  USE TH_MOD,       ONLY : ThVar
  
  
  implicit none
  LOGICAL, SAVE :: lFirst = .TRUE.  
  logical, save :: lfinalmatra=.false.
  INTEGER, PRIVATE :: myzb, myze, nz, nz_act
  
!For channel geometry information  
  type cell_info
    integer :: globidx=0, radidx=0
    real*8 :: delx=0., dely=0., delz=0., area=0., vol=0.
  end type cell_info
  type, extends(cell_info) :: rod_type
    real*8 :: diam=0., radius, peri, relpow
    integer ::  asy_type
    integer :: nbch(6)=0

  end type rod_type
  
  type, extends(cell_info) :: vel_type
    real*8 :: l=0. ,slit=0., w(2)=0., ang=0.
    integer :: nbch(6)=0, nbvec(2)=0

  end type
  
  
  type nbrod_type
    class(rod_type),pointer :: rod
  end type nbrod_type
  
  type, extends(cell_info) :: ch_type
    real*8 :: wetperi=0., hyddiam=0., gap(4)
    integer:: nbch(6)=0, nbvec(4)=0, nnbch=0
    type(nbrod_type) :: nbrod(4)

  end type ch_type
  type(ch_type), allocatable :: ch(:)
  type(rod_type), allocatable, target :: rod(:)
  type(vel_type), allocatable :: vec(:)
  
  integer :: nchn, nrod, ngap, nchcol, nchrow
  real*8, allocatable :: delz(:)
  real*8 :: htot, areatot
  integer, allocatable :: chxy2rad(:,:)
  
  type zinfo_type
    real*8 :: delz=0.
    integer :: n=1
  end type
  type(zinfo_type), allocatable :: zinfo(:)
  
  
  type :: cell2Dinfo_type
    real*8 :: relpow=1.0
    real*8, allocatable :: pinrelpow(:)
    integer :: asy_type=0, globidx=0
    
  end type cell2Dinfo_type
  
  type bdcell_type
    integer :: istart=0, iend=0
  end type

  type(cell2Dinfo_type), allocatable :: core2Dinfo(:,:), pin2dinfo(:,:)
  type(bdcell_type), allocatable :: asyrow(:), asycol(:), pinrow(:), pincol(:), chrow(:), chcol(:)
  
  
  contains
  
  SUBROUTINE MATRA_TH(Core, FmInfo, CmInfo, ThInfo, GroupInfo, nTracerCntl,ToutAvg,tfmax, PE)
#ifdef MPI_ENV
  USE MPIComm_Mod, ONLY : BCAST, MPI_SYNC
#endif
  IMPLICIT NONE
#ifdef MPI_ENV
  include 'mpif.h'
#endif
  TYPE(CoreInfo_Type) :: Core
  TYPE(FmInfo_Type) :: FmInfo
  TYPE(CmInfo_Type) :: CmInfo
  TYPE(ThInfo_Type) :: ThInfo
  TYPE(GroupInfo_Type) :: GroupInfo
  TYPE(nTracerCntl_Type) :: nTracerCntl
  TYPE(PE_Type) :: PE
  REAL*8 :: ToutAvg, tfmax
  LOGICAL :: lidxmatch=.TRUE.
  integer :: iz
#ifdef HAVE_MATRA  
  if (lfirst) then
    nz_act =0
    DO iz = 1, Core%nz
      IF(.NOT. Core%lFuelPlane(iz)) CYCLE
      nz_act = nz_act + 1
    ENDDO
    
    IF(PE%CmfdMaster) THEN
      call MAKE_MATRA_INPUT(ThInfo, nTracerCntl)
      CALL MATRA_Initialize ()
    endif
    lfirst=.FALSE.
  ENDIF
#ifdef MPI_ENV
  CALL MPI_SYNC(PE%MPI_CMFD_COMM)	
#endif
  !CALL Set_MATRA_DATA()
  IF(PE%CmfdMaster) CALL MATRA_Standalone(lfinalmatra,lidxmatch)
  !CALL Get_MATRA_DATA()
#ifdef MPI_ENV
  CALL MPI_SYNC(PE%MPI_CMFD_COMM)	
#endif  
#endif  
  END SUBROUTINE
  

  !SUBROUTINE Get_MATRA_DATA()
  !INTEGER :: nrod, iz, i, izmatra
  !REAL*8 :: qflux
  !EXTERNAL :: MATRA_get_coolant_temp,MATRA_get_coolant_dens
  !nrod=nchan
  !
  !izmatra=1
  !DO iz=kfbega,kfenda
  !  
  !  DO i=1, nrod
  !    call MATRA_get_coolant_temp(i,izmatra,tcool(iz,i))
  !    call MATRA_get_coolant_dens(i,izmatra,dcool(iz,i))
  !  ENDDO
  !  izmatra=izmatra+1
  !ENDDO
  !
  !
  !
  !DO iz=kfenda+1,nzth+1
  !  DO i=1, nrod
  !    call MATRA_get_coolant_temp(i,izmatra,tcool(iz,i))
  !    call MATRA_get_coolant_dens(i,izmatra,dcool(iz,i))
  !  ENDDO
  !ENDDO
  !
  !
  !END SUBROUTINE
  !
  !SUBROUTINE Set_MATRA_DATA()
  !INTEGER :: nrod, iz, i, izmatra
  !REAL*8 :: qflux,roddiam
  !nrod=nchan
  !roddiam=2.*sqrt(npint*rw**2.+ngt*rgt**2.)
  !izmatra=0
  !DO iz=kfbega,kfenda
  !  izmatra=izmatra+1
  !  DO i=1, nrod
  !    qflux=plevel*powlin*relp(iz,i)/(2*rw*pi)*0.001
  !    !qflux=plevel*powlin*relp(iz,i)/(roddiam*pi)*0.001
  !    call MATRA_set_rodpower_kW_m2(i,izmatra,qflux)
  !  ENDDO
  !ENDDO
  !
  !
  !END SUBROUTINE
  
  
  SUBROUTINE MAKE_MATRA_INPUT(ThInfo, nTracerCntl)
  use files,    only: caseid
  use ioutil,   only: newunit
  INTEGER :: inpmatra
  INTEGER :: i, j, nline, m, nonunihz, k, ie, ise, isw,l, iz, nbgapidx, nbchidx
  INTEGER :: nCASE=1, J1=0
  INTEGER :: maxt=99999, iiunit=1, jiunit=1, jounit=1
  CHARACTER*10  IX1(1000)
  CHARACTER*10  X1(1000),X2(1000),X3(1000)
  CHARACTER*10 ndum, idum, kdum
  REAL*8 :: turfricfac1, turfricfac2, turfricfac3, lamfricfac1, lamfricfac2, lamfricfac3, roddiam
	REAL*8 :: pexit, ABETA, ttime,  flowarea, heatperi,wetperi, efflen, gapsize, BBETA
  real*8 :: pwfrac(4)=0.
  REAL*8, allocatable :: relhz(:), relzpower(:), unicellrelz(:), relcumhz(:), angle(:)
  real*8 :: eps(3), kij=0.5 , sl=0.5 , ftm=0., theta=0. , usdon=0. , dampng=0.8, accel(2)
  data eps / 0.1, 0.001, 0.001 /
  data accel / 1.6, 0.7 /
  integer, allocatable ::  unicelln(:)
  integer :: ndt, ntreis, itry, itrym
  INTEGER :: matracard1(8), matracard2(8), matracard3(2), matracard4(7), matracard7(7), matracard8(10), &
            matracard9(10), matracard10(7), matracard11(9),matracard12(9), matracard13(2)
  TYPE(ThInfo_Type) :: ThInfo
  TYPE(nTracerCntl_Type) :: nTracerCntl
  data matracard1 /1, 0, 0, 0, 1, 0, 0, 0/
  data matracard2 /2, 1, 3, 0, 0, 1, 0, 0/
  data matracard3 /3, 0/
  data matracard4 /4, 0, 0, 0, 0, 0, 0/
  data matracard7 /7, 1, 0, 0, 0, 0, 0/
  data matracard8 /8, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data matracard9 /9, 0, 0, 0, 0, 0, 0, 4, 0, 0/
  data matracard10 /10, 0, 0, 0, 0, 0, 0/
  data matracard11 /11, 1, 0, 0, 0, 0, 0, 0, 0/
  data matracard12 /12, 3, 0, 0, 0, 0, 0, 0, 0/
  data matracard13 /13, 1/


  turfricfac1=0.204
  turfricfac2=-0.20
  turfricfac3=0.
  
  lamfricfac1=64.
  lamfricfac2=-1.
  lamfricfac3=0.
  
  
  call init_SubChn_geo()
  
  matracard3(2)=nz_act+1
  
  
  allocate(relcumhz(0:nz_act),relzpower(0:nz_act))
  relcumhz=0.; relzpower=1.
  relcumhz(1)=delz(1)

  do iz=1, nz_act
    relcumhz(iz)=delz(iz)+relcumhz(iz-1)
  enddo

  
  
  relcumhz=relcumhz/THvar%Hact

  !do i=kfbega-1,kfenda
  !  !relzpower(i)=pi/2.*sin(pi*relcumhz(i))
  !  relzpower(i)=1.0
  !enddo
  
  !cal. or init. card no. 4 parameter
  matracard4(2)=nchn
  matracard4(3)=nchn
  matracard4(4)=0  

    
  !cal. or init. card no. 8 parameter
  matracard8(2)=nrod
  matracard8(3)=nrod
  
  
  nonunihz=0
  DO k=2,nz_act
    if(delz(k).EQ.delz(k-1)) cycle
    nonunihz=nonunihz+1
  ENDDO
  
  
  
  
  if(nonunihz.NE.0) THEN
    nonunihz=nonunihz+1
    allocate(unicelln(nonunihz))
    allocate(unicellrelz(nonunihz))
    k=1
    m=1
    DO i=2,nz_act
      if(delz(i).EQ.delz(i-1)) THEN
        m=m+1
      else
        unicelln(k)=m
        unicellrelz(k)=delz(i-1)/THvar%Hact
        m=1
        k=k+1
      ENDIF
    ENDDO
    unicelln(k)=m
    unicellrelz(k)=delz(i-1)/THvar%Hact
  ENDIF
  
  
  matracard9(7)=nonunihz
  
  ABETA=0.05
  BBETA=0.00

    
  ttime=0.

  ndt=0; ntreis=200; itry=150; itrym=0
  
  pexit = ThInfo%Pexit*1.E-6 !Mpa
  ! card1(6)   Heat transfer coeff. and friction factor
  !IFLUID=0 : Dittus-Boelter (for water)
  !                 IFLUID=1 : Schad-modified (CRBRP), (Nuclear System I pp. 451)
  !                 IFLUID=2 : Westinghouse (FFTF)
  !                 IFLUID=3 : Lyon-Martinelli (El-Wakil pp. 268)
      ! card1(7)   Flow split and de lp
  !                 NFSPLT=0 : Novendstern
  !                 NFSPLT=1 : Chiu-rohsenow-todreas
  !                 NFSPLT=2 : Cheng-todreas


   open(unit=newunit(inpmatra),file='matra.inp',status='unknown')
   write(inpmatra,'(4I10)') maxt,iiunit, jiunit, jounit
    write(inpmatra,'(2I10,4x,17A)') nCASE, J1, trim(caseid)
    write(inpmatra,'(8I10)') matracard1
    write(inpmatra,'(8I10)') matracard2
    write(inpmatra,'(3f10.3)') turfricfac1, turfricfac2, turfricfac3
    write(inpmatra,'(3f10.3)') lamfricfac1, lamfricfac2, lamfricfac3
  
    write(inpmatra,'(2I10)') matracard3
  
    nline=ceiling((nz_act+1)/5.)
    do j=1,nline-1
      write(inpmatra,'(10f10.5)') (relcumhz(i), relzpower(i), i=5*(j-1),5*j-1)
    enddo

    write(inpmatra,'(10f10.5)') (relcumhz(i), relzpower(i), i=5*(nline-1),nz_act)
    write(inpmatra,'(7I10)') matracard4
    do i=1,nchn
      flowarea=ch(i)%area*1.e+06
      wetperi=ch(i)%wetperi*1.e+03
      heatperi=wetperi
      write(inpmatra,"(I10,f10.4,2f10.4)",advance='no')  i, flowarea, wetperi, heatperi
      do k=1,2
        nbchidx=ch(i)%nbch(k)
        if (nbchidx/=0) then
          gapsize=ch(i)%gap(k)*1.e+03
          nbgapidx=ch(i)%nbvec(k)
          efflen=vec(nbgapidx)%l*1.e+03
          write (inpmatra,"(I10,f10.4,f10.4)",advance='no') nbchidx, gapsize, efflen
        endif
      enddo
      write (inpmatra,*)

    enddo

    write(inpmatra,'(10I10)') matracard8
    DO I=1,nrod
      do k=1,4
      if (rod(I)%nbch(k)==0) then
        pwfrac(k)=0.
      else
        pwfrac(k)=0.25
      endif
    enddo
      write(inpmatra,2280)  I, rod(I)%diam*1.E+03, 1.,(rod(I)%nbch(k), pwfrac(k), k=1,4)
    ENDDO

    write(inpmatra,'(10I10)') matracard9
    write(inpmatra,'(2f10.3,3es10.1e2,8f10.3)')THvar%Hact*1000., ttime, eps, kij, sl, ftm, theta, usdon, dampng, accel
    write(inpmatra,2340) nz_act, ndt, ntreis, itry, itrym
    if (matracard9(7).GT.1) THEN
      WRITE (inpmatra,2370) (unicelln(L),unicellrelz(L),L=1,matracard9(7))
      deallocate(unicelln,unicellrelz)
    ENDIF
    
    write(inpmatra,'(7I10)') matracard10
    write(inpmatra,'(2f10.3)') ABETA, BBETA

    
    write(inpmatra,'(9I10)') matracard11
    write(inpmatra,'(f10.3,3f10.3)') pexit, ThInfo%tin, nTracerCntl%fMdotFA/(ThVar%acf*ThVar%nAsyCh),  &
                                      ThInfo%PowLv*ThInfo%powlin/(Thvar%zeta)*0.001
    write(inpmatra,'(9I10)') matracard12
    write(inpmatra,'(2I10)') matracard13
    write(inpmatra,'(/)')
    close(inpmatra)

    deallocate(relzpower,relcumhz)
    
    

2140 FORMAT (I10,f10.2, 2f10.3)
     

2250 FORMAT (10(I10))
2280 FORMAT (I10,f10.4,f10.4,4(I10,f10.4))
2330 FORMAT (f10.3,f10.5,10(10X),f10.5)
2340 FORMAT (5I10)
2370 FORMAT (10(I10,f10.7))
  END SUBROUTINE

  subroutine init_SubChn_geo
  implicit none
  
    call gen_scalarcell
    call gen_vectorcell
    !call gen_powinfo
  end subroutine




subroutine gen_scalarcell

  use geom,         only: nCellX0, core, asypitch, cellpitch
  use TH_MOD,       only: ThVar
  use ioutil,       only: newunit
  implicit none
  integer :: i,j,m,n, output, iyc, ixc, ixr, iyr, ixa,iya,output2, k, nass=0, flag
  integer :: iycs, iyce, ixcs, ixce, ixrs, ixre, nnbrod

  real*8 :: rodarea, asygap, ppitch, drod, asypitch_m
  integer ::  nxa, nya, npin, nrodrow, nrodcol, npair, ixya, iz
  logical :: lasygap, lwall=.false.

  npin=nCellX0
  nxa=core%nxa
  nya=core%nya
  drod=2.*ThVar%rw
  lasygap=core%lgap
  ppitch=cellpitch*1.E-2
  asypitch_m=asypitch*1.E-2
  nrodrow=npin*nxa
  nrodcol=npin*nya
  
  nchrow=npin*nxa+1
  nchcol=npin*nya+1

  allocate(delz(0:nz_act+1))
  delz=0.
  k=0
  do iz=1, Core%nz
    IF(.NOT. Core%lFuelPlane(iz)) CYCLE
    k=k+1
    delz(k)=Core%hz(iz)*0.01_8
  enddo
  
  allocate (core2Dinfo(0:nxa+1, 0:nya+1))


  m=0
  do j=1,nya
    do i=1,nxa
      ixya=core%coreIdx(nxa, nya)
      if (ixya == 0) cycle
      core2Dinfo(i,j)%asy_type=core%coremap(ixya)
!      if(core2Dinfo(i,j)%asy_type == 0) cycle
      m=m+1
      core2Dinfo(i,j)%globidx=m
      nass=max(nass,core2Dinfo(i,j)%globidx)
    enddo
  enddo
  

  allocate(asyrow(nya),asycol(nxa))
  
  flag=0
  do j=1,nya
    i=1
    do while (flag==0)
      
      if(core2Dinfo(i,j)%asy_type /= 0) then
        asyrow(j)%istart=i
        flag=flag+1
      endif
      i=i+1
    enddo
    flag=0
  enddo
  
  
  do j=1,nya
    i=nxa
    do while (flag==0)
      
      if(core2Dinfo(i,j)%asy_type /= 0) then
        asyrow(j)%iend=i
        flag=flag+1
      endif
      i=i-1
    enddo
    flag=0
  enddo
  
  do i=1,nxa
    j=1
    do while (flag==0)
      
      if(core2Dinfo(i,j)%asy_type /= 0) then
        asycol(i)%istart=j
        flag=flag+1
      endif
      j=j+1
    enddo
    flag=0
  enddo
  
  do i=1,nxa
    j=nya
    do while (flag==0)
      
      if(core2Dinfo(i,j)%asy_type /= 0) then
        asycol(i)%iend=j
        flag=flag+1
      endif
      j=j-1
    enddo
    flag=0
  enddo
  
  
  allocate(pinrow(nrodcol),pincol(nrodrow))
  nrod=0
  do iyr=1,nrodcol
    iya=(iyr-1)/npin+1
    pinrow(iyr)%istart=(asyrow(iya)%istart-1)*npin+1
    pinrow(iyr)%iend=(asyrow(iya)%iend)*npin
    nrod=nrod+pinrow(iyr)%iend-pinrow(iyr)%istart+1
  enddo
  
  do ixr=1,nrodrow
    ixa=(ixr-1)/npin+1
    pincol(ixr)%istart=(asycol(ixa)%istart-1)*npin+1
    pincol(ixr)%iend=(asycol(ixa)%iend)*npin
  enddo
  
  
  
  allocate(pin2Dinfo(nrodrow, nrodcol))
  do iyr=1,nrodcol
    do ixr=1,nrodrow
      i=(ixr-1)/npin+1
      j=(iyr-1)/npin+1
      pin2Dinfo(ixr,iyr)%asy_type=core2Dinfo(i,j)%asy_type
    enddo
  enddo
  
  
  
  
  
  allocate(chrow(nchcol),chcol(nchrow))
  
  
  
  
    
  iyce=0
  do iya=1,nya
    iycs=iyce+1
    if (core2Dinfo(1,iya)%asy_type/=0 .and. core2Dinfo(1,iya+1)%asy_type==0) then
      iyce=iyce+npin+1
    else
      iyce=iyce+npin
    endif
    chrow(iycs:iyce)%istart=(asyrow(iya)%istart-1)*npin+1
    chrow(iycs:iyce)%iend=(asyrow(iya)%iend)*npin+1
  enddo
  
  
  ixce=0
  do ixa=1,nxa
    ixcs=ixce+1
    if (core2Dinfo(ixa,1)%asy_type/=0 .and. core2Dinfo(ixa+1,1)%asy_type==0) then
      ixce=ixce+npin+1
    else
      ixce=ixce+npin
    endif
    chcol(ixcs:ixce)%istart=(asycol(ixa)%istart-1)*npin+1
    chcol(ixcs:ixce)%iend=(asycol(ixa)%iend)*npin+1
  enddo
  
  
  
  nchn=0
  do iyc=1,nchcol
    nchn=nchn+chrow(iyc)%iend-chrow(iyc)%istart+1
  enddo
  
  
  allocate (ch(nchn), chxy2rad(0:nchrow+1, 0:nchcol+1))
  allocate(rod(0:nrod))
  
  chxy2rad=0
  
  
  if (lasygap) then
    asygap=asypitch_m-ppitch*npin
  else
    asygap=0.
  endif
  
  I=0
  do iyc=1,nchcol
    ixcs=chrow(iyc)%istart
    ixce=chrow(iyc)%iend
    do ixc=ixcs,ixce
      I=I+1
      iycs=chcol(ixc)%istart
      iyce=chcol(ixc)%iend
      ch(I)%globidx=I
      chxy2rad(ixc,iyc)=I
      if (ixc==ixcs .or. ixc==ixce) then
        if (iyc>=iycs .and. iyc<=iyce) then
          ch(I)%delx=(ppitch+asygap)/2.
        else
          ch(I)%delx=ppitch+asygap          
        endif
    
      elseif  (mod(ixc,npin)==1) then
         ch(I)%delx=ppitch+asygap
      
      else
        ch(I)%delx=ppitch
      endif
  
    enddo
  enddo
  
  
  
  do ixc=1,nchrow
    iycs=chcol(ixc)%istart
    iyce=chcol(ixc)%iend
    do iyc=iycs,iyce
      I=chxy2rad(ixc,iyc)
      ixcs=chrow(iyc)%istart
      ixce=chrow(iyc)%iend
      if (iyc==iycs .or. iyc==iyce) then
        if (ixc>=ixcs .and. ixc<=ixce) then
          ch(I)%dely=(ppitch+asygap)/2.
        else
          ch(I)%dely=ppitch+asygap
        endif
    
      elseif  (mod(iyc,npin)==1) then
         ch(I)%dely=ppitch+asygap
      
      else
        ch(I)%dely=ppitch
      endif
      ch(I)%nbrod(1)%rod=>rod(0)
      ch(I)%nbrod(2)%rod=>rod(0)
      ch(I)%nbrod(3)%rod=>rod(0)
      ch(I)%nbrod(4)%rod=>rod(0)
    enddo
  enddo
  
  
  
  I=0
  
  do iyc=1,nchcol
    ixcs=chrow(iyc)%istart
    ixce=chrow(iyc)%iend
    do ixc=ixcs,ixce
      I=I+1
      ch(I)%nbch(1)=chxy2rad(ixc+1,iyc)
      ch(I)%nbch(2)=chxy2rad(ixc,iyc+1)
      ch(I)%nbch(3)=chxy2rad(ixc-1,iyc)
      ch(I)%nbch(4)=chxy2rad(ixc,iyc-1)
      n=0
      do k=1,4
        if (ch(I)%nbch(k)==0) cycle
        n=n+1
      enddo
      ch(I)%nnbch=n
    enddo
  enddo
  


  
  I=0
  do iyr=1,nrodcol
    ixrs=pinrow(iyr)%istart
    ixre=pinrow(iyr)%iend
    do ixr=ixrs,ixre
      I=I+1
      rod(I)%globidx=I
      rod(I)%diam=drod
      rod(I)%radius=drod/2.
      rod(I)%area=rod(I)%radius**2.*pi
      rod(I)%peri=rod(I)%diam*pi
      rod(I)%asy_type=pin2Dinfo(ixr,iyr)%asy_type
      rod(I)%nbch(1)=chxy2rad(ixr,iyr)
      rod(I)%nbch(2)=chxy2rad(ixr+1,iyr)
      rod(I)%nbch(3)=chxy2rad(ixr,iyr+1)
      rod(I)%nbch(4)=chxy2rad(ixr+1,iyr+1)
      ch(chxy2rad(ixr,iyr))%nbrod(2)%rod=>rod(I)
      ch(chxy2rad(ixr+1,iyr))%nbrod(3)%rod=>rod(I)
      ch(chxy2rad(ixr,iyr+1))%nbrod(1)%rod=>rod(I)
      ch(chxy2rad(ixr+1,iyr+1))%nbrod(4)%rod=>rod(I)
    enddo
  enddo

  
 
  do I=1,nchn
    nnbrod=0
    rodarea=0.
    ch(I)%wetperi=0.
    do k=1,4
      rodarea=rodarea+0.25*ch(I)%nbrod(k)%rod%area
      ch(I)%wetperi=ch(I)%wetperi+0.25*ch(I)%nbrod(k)%rod%peri
      if (ch(I)%nbrod(k)%rod%globidx/=0) nnbrod=nnbrod+1
    enddo


    ch(I)%gap(1)=ch(I)%dely-ch(I)%nbrod(1)%rod%radius-ch(I)%nbrod(2)%rod%radius !east
    ch(I)%gap(2)=ch(I)%delx-ch(I)%nbrod(2)%rod%radius-ch(I)%nbrod(3)%rod%radius !south
    ch(I)%gap(3)=ch(I)%dely-ch(I)%nbrod(3)%rod%radius-ch(I)%nbrod(4)%rod%radius !west
    ch(I)%gap(4)=ch(I)%delx-ch(I)%nbrod(4)%rod%radius-ch(I)%nbrod(1)%rod%radius !north
    ch(I)%area=ch(I)%delx*ch(I)%dely-rodarea
    if (nnbrod==3) then
      if (ch(I)%nbrod(1)%rod%globidx==0) then
        ch(I)%gap(1)=0.5*ch(I)%dely-ch(I)%nbrod(2)%rod%radius !east
        ch(I)%gap(4)=0.5*ch(I)%delx-ch(I)%nbrod(4)%rod%radius !north
      elseif (ch(I)%nbrod(2)%rod%globidx==0) then
        ch(I)%gap(1)=0.5*ch(I)%dely-ch(I)%nbrod(1)%rod%radius !east
        ch(I)%gap(2)=0.5*ch(I)%delx-ch(I)%nbrod(3)%rod%radius !south
      elseif (ch(I)%nbrod(3)%rod%globidx==0) then
        ch(I)%gap(2)=0.5*ch(I)%delx-ch(I)%nbrod(2)%rod%radius !south
        ch(I)%gap(3)=0.5*ch(I)%dely-ch(I)%nbrod(4)%rod%radius !west      
      else
        ch(I)%gap(3)=0.5*ch(I)%dely-ch(I)%nbrod(3)%rod%radius !west
        ch(I)%gap(4)=0.5*ch(I)%delx-ch(I)%nbrod(1)%rod%radius !north
      endif
      
      ch(I)%area=0.75*ch(I)%delx*ch(I)%dely-rodarea
 
     endif
    if (lwall) then
      select case (nnbrod)
        case (1)
          ch(I)%wetperi=ch(I)%wetperi+ch(I)%delx+ch(I)%dely
        case (2)
          ch(I)%wetperi=ch(I)%wetperi+max(ch(I)%delx,ch(I)%dely)
        case (3)
          ch(I)%wetperi=ch(I)%wetperi+(ch(I)%delx+ch(I)%dely)*0.5
      
      end select
    
    endif
    
      
  enddo

  
  open(unit=newunit(output),file='scalarcell',status='unknown')
  write(output,'(A)') 'CH no, x id. y id.       delx,        dely        area         wetper       '
  do iyc=1,nchcol
    ixcs=chrow(iyc)%istart
    ixce=chrow(iyc)%iend
    do ixc=ixcs,ixce
      I=chxy2rad(ixc,iyc)
      write(output,'(3I5,4ES21.10)') I, ixc, iyc, ch(I)%delx, ch(I)%dely, ch(I)%area, ch(I)%wetperi
    enddo
  enddo
  do J=1,nz_act
    write(output, '(I5,E17.9)') J, delz(J)
  enddo
  
  close(output)
  
  npair=1
  do J=2,nz_act
    if(delz(J-1)==delz(J)) cycle
    npair=npair+1
  
  enddo
  
  allocate(zinfo(npair))
  zinfo(1)%delz=delz(1)
  npair=1
  do J=2,nz_act
    if(delz(J-1)==delz(J)) then
      zinfo(npair)%delz=delz(J)
      zinfo(npair)%n=zinfo(npair)%n+1
      cycle
    endif
    npair=npair+1
  enddo
  
  
  htot=0.
  do J=1,nz_act
    htot=htot+delz(J)
  enddo
  
  areatot=0.
  do I=1,nchn
    areatot=areatot+ch(I)%area
  enddo  
  
  
  


end subroutine gen_scalarcell


subroutine gen_vectorcell
  
  implicit none
  integer :: iyc, ixc, ixce, ixcs, icell, icellp1, iyce, iycs, igap
!  
  ngap=0
  do iyc=1,nchcol
    ngap=ngap+chrow(iyc)%iend-chrow(iyc)%istart
  enddo
  
  do ixc=1,nchrow
    ngap=ngap+chcol(ixc)%iend-chcol(ixc)%istart
  enddo
  
  allocate(vec(ngap))
  igap=0
  do iyc=1,nchcol
    ixcs=chrow(iyc)%istart
    ixce=chrow(iyc)%iend
    do ixc=ixcs,ixce-1
      igap=igap+1
      vec(igap)%globidx=igap
      icell=chxy2rad(ixc,iyc)
      icellp1=chxy2rad(ixc+1,iyc)
      vec(igap)%l=0.5*(ch(icell)%delx+ch(icellp1)%delx)
      vec(igap)%w(1)=0.5*ch(icell)%delx
      vec(igap)%w(2)=0.5*ch(icellp1)%delx
      vec(igap)%slit=ch(icell)%gap(1)
      vec(igap)%nbch(1)=icell
      vec(igap)%nbch(2)=icellp1
      vec(igap)%ang=0.
      ch(icell)%nbvec(1)=igap
      ch(icellp1)%nbvec(3)=igap
    enddo
  enddo
    
  
  
  do ixc=1,nchrow
    iycs=chcol(ixc)%istart
    iyce=chcol(ixc)%iend
    do iyc=iycs,iyce-1
      igap=igap+1
      vec(igap)%globidx=igap
      icell=chxy2rad(ixc,iyc)
      icellp1=chxy2rad(ixc,iyc+1)
        
      vec(igap)%l=0.5*(ch(icell)%dely+ch(icellp1)%dely)
      vec(igap)%w(1)=0.5*ch(icell)%dely
      vec(igap)%w(2)=0.5*ch(icellp1)%dely
      vec(igap)%slit=ch(icell)%gap(2)
      vec(igap)%nbch(1)=icell
      vec(igap)%nbch(2)=icellp1
      vec(igap)%ang=270.
      ch(icell)%nbvec(2)=igap
      ch(icellp1)%nbvec(4)=igap
    enddo
  enddo


  

end subroutine gen_vectorcell
  
  
  
  
  
  
  
END MODULE