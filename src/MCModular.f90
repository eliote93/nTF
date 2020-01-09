#include <defines.h>
module MCModular
    use MCDefine
    USE TYPEDEF, ONLY : FxrInfo_Type
    implicit none
    !private
    !public :: InitMCModular, MCSimulation, SetMCEnv, SetMCCMFDEnv
    integer :: mcng
    logical :: lic, lgeninitsrc, lomp, lonfdb, ladjusted
    real(8) :: trk_gk, wtnorm
    type(MCCTRL) :: ctrlMC
    integer(8) :: gseed
    TYPE(FxrInfo_Type), POINTER :: FXR(:,:)
    INTEGER :: nfxr, ncorez

    !logical :: lbenchXS
contains

subroutine SetMCEnv(nht, ninact, nact)
    implicit none
    integer, intent(in) :: ninact, nact, nht
    ctrlMC%ninact=ninact
    ctrlMC%nact=nact
    ctrlMC%nht=nht
end subroutine

subroutine SetMCCMFDEnv(lcmfd, lfdb)
    implicit none
    logical, intent(in) :: lcmfd, lfdb(2)
    ctrlMC%lcmfd=lcmfd
    ctrlMC%lfdb=lfdb
end subroutine

!subroutine ReadMC_CARD(io)
!    implicit none
!    integer, intent(in) :: io
!    character*512 oneline
!    character(10) :: astring
!    read(indev, '(a256)') oneline
!    read(oneline), astring
!    select case(astring)
!        case("MCENV")
!end subroutine
!
subroutine InitMCModular()
!$  use OMP_LIB
    use typedef
    use rng
    use MCDefine
    USE Core_mod,   ONLY: FMInfo, THInfo, GroupInfo
    use Geom,       only: Core, AsyInfo, CellInfo, PinInfo, asypitch, albedo, ledge, ng, hz
    use BenchXs,    only: ngben
    use MCBasics,   only: InitMCBasics
    use MCTally,    only: InitMCTally
    use CMFDTally,  only: InitCMFDTally
    use MCMap,      only: InitMCMap, GetTotNMod
    use Xsec4MC,    only: InitXsec4MC
    use MCFuelRect, only: InitMCFuelRect
    use MCGapRect,  only: InitMCGapRect
    use MCCMFD,     only: InitMCCMFD, SetCMFDSet
    use MCIterLog,  only: InitMCIterLog
    USE CNTL,             ONLY : nTracerCntl
    USE PE_MOD,           ONLY : PE
    USE RAYS,             ONLY : RayInfo
    USE RNG,            ONLY : InitRN
    USE MCXsLib
    USE InitMC_MOD
    USE SubGrp_mod,      ONLY : SubGrpFsp, SubGrpEffXsGen
    USE Files,     ONLY : io8
    USE IOUTIL,    ONLY : message
    implicit none
    !type(coreinfo_type), intent(in) :: core
    !type(asyinfo_type), pointer, intent(in), dimension(:) :: asyinfo
    !type(pininfo_type), pointer, intent(in), dimension(:) :: pininfo
    !type(cell_type), pointer, intent(in), dimension(:) :: cellinfo
    !type(FmInfo_type) :: Fminfo
    LOGICAL :: lbenchXS, lmicxs

  IF(PE%Master) THEN
    mesg = 'Initialize the MC module'
    CALL message(io8,TRUE, TRUE, mesg)
  ENDIF
    ctrlMC%nthrd=PE%nthread
    IF(ctrlMC%nthrd .EQ. 1)THEN
        lomp=.FALSE.
    ELSE
        lomp=.TRUE.
    ENDIF
    ctrlMC%lomp=lomp
    ctrlMC%lmicXS=.FALSE.
    ctrlMC%lmicXS=.TRUE.


!$  call omp_set_dynamic(.FALSE.)
!$  call omp_set_num_threads(ctrlMC%nthrd)
    lbenchXs=nTracerCntl%lbenchXs
    lmicxs=ctrlMC%lmicXS

    CALL InitRN()
      CALL SetMCEnv(nht, ninact, nact)
      IF( lMCCmfd )THEN
        CALL SetMCCMFDEnv(lmccmfd, lfdb)
        IF( lMCcmfdSet )THEN
            CALL SetCMFDSet(idx, period, skip, flag)
        ENDIF
      ENDIF

    gseed=getseed()

    ctrlMC%lXslib=.not.lBenchXs
    IF(lBenchXs)THEN
        ctrlMC%ngMC=ngben
    ELSE
        ctrlMC%ngMC=ng
    ENDIF
    ctrlMC%ngCMFD=2
    if (ctrlMC%ngMC .eq. 7) then
        ctrlMC%gThrm=4
    elseif (ctrlMC%ngMC .eq. 47) then
        ctrlMC%gThrm=35
    elseif (ctrlMC%ngMC .eq. 281)THEN
        ctrlMC%gthrm=142
    else
        ctrlMC%gThrm=ctrlMC%ngMC
    endif
    FXR=>FmInfo%FXR
    mcng=ctrlMC%ngMC
    ctrlMC%lsgfsp=nTracerCntl%lrestrmt
    IF(ctrlMC%lXslib)THEN
      IF(PE%Master) THEN
        mesg = 'Preparing XS for MC using nTRACER library'
        CALL message(io8,TRUE, TRUE, mesg)
      ENDIF
        IF(nTracerCntl%lrestrmt) THEN
            write(*,*) 'Entering SubgrpFSP'
            CALL SubGrpFsp(Core, FmInfo%Fxr, THInfo, RayInfo, GroupInfo, nTracerCntl, PE)
            write(*,*) 'Entering SubgrpEffXsGen'
            CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, 1._8, GroupInfo, nTracerCntl, PE)
        ENDIF
      IF(PE%Master) THEN
        mesg = 'Setting Macroscopic XS for MC'
        CALL message(io8,TRUE, TRUE, mesg)
      ENDIF
        IF( ctrlMC%lmicXS )THEN
            CALL SetMcXslib(Core, FmInfo%Fxr, PE, ctrlMc)
        ELSE
            CALL SetMCXsLib_old(Core, FmInfo%Fxr, PE, ctrlMc)
        ENDIF

    ELSE
      IF(PE%Master) THEN
        mesg = 'Setting Macroscopic XS for MC'
        CALL message(io8,TRUE, TRUE, mesg)
      ENDIF
        call InitXsec4MC(ctrlMC%ngMC, ctrlMC%gThrm)
    ENDIF
      IF(PE%Master) THEN
        mesg = 'Setting MC variables'
        CALL message(io8,TRUE, TRUE, mesg)
      ENDIF
    call InitMCMap(Core, Asyinfo, PinInfo, Cellinfo, FmInfo, ctrlMC, asypitch, albedo, ledge, lbenchXS, lmicxs, ng)
    call InitMCBasics(ctrlMC%nht, ctrlMC%nxyz, albedo, lomp, ctrlMC%nthrd)
    call InitMCTally(ctrlMC%nxyz,nfxr,ng,ctrlMC%nthrd, ctrlMC%lXsLib, ctrlMC%lmicxs)
    IF (ctrlMC%lCMFD) THEN
      call InitCMFDTally(ctrlMC%nxyzCMFD, ctrlMC%ngCMFD, ctrlMC%nthrd, ctrlMC%ngMC, asypitch, hz, ctrlMC%CoreMap, ctrlMC%cutxy)
      call InitMCCMFD(ctrlMC%nxyzCMFD, albedo, ctrlMC%lfdb)
    ENDIF
    call InitMCGapRect(ctrlMC%lcmfd)
    call InitMCFuelRect(ctrlMC%lcmfd, ctrlMC%lXslib, ctrlMC%lMicXs)
    call InitMCIterLog(ctrlMC%lcmfd)

      IF(PE%Master) THEN
        mesg = 'Setting MC END'
        CALL message(io8,TRUE, TRUE, mesg)
      ENDIF
    lonfdb=.false.


end subroutine

subroutine MCSimulation
    use MCDefine
    use MCMap
    use MCBasics
    use MCTally
    use CMFDTally
    use MCIterLog
    type(STR_RectFuel), pointer :: mpFR
    type(nfs) :: nf
    integer :: icyc, nsrc
    real(8) :: ts, te
    real(8), pointer, dimension(:,:,:) :: pPPow, pFSD

    call AdrPPow(pPPow)
    !if (.not. lgeninitsrc) then
    !    call GenInitialSource()
    !    lgeninitsrc=.false.
    !endif

    call TallyTurnOff()
    ctrlMC%lFirstCycle=.TRUE.
    do icyc=1, ctrlMC%ninact
        ts=timecheck()
        call SimulateSingleCycle()
        te=timecheck()
        call AssignLogMC()
        call PrepareNextCycle(.false., icyc)
        call WriteLogMC()
    enddo

    call TallyTurnOn()
    IF (ctrlMC%lCMFD) call ResetCMFDTally()
    call ResetLogK()
    print *, "Active Cycle"
    do icyc=1, ctrlMC%nact
        ts=timecheck()
        call SimulateSingleCycle()
        te=timecheck()
        call AssignLogMC()
        call AcmPinPow(pPPow)
        call PrepareNextCycle(.true., icyc)
        call WriteLogMC()
    enddo
    call WritePinPow(ctrlMC)
contains

subroutine AssignLogMC()
    call SetLogCycle(icyc)
    call SetLogNSSRC(getNQueueCur())
    call SetLogTCYC(te-ts)
    call SetLogKeff(GetMCK())
!$  call CollectTallies()
    call SetLogShn(ShannonEntropy())
end subroutine
end subroutine


subroutine PrepareNextCycle(lact, icyc)
    use rng
    use MCBasics
    use CMFDTally
    use MCCMFD
    use MCIterLog
    USE MCTally, ONLY : CollectSpecTally
    USE MCXsLib, ONLY : UpdtEffchi
    implicit none
    logical, intent(in) :: lact
    integer, intent(in) :: icyc
    integer :: itercmfd, icmfd
    real(8) :: ts, te, errcmfd, kcmfd, shncmfd
    logical :: lreset

    ctrlMC%lFirstCycle=.FALSE.
    if (ctrlMC%lcmfd) then
        ladjusted=lonfdb
        ts=timecheck()
        call SetKeff4NBal(GetMCK())
        call CollectTally()
        call GenHomXsec()
        call CMFD4MC(kcmfd, itercmfd, errcmfd, shncmfd)
        call CMFDCTRL(lact, icyc, icmfd, lonfdb, lreset) ! lonfdb is module variable
        if (lonfdb) call CalWtAdj(ctrlMC%nht)
        if (lreset) call ResetCMFDTally()
        call SetZeroCMFDTally()
        te=timecheck()
        call SetLogCMFD(icmfd, lonfdb, te-ts, kcmfd, shncmfd)
    endif
    if( ctrlMC%lxslib .AND. .NOT. ctrlMC%lmicXS )THEN
        CALL CollectSpecTally()
        CALL UpdtEffChi()
    ENDIF

    gseed=strideN(gseed,GetNqueueCur())
    call ResetKeffTally()
    call ResetFSD()
    call ResetQueue()
!    call NormalizeQueue()
end subroutine

subroutine SimulateSingleCycle()
    use rng
    use MCDefine
    use MCBasics
    USE Cntl, ONLY : nTracerCntl
!$  use omp_lib
    implicit none
    integer :: i,j
    type(nfs) :: nf
    type(nstat) :: nst
    integer(8) :: lseed
    INTEGER :: ifxr, iz
    TYPE(FxrInfo_Type), POINTER :: myFXR

    lseed=gseed
    ! Fission neutron weight re-normalization
    wtnorm=1.*ctrlMC%nht/getNQueueCur()

!$  if (lomp) then
!$omp parallel default(shared) &
!$omp private(i, nf, nst)
!$omp do schedule (static)
!$  do i=1, getNQueueCur()
!$      nf=GetFSNXT(i)
!$      nst%tid= omp_get_thread_num()+1
!$      nst%seed=strideN(lseed, i-1)
!$      nst%lmn(:)=nf%lmn(:)
!$      nst%surf=0
!$      nst%axis=0
!$      nst%pr=nf%pr
!$      nst%ir=nf%ir
!$      nst%wt=nf%wt
!$      nst%iso=nf%iso
!$      call SimulateSingleNeutron(nst)
!$  enddo
!$omp end do
!$omp end parallel
!$ else
    do i=1, getNQueueCur()
       !IF( i.EQ.34 )THEN
       !    j=i
       !ENDIf

        nf=GetFS(i)
        nst%tid=1
        nst%seed=strideN(lseed, i-1)
        nst%lmn(:)=nf%lmn(:)
        nst%surf=0
        nst%axis=0
        nst%pr=nf%pr
        nst%ir=nf%ir
        nst%wt=nf%wt
        nst%iso=nf%iso
!        IF (i .eq. 1) THEN
!          nTracerCntl.lMCcycle = .TRUE.
!        ELSE
!!          nTracerCntl.lMCcycle = .FALSE.
!        END IF
       !write(*,*) i, nf%lmn(:)
        call SimulateSingleNeutron(nst)
    enddo
!$ endif
end subroutine

subroutine SimulateSingleNeutron(nst)
    use MCDefine
    use MCMap
    use MCBasics
    use MCFuelRect
    use MCGapRect
    use MCVacRect
    use CMFDTally
    use MCCMFD
    USE HighOrderSC
    USE MCXsLib, ONLY : fxr2xsmap, xs2fxrmap
    implicit none
    type(nstat) :: nst
    type(STR_RectFuel), pointer :: mpFR
    type(STR_RectGap), pointer :: mpGR
    type(STR_RectVac), pointer :: mpVR
    type(MCtallyset), pointer :: tallyset
    type(XsecSet), pointer :: xs
    integer :: iret, dcyl, dirneigh, axisneigh, tneigh
    integer :: surfrpt
    integer :: lmnorg(3)
    logical :: term, lcmfd
    real(8) :: dts
    TYPE(FxrInfo_Type), POINTER :: myFXR
    INTEGER :: ifxr, iz, ct

    term=.false.
    lcmfd=ctrlMC%lcmfd

    if(TypePin(nst%lmn) .eq. FRECT) then

        call AdrRectFuel(nst%lmn, mpFR)
        ! Cross sections are required for sampling neutron energy group and direction

        xs=>mpFR%xs(nst%ir)
        nst%pa=mpFR%cent+nst%pr

        ! fission neutron's weight normailzation
        if (lonfdb) then
            nst%wt=nst%wt*GetAdjFactor(mpFR%infoCMFD%lmn)
        else
            nst%wt=nst%wt*wtnorm
        endif
        ! Adjust FS's position
        if (mpFR%cut) then
            if (mpFR%cutxy(1)) nst%pr(1)=abs(nst%pr(1))  ! cutline = y-axis
            if (mpFR%cutxy(2)) nst%pr(2)=-abs(nst%pr(2)) ! cutline = x-axis
            nst%pa=mpFR%cent+nst%pr
        endif
    endif

    IF( ctrlMC%lXslib .AND. .NOT. ctrlMC%lFirstCycle .AND. ctrlMC%lmicXs )THEN  !microscopic
        nst%g=SmpGroupIso(nst, xs)
        !nst%g=SmpGroup(nst, xs.chi)
    ELSE  !macroscopic chi
        nst%g=SmpGroup(nst, xs%chi)
    ENDIF

    !nst%g=1
    call SmpDirIso(nst, nst%dir)

    !CALL SmpDirAnIsoGaussian(nst, nst%dir, xs)
    nst%rd=OUTDIR
    ct=0
    do while(.true.)
        ! MC simulation in modular cells
        if (TypePin(nst%lmn) .eq. FRECT) then ! fuel
            if (lcmfd) call CurrentTallyIn(nst, mpFR%infoCMFD) ! CMFD incoming surface current Tally
            term=MCRectFuel(nst, mpFR)
            if (lcmfd .and. .not. term) call CurrentTallyOut(nst, mpFR%infoCMFD) ! CMFD outgoing surface current Tally
        elseif (TypePin(nst%lmn) .eq. GRECT) then ! gap
            if (lcmfd) call CurrentTallyIn(nst, mpGR%infoCMFD)
            term=MCRectGap(nst, mpGR)
            if (lcmfd .and. .not. term) call CurrentTallyOut(nst, mpGR%infoCMFD)
        elseif (TypePin(nst%lmn) .eq. VRECT) then ! Vacuum
            call MCRectVac(nst, mpVR)
            term=.false.
        else

        endif
        if (term) then
            exit
        else
            ! find new modular cell
            lmnorg=nst%lmn
            tneigh=MCNeighbor(nst)
            if (tneigh > 0) then ! move to next
                if(TypePin(nst%lmn) .eq. FRECT) then
                    call AdrRectFuel(nst%lmn, mpFR)
                elseif(TypePin(nst%lmn) .eq. GRECT) then
                    call AdrRectGap(nst%lmn, mpGR)
                elseif(TypePin(nst%lmn) .eq. VRECT) then
                    call AdrRectVac(nst%lmn, mpVR)
                else
                    print *, "[ERR] Next rectangular fuel cannot be found!"
                endif
            ! treatment for boundary condition
            elseif (tneigh .eq. REFL) then ! reflective
                nst%lmn=lmnorg
                call SetDirRefl(nst)
                nst%surf=SurfOp(nst%surf, nst%axis)
            elseif(tneigh .eq. RPT) then ! repeated geometry for checkerboard
                nst%lmn=lmnorg
                tneigh=MCNeighborRepeated(nst%lmn, nst%surf, nst%axis)
                if (tneigh > 0) then ! move to next
                    if(TypePin(nst%lmn) .eq. FRECT) then
                        call AdrRectFuel(nst%lmn, mpFR)
                        surfrpt=GetSurfOp(nst%surf, nst%axis)
                        call Move2SurfRect(nst%pa, mpFR%cent, mpFR%slen, surfrpt, nst%axis)
                    elseif(TypePin(nst%lmn) .eq. GRECT) then
                        call AdrRectGap(nst%lmn, mpGR)
                        surfrpt=GetSurfOp(nst%surf, nst%axis)
                        call Move2SurfRect(nst%pa, mpGR%cent, mpGR%slen, surfrpt, nst%axis)
                    endif
                else
                    print *, "Inappropriate boundary in MC calculation"
                endif
            else ! leak out
                term=.true.
                exit
            endif
        endif
        nst%rd=INDIR
    enddo
end subroutine




function timecheck() result(time)
!$  use OMP_LIB
    implicit none
    real(8) :: time
!$  if(lomp) then
!$      time=omp_get_wtime()
!$  else
        call cpu_time(time)
!$  endif
end function



end module MCModular

