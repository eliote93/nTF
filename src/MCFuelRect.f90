module MCFuelRect
    logical, private :: lcmfd, lxslib, lmicxs
contains

subroutine InitMCFuelRect(lcmfd_, lxslib_, lmicXs_)
    implicit none
    logical :: lcmfd_, lxslib_, lmicxs_
    lcmfd=lcmfd_
    lxslib=lxslib_
    lmicxs=lmicxs_
end subroutine

function MCRectFuel(nst, mpFR) result(term)
    use MCDefine
    use MCBasics
    use CMFDTally
    USE MCTally, ONLY : SpecTally
    !USE MCModular !, ONLY : FXR, ctrlMC
    USE TYPEDEF, ONLY : FxrInfo_Type
    implicit none
    type(nstat), intent(inout) :: nst
    type(STR_RectFuel), intent(in) :: mpFR
    type(XsecSet), pointer :: xs
!    integer :: surf, leakaxis
    integer :: nadd, naddsum
    real(8) :: dts, xtot, DTSR, DTSC
    logical :: term, lrect
    INTEGER :: ifxr, iz, surf, axis, io, ct
    TYPE(FxrInfo_Type), POINTER :: myFXR

    naddsum=0

    nst%pr=nst%pa-mpFR%cent
    if (nst%axis .eq. AXIAL) then
        nst%ir=findring(mpFR%nr, mpFR%rrsq, nst%pr(1), nst%pr(2))
    elseif(nst%axis .eq. RADIAL) then
        IF( mpFR%nor .EQ. 0 )THEN
            nst%ir=0
        ELSE  ! spacer grid modeling
            nst%ir=findring(mpFR%nr, mpFR%rrsq, nst%pr(1), nst%pr(2))
        ENDIF
    endif

    xs=>mpFR%xs(nst%ir)
    nst%dtc=SmpDTC(nst, xs)
    !io=109
    !OPEN(unit=io, file='mc.dbg', status='replace')

    ct=0
    term=.false.
    do while (.NOT. term)
        !ct=ct+1;
        !write(*,*) ct ! dummy counter
        ! BYS edit 16/12/06 r557 for spacer grid
        if( (mpFR%nor .EQ. 0) .OR. (nst%ir.GT.mpFR%nor) )THEN
        if (nst%rd .eq. INDIR) then
            if (nst%ir+1<=mpFR%nr) then
                dts=dtcl(nst%pr, nst%dir, mpFR%rr(nst%ir+1), mpFR%slen(3), nst%rd, nst%surf, nst%axis)
            else
                nst%rd=OUTDIR
                cycle
            endif

            if (dts .eq. 0.) then
                nst%rd=OUTDIR
                cycle
            else
                if (nst%dtc>dts) then ! move to next region
                    nst%dtc=nst%dtc-dts
                    call MoveWithinCell()
                    if (nst%axis .eq. RADIAL) then
                        nst%ir=nst%ir+1
                        xs=>mpFR%xs(nst%ir)
                        nst%dtc=SmpDTC(nst, xs)
                    else
                        exit
                    endif
                else    ! collision
                    call proc_collision()
                endif
            endif
        else ! OUT DIR
            if (nst%ir .eq. 0) then
                dts=dtsRect(nst, mpFR%slen)
                if (nst%dtc>dts) then ! leak out to fuel pin-cell surface
                    nst%dtc=nst%dtc-dts
                    call MoveWithinCell()
                    exit ! continue to neighbor node
                else ! collision
                    call proc_collision()
                endif
            else ! leak out to annular geometry
                dts=dtcl(nst%pr, nst%dir, mpFR%rr(nst%ir), mpFR%slen(3), nst%rd, nst%surf, nst%axis)
                if (nst%dtc>dts) then
                    nst%dtc=nst%dtc-dts
                    call MoveWithinCell()
                    if (nst%axis .eq. RADIAL) then
                        nst%ir=nst%ir-1
                        xs=>mpFR%xs(nst%ir)
                        nst%dtc=SmpDTC(nst, xs)
                    else
                        exit ! continue to axial neighbor node
                    endif
                else    ! collision
                    call proc_collision()
                endif
            endif
        endif !end or normal circle
        else ! start of outer circle\

            DTSR=dtsRect(nst, mpFR%slen)

            surf=nst%surf
            axis=nst%axis
            lrect=0
            IF( nst%ir .EQ. 0 )THEN ! outer most ring
                nst%rd=INDIR
                DTSC=dtcl(nst%pr, nst%dir, mpFR%rr(+1), mpFR%slen(3), nst%rd, nst%surf, nst%axis)
                IF( DTSC .EQ. 0 )THEN ! OUTDIR
                    dts=DTSR
                    lrect=1
                ELSE ! INDIR
                    IF( DTSC .GT. DTSR )THEN
                        dts=DTSR
                        lrect=1
                    ELSE
                        dts=DTSC
                    ENDIF
                ENDIF
            ELSE !(nst%ir .NE. 0)THEN normal ring#1~nor

                DTSC=dtcl(nst%pr, nst%dir, mpFR%rr(nst%ir+1), mpFR%slen(3), nst%rd, nst%surf, nst%axis)
                IF( DTSC.EQ.0 )THEN !OUTDIR
                    nst%rd=OUTDIR
                    DTSC=dtcl(nst%pr, nst%dir, mpFR%rr(nst%ir), mpFR%slen(3), nst%rd, nst%surf, nst%axis)
                    IF( DTSC .GT. DTSR )THEN
                        dts=DTSR
                        lrect=1
                    ELSE
                        dts=DTSC
                    ENDIF
                ELSE ! INDIR
                    nst%rd=INDIR
                    IF( DTSC .GT. DTSR )THEN
                        dts=DTSR
                        lrect=1
                    ELSE
                        dts=DTSC
                    ENDIF
                ENDIF
            ENDIF
            IF( lrect )THEN ! headed to rectangular boundary
                IF( nst%dtc .GT. dts )THEN ! migration
                    nst%dtc=nst%dtc-dts
                    nst%surf=surf
                    nst%axis=axis
                    call MoveWithinCell()
                    exit ! continue to neighbor node
                else ! collision
                    call proc_collision()
                endif
            ELSE ! move to other circle
                if (nst%dtc>dts) then ! move to next region
                    nst%dtc=nst%dtc-dts
                    call MoveWithinCell()
                    if (nst%axis .eq. RADIAL) then
                        IF( nst%rd .EQ. INDIR )THEN
                            nst%ir=nst%ir+1
                        ELSE ! ( nst%rd .EQ. OUTDIR )
                            nst%ir=nst%ir-1
                        ENDIF
                        xs=>mpFR%xs(nst%ir)
                        nst%dtc=SmpDTC(nst, xs)
                    else
                        exit
                    endif
                else    ! collision
                    call proc_collision()
                endif
            ENDIF
        endif
        !WRITE(io,'(100E15.7)') nst%dir
        !WRITE(io,'(100E15.7)') nst%pr
        !continue


    enddo
    if (lcmfd .and. naddsum>0)  call SrcDistTally(mpFR%infoCMFD%lmn, naddsum, nst%tid)
    if (mpFR%cut) then
        call ReflCent(mpFR%cutxy, nst)
!        nst%pa=nst%pr+mpFR%cent
    endif
    !close(io)
    continue

contains

subroutine MoveWithinCell
    call move(dts, nst)
!    call trackFSD(mpFR%lmn, xs.xsnf(nst%g), dts*nst%wt, nst%tid)
    call trackFSD(mpFR%lmn, nst%g, xs, dts*nst%wt, nst%tid)
    call tally(nst, xs, dts, mpFR%tally)
    if (lcmfd) call CMFDRRTally(mpFR%infoCMFD%lmn, nst%g, xs, dts*nst%wt, nst%tid)
    if (lxslib .AND. .NOT. lmicXs) call SpecTally(mpFR%infoCMFD%lmn, nst%g, xs, dts*nst%wt, nst%tid)
end subroutine

subroutine proc_collision
    USE HighOrderSC
    !USE Cntl, ONLY : nTracerCntl
    call move(nst%dtc, nst)
!    call trackFSD(mpFR%lmn, xs.xsnf(nst%g), nst%dtc*nst%wt, nst%tid)
    call trackFSD(mpFR%lmn, nst%g, xs, nst%dtc*nst%wt, nst%tid)
    call tally(nst, xs, nst%dtc, mpFR%tally)
    if (lcmfd) call CMFDRRTally(mpFR%infoCMFD%lmn, nst%g, xs, nst%dtc*nst%wt, nst%tid)
    if (lxslib .AND. .NOT. lmicXs) call SpecTally(mpFR%infoCMFD%lmn, nst%g, xs, nst%dtc*nst%wt, nst%tid)
    call track_keff(nst, xs)
    nadd=SmpNeutron(nst, xs)
    IF( lxslib .AND. lmicXs )THEN
        call BankNeutronIso(nst, xs, nadd)
    ELSE
        call BankNeutron(nst, nadd)
    ENDIF



    naddsum=naddsum+nadd
    term=collision(nst, xs)
    if (term) then
        nst%surf=0
        nst%axis=0
        !write(*,*) 'die'
        continue
    else
        nst%rd=INDIR
    endif
    nst%dtc=SmpDTC(nst, xs)
    IF(nTracerCntl%scatod .eq. 0) THEN
      CALL SmpDirIso(nst, nst%dir)
    ELSE
      CALL SmpDirAnIsoGaussian(nst, nst%dir, xs)
      !CALL SmpDirAnIsoGaussian_OLD(nst, nst%dir, xs)
      !CALL SmpDirAnIsoGaussian_OLDmic(nst, nst%dir, xs)
    ENDIF
    !write(*,*) 'scat'
    !call SmpDirIso(nst, nst%dir)
!    CALL SmpDirAnIsoGaussian(nst, nst%dir, xs)
end subroutine

end function

subroutine Move2SurfaceFuelRect(pa, FR, surf, axis)
    use PARAM
    use MCDefine
    implicit none
    real(8), intent(inout) :: pa(3)
    type(STR_RectFuel), intent(in) :: FR
    integer, intent(in) :: surf, axis

    if (axis .eq. RADIAL) then
        select case(surf)
            case(SOUTH)
                pa(2)=FR%cent(2)-FR%slen(2)/2
            case(WEST)
                pa(1)=FR%cent(1)-FR%slen(1)/2
            case(NORTH)
                pa(2)=FR%cent(2)+FR%slen(2)/2
            case(EAST)
                pa(1)=FR%cent(1)+FR%slen(1)/2
        end select
    else ! axial
        select case(surf)
            case(TOP)
                pa(3)=FR%cent(3)+FR%slen(3)/2
            case(BOTTOM)
                pa(3)=FR%cent(3)-FR%slen(3)/2
        end select
    endif
end subroutine

end module
