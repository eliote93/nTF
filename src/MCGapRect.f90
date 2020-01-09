module MCGapRect
    logical, private :: lcmfd
contains

subroutine InitMCGapRect(lcmfd_)
    implicit none
    logical :: lcmfd_
    lcmfd=lcmfd_
end subroutine

function MCRectGap(nst, mpGap) result(term)
    use MCDefine
    use MCBasics
    use CMFDTally
    USE HighOrderSC
    USE Cntl, ONLY : nTracerCntl
    implicit none
    type(nstat), intent(inout) :: nst
    type(STR_RectGap), intent(in) :: mpGap
    type(XsecSet), pointer :: xs
    integer :: surf, leakaxis
    integer :: nadd
    real(8) :: dts, xtot
    logical :: term

    nst%ir=1
    xs=>mpGap%xs(nst%ir)
    nst%dtc=SmpDTC(nst, xs)
    nst%pr=nst%pa-mpGap%cent

    term=.false.
    do while (.NOT. term)
        dts=dtsRect(nst, mpGap%slen)
        if (nst%dtc>dts) then ! leak out to fuel pin-cell surface
            nst%dtc=nst%dtc-dts
            call move(dts, nst)
            call tally(nst, xs, dts, mpGap%tally)
            if (lcmfd) call CMFDRRTally(mpGap%infoCMFD%lmn, nst%g, xs, dts*nst%wt, nst%tid)
            exit ! continue to neighbor node
        else ! collision
            call move(nst%dtc, nst)
            call tally(nst, xs, nst%dtc, mpGap%tally)
            if (lcmfd) call CMFDRRTally(mpGap%infoCMFD%lmn, nst%g, xs, nst%dtc*nst%wt, nst%tid)
            call track_keff(nst, xs)
            term=collision(nst, xs)
            if (term) then
                nst%surf=0
                nst%axis=0
                continue
            endif
            nst%dtc=SmpDTC(nst, xs)
            IF(nTracerCntl%scatod .eq. 0) THEN
              CALL SmpDirIso(nst, nst%dir)
            ELSE
              CALL SmpDirAnIsoGaussian(nst, nst%dir, xs)
              !CALL SmpDirAnIsoGaussian_OLD(nst, nst%dir, xs)
              !CALL SmpDirAnIsoGaussian_OLDmic(nst, nst%dir, xs)
            END IF
            !call SmpDirIso(nst, nst%dir)
            !CALL SmpDirAnIsoGaussian(nst, nst%dir, xs)
        endif
    enddo
    if (mpGap%cut) call ReflCent(mpGap%cutxy, nst)
end function

subroutine Move2SurfaceGapRect(pa, GR, surf, axis)
    use PARAM
    use MCDefine
    implicit none
    real(8), intent(inout) :: pa(3)
    type(STR_RectGap), intent(in) :: GR
    integer, intent(in) :: surf, axis

    if (axis .eq. RADIAL) then
        select  case(surf)
            case(SOUTH)
                pa(2)=GR%cent(2)-GR%slen(2)/2
            case(WEST)
                pa(1)=GR%cent(1)-GR%slen(1)/2
            case(NORTH)
                pa(2)=GR%cent(2)+GR%slen(2)/2
            case(EAST)
                pa(1)=GR%cent(1)+GR%slen(1)/2
        end select
    else ! axial
        select  case(surf)
            case(TOP)
                pa(3)=GR%cent(3)+GR%slen(3)/2
            case(BOTTOM)
                pa(3)=GR%cent(3)-GR%slen(3)/2
        end select
    endif
end subroutine

end module
