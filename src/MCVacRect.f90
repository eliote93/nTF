module MCVacRect

    contains

subroutine MCRectVac(nst, mpVac) 
    use MCDefine
    use MCBasics
    implicit none
    type(nstat), intent(inout) :: nst
    type(STR_RectVac), intent(in) :: mpVac
    type(XsecSet), pointer :: xs
    integer :: surf, leakaxis
    integer :: nadd
    real(8) :: dts, xtot
    logical :: term
    
    nst%pr=nst%pa-mpVac%cent

    dts=dtsRect(nst, mpVac%slen) 
    call move(dts, nst)
    if (mpVac%cut) call ReflCent(mpVac%cutxy, nst)
end subroutine

end module
