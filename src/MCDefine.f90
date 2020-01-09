module MCDefine
    use PARAM
    integer, parameter :: RPT=-2, VAC=-1, REFL=0, FRECT=1, GRECT=2, VRECT=3
    integer, parameter :: CYLB=1, CYLT=2, CYLR=3
    integer, parameter :: RADIAL=1, AXIAL=2
    integer, parameter :: INDIR=1, OUTDIR=2, OUTPLANE=3
    integer, parameter :: SurfOp(4,2)=RESHAPE((/NORTH, EAST, SOUTH, WEST, TOP, BOTTOM, -1, -1/), SHAPE(SurfOp))
    
    type MCCTRL
        integer :: nxyz(3), nxyzCMFD(3)
        integer :: ngMC, ngCMFD, gThrm
        integer :: ninact, nact, nht
        logical :: lomp=.false., lcmfd=.false., lfdb(2)=.false.
        integer :: nthrd=1
        integer, pointer, dimension(:,:) :: CoreMap
        logical, pointer, dimension(:,:,:) :: cutxy
        logical, pointer, dimension(:,:) :: lwrt
        logical, pointer, dimension(:,:) :: lgap
        integer :: lwrtx(2), lwrty(2)
        LOGICAL :: lXsLib, lsgfsp, lfirstcycle, lmicxs
        
    end type
    
    type CMFDTallyHelper
        integer :: icmfd !, nsurf(2)
        integer :: lmn(3)
        logical :: lsurf(4,2)
    end type
    
    type MCTallySet
        real(8), pointer, dimension(:,:) :: trkl
        real(8) :: pinpow
    end type
    
    type XsecSet
        integer :: idx, ifxr, iz, ixs
        real(8), pointer, dimension(:) :: xst, xstr, xsnf, xskf, xsa, chi, xss, xsrmv
        LOGICAL :: lfuel
        INTEGER :: niso
        INTEGER, POINTER :: idiso(:)
        REAL(8), POINTER, DIMENSION(:) :: pnum        
        real(8), pointer, dimension(:,:) :: FisCDF
        real(8), pointer, dimension(:,:) :: xssm, xssm2g
        real(8), pointer, dimension(:,:) :: SctCDF
        real(8), pointer, dimension(:) :: xsaa ! for fixing negative scattering cross section
        REAL(8), POINTER, DIMENSION(:,:) :: xss0, xss1, xss2, xss3
        REAL(8), POINTER, DIMENSION(:,:) :: rt1, rt2, w1, w2
    end type
    
    type STR_RectFuel
        integer :: idx, idxt, lmn(3)
        integer :: FxrIdxSt
        integer :: nr, nor
        real(8) :: slen(3), cent(3), vol
        real(8), pointer, dimension(:) :: rr
        real(8), pointer, dimension(:) :: rrsq
        type(XsecSet), pointer, dimension(:) :: xs
        type(MCTallySet) :: tally
        type(CMFDTallyHelper) :: infoCMFD
        logical :: cut, cutxy(2)
    end type
    
    type STR_RectGap
        integer :: idx, idxt, lmn(3)
        integer :: FxrIdxSt
        real(8) :: slen(3), cent(3), vol
        type(XsecSet), pointer, dimension(:) :: xs
        type(MCTallySet) :: tally
        type(CMFDTallyHelper) :: infoCMFD
        logical :: cut, cutxy(2)
    end type

    type STR_RectVac
        integer :: idx, lmn(3)
        integer :: FxrIdxSt
        real(8) :: slen(3), cent(3)
        logical :: cut, cutxy(2)
    end type

    
    type cellmap_elmt
        integer :: telmt
        integer :: ielmt
    end type cellmap_elmt
    
    type nstat
        real(8) :: pr(3), pa(3), dir(3) ! relative position, absolute position, direction
        integer :: g, lmn(3), ir, idx ! group, cell index for xyz, index in geom
        integer :: rd ! OUTDIR | INDIR
        real(8) :: dtc, dts, wt
        integer :: surf, axis ! leak out info
        integer(8) :: seed
        integer :: tid
        INTEGER :: iso ! reacted isotope in previous cycle
    end type nstat
    
!    integer, parameter :: FLG_NONE=0, FLG_LEAKAXIAL=1
    type nfs
        integer :: lmn(3), idx, ir
        real(8) :: pr(3), wt
        integer :: iseq ! for debugging
        INTEGER :: iso
    end type nfs
contains 

subroutine CopyNFS(dest, org)
    type(nfs), intent(out) :: dest
    type(nfs), intent(in) :: org
    dest%lmn=org%lmn
    dest%idx=org%idx
    dest%ir=org%ir
    dest%pr=org%pr
    dest%wt=org%wt
end subroutine

end module MCDefine
