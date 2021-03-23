MODULE XSLIB_MOD
    USE PARAM
    IMPLICIT NONE

    type scatmat
        integer :: ib, ie, ioutsb, ioutse
        real, pointer :: from(:)
    end type

    type riftype
        integer :: nid,nrat
        real,pointer,dimension(:) :: rat,ratlog
        real(4),pointer,dimension(:,:,:,:) :: abs,sct,ssct,fis
    end type

    TYPE ppm_TYPE ! photon production matrices
      LOGICAL :: lexist
      INTEGER :: iglow, igup
      INTEGER, POINTER, DIMENSION(:,:) :: ioutpb, ioutpe
      TYPE(SCATMAT), POINTER, DIMENSION(:, :) :: mat
    END TYPE

    TYPE GAMLIBDATA ! Photoatomic Library Data
      INTEGER :: nid
      CHARACTER*8  :: aid
      REAL :: aw
      ! PHOTON REACTION DATA                          !--- JSU EDIT 20170720       FOR ITYP ~= 1 (0, 2)
      !        ONLY FOR ELEMENT **000 (PHOTO-ATOMIC REACTION)
      REAL, POINTER, DIMENSION(:) :: SIGA,SIGTR,SIGS,SIGSTR,KERMA,SIGSS  ! TRANSPORT CORRECTED XS FOR PHOTON TRANSPORT
      TYPE(scatmat), POINTER, DIMENSION(:) :: SM, SMp1, SMp2, SMp3   ! SCATTERING MATRIX UPTO 3rd ORDERS
      REAL, POINTER, DIMENSION(:) :: SIGSP1, SIGSP2, SIGSP3
    END TYPE
    
    TYPE pwxs_type ! point-wise XS file
      SEQUENCE
      INTEGER :: nid, ntemp, nbyte, location
      REAL(4) :: awr, xspot, emin, emax
      LOGICAL :: lfis
      REAL(4), POINTER :: temp(:), xsa(:,:), xss(:,:), xsf(:,:), xst(:,:), tempsq(:)
    END TYPE

    type libdata
      integer :: nid, ityp, ichi, ifis, ibur, inmn
      logical :: lreso=.false.,lfuel=.false.,lclad=.false.
      integer :: ntemp, np1temp, nrtemp, npot
      character(20) :: aid
      real :: aw, mu, kappa, kappa0, dcy, crit_nd, sigp, lamsigp1G
      real, pointer, dimension(:) :: temp, p1temp, rtemp, rtempsq, sig0sq, lamsigp
      real, pointer, dimension(:) :: chi, beta, yield
      real, pointer, dimension(:) :: chip, chid
      real, pointer, dimension(:, :) :: chidk
      real, pointer, dimension(:, :) :: siga, sigf, signf, sigtr, sigs, sigstr, sigss ! no transport corrected self-scattering
      real, pointer, DIMENSION(:, :) :: sigsp1, sigsp2, sigsp3, sigsp4, sigsp5  
      type(scatmat), pointer, dimension(:, :) :: sm, smp1, smp2, smp3, smp4, smp5
      real, pointer, dimension(:) :: sign2n, sign3n
      !!! Resonance Data !!!
      INTEGER :: iresoreg = -1
      ! Effective XSs instead of Resonance Integral
      real(4), pointer, dimension(:, :, :) :: xsalog, xss, ri_a
      real(4), pointer, dimension(:) :: ssf
      ! Subgroup Data
      integer :: nlv,nlvflx
      real, pointer, dimension(:, :) :: lvabs, lvfis, lvflx
      real, pointer, dimension(:, :) :: lvabs_log, lvflx_log
      real, pointer, dimension(:, :, :) :: wgtabs, wgtfis
      ! Resonance Interference Data
      integer :: nsig0,nrif
      integer,pointer,dimension(:) :: rifid
      type(riftype),pointer,dimension(:) :: rif
      integer :: icat = 0
      ! Photon Production Data
      INTEGER :: ipp
      LOGICAL :: lphoton(1:4)
      TYPE(ppm_TYPE), POINTER :: ppm(:)
      ! exact kappa calculation data
      REAL,POINTER :: exkappa(:), exkappa0(:)                                ! 1:FissionProducts, 2:PromptNeutron, 3:DelayedNeutron, 4:PromptGamma, 5: DelayedGamma, 6: DelayedBeta
      REAL :: capkappa, capQ, n2nQ, n3nQ                                     !
      REAL,POINTER,DIMENSION(:,:) :: kerma_t,kerma_s,kerma_d,kerma_p,kerma_f ! Neutron induced KERMA (total, scattering, disappearance(absorption except fission), photon, fission
      ! Photoatomic data (pointer)
      TYPE(GAMLIBDATA), POINTER :: Phatom
      ! Parameters used in Transient Calculation ** Valid only with nver>=101 library             ! EDIT JSU 2020/02/10
      REAL, POINTER, DIMENSION(:) :: dcy_del   ! decay constant of delayed neutron precursor [1/s]
      REAL, POINTER, DIMENSION(:, :) :: InvV   ! Inverse Speed for transient [s/m]  
      ! Parameters for PSM
      integer :: pxsid = 0, icatf, icatm
    end type

    INTEGER :: libtyp
    !INTEGER :: ihel0, ihel1, ihel2, ihel3, ihel4, ihel5, ihel7, ihel8, ihel10
    !Library
    INTEGER :: indreshel, indxsdhel, indburhel, indp1mhel
    INTEGER :: indn2nhel, indn3nhel, n2nhel, n3nhel
    INTEGER :: indchxhel, nchixhel, indp2mhel, indp3mhel, indsghel
    !Library Index
    INTEGER :: noghel, nofghel, norghel, notghel, nchihel
    !         #Group, #fast Group(above Resonance), # Resonance Grp,
    !         # of grp with up-scat, # of group with fission source
    INTEGER :: nelthel, nelrhel, nreshel, nburhel, nbahel, nfishel, np1hel
    !          # of element, # isotope with complete CX(P0 scat), # of istope with RI
    !          # of burn up isotop, # of burnable - absorber isotopes, # fissionable isotope, # isotope with P1
    INTEGER :: nxsghel,nwsghel,nwasgfhel,moghel, ntcateg, nyldhel, nflxhel, nsubghel
    !
    INTEGER :: nmlver = 0, nrlver = 0, nslver = 0
    ! xs library (helios) data variables
    INTEGER, POINTER, DIMENSION(:) :: idcodehel,nxsthel,                                         &
                                      ntemphel,npothel,ntabhel,idnp1hel,                         &
                                      np1thel,nwaddhel,idburhel,nuclidhel,idfishel,              &
                                      idn2nhel,idn3nhel

    INTEGER, POINTER, DIMENSION(:, :) :: istarthel, ip1stahel, ip2stahel, ip3stahel, infchxhel

    REAL, POINTER, DIMENSION(:) :: chihel, chix

    REAL(4), POINTER, DIMENSION(:) :: xstemphel, rstab1hel, p1temphel, uhel, enbhel,             &
                                      awhel, xenhel, xen0hel, decayhel, enbgam

    REAL(4), POINTER, DIMENSION(:, :) :: xsdatahel,  resdathel, rstab2hel,                          &
                                      p1datahel, p2datahel, p3datahel,                           &
                                      xsghel, xsglog, wsghel, xasgfhel, xasgflog, wasgfhel,      &
                                      chixhel, betadnhel, yieldhel, xsn2nhel, xsn3nhel

    ! mapping of nuclide ID to helios library order
    INTEGER, PARAMETER :: nxtempmap=5000
    INTEGER, POINTER :: itempmap(:, :), itempmapp1(:, :), mapn2n(:), mapn3n(:)
    INTEGER,PARAMETER :: maxnid=110000,nmaxz=100
    INTEGER :: mapnucl(maxnid),  mapnuclp13(maxnid), mapfis(maxnid)
    INTEGER, POINTER :: mapnuclRes(:,:)
    !! Variables for Gamma Calculation  |- JSU EDIT 20190509
    INTEGER :: nogghel, nelmGAM, noggphl
    INTEGER :: mapnuclpp(maxnid,1:4)  ! Mapping nuclides for Photon Production...
    INTEGER :: nphprod(1:4) ! PHoton PRODuction
    INTEGER, POINTER, DIMENSION(:,:) :: idphprod
    !INTEGER, POINTER, DIMENSION(:) :: nuclidgam,idfisGAM,idradGAM,idinelGAM,idelmGAM
    INTEGER, POINTER, DIMENSION(:) :: elementGAM,iso2elm
    INTEGER :: mapEleGAM(maxnid)
    ! indices of U-238 for resonance cal.
    !INTEGER :: it1r, it2r, jt1r, jt2r, kt1r
    ! library data by isotope
    TYPE(libdata), POINTER, dimension(:) :: ldiso
    TYPE(GAMLIBDATA), POINTER, DIMENSION(:) :: phatom
    ! Variables for PSM
    INTEGER :: npwxs
    TYPE(pwxs_type), POINTER :: pwxs(:)
    !
    
    type mlgtype ! Macro-Level Grid type
    ! The data for Macro-level grid of each plane 'iz'
    !~ Since normal MLG is done for whole plane, max and min radii are standard of MLG
    ! Fuel : f_nmaclv : # of MLG   / f_maclv : MLG   / f_maclv_log : MLG_log
    ! AIC  : f_nmaclv1G : # of MLG / f_maclv1G : MLG / f_maclv1G_log : MLG_log
    ! Clad : c_nmaclv1G : # of MLG / c_maclv1G : MLG / c_maclv1G_log : MLG_log => Identical in whole planes
        SEQUENCE
        integer :: f_nmaclv,f_nmaclv1G,c_nmaclv1G ! # of macro level grid for fuel and AIC(f~1G) (DEF:8) and cladding(DEF:5)
        real,pointer,dimension(:) :: f_maclv,f_maclv1G,c_maclv1G
        real,pointer,dimension(:) :: f_maclv_log,f_maclv1G_log,c_maclv1G_log
        ! EDIT JSU 2020/08/12 Unused variables
!       INTEGER :: c_nmaclv        
!        REAL,POINTER,DIMENSION(:) :: c_maclv, c_maclv_log
    ! Variables for Equivalent Dancoff Cell
    !~ Since EDC is done for individual pin, radius of pins themselves are standard of MLG
    ! ngeomtype : # of fuel pellets in the plane
    ! f_maclv_pin(_log) : MLG(_log) for eqch fuell pellets (nMLG, nGeomtype)
        integer :: ngeomtype
        real,pointer,dimension(:,:) :: f_maclv_pin       ! 1st index : macro level, 2nd index : geom type
        real,pointer,dimension(:,:) :: f_maclv_pin_log   ! 1st index : macro level, 2nd index : geom type
        REAL :: del_c1g, del_fmg, del_f1g
    end type
    TYPE(mlgtype),POINTER :: mlgdata(:) ! Allocated and prepared in PrepareXS.F90>SetMLGrid
    TYPE(mlgtype) :: mlgdata0 ! variable for cladding (unchanged..)
    
    INTEGER :: igresb, igrese
    ! temporary variable
    REAL(4), POINTER :: chitem(:)
    REAL(8), POINTER :: sigbeta(:)
    ! Positions of I,  Xe,  Pm in Burnable Isotopes
    INTEGER :: iphel_i, iphel_xe, iphel_pm

    CHARACTER(120) :: LibInfo(3)

    INTEGER, POINTER :: IDRES_COR(:,:),NRES_COR(:)
    LOGICAL :: lResCore(10000) ! not used
    INTEGER :: nlvmax = 0
    INTEGER,PARAMETER :: nlvflxmax=30
    
    integer :: nCat
    integer,pointer :: nActiveCat(:)
    type cattype
        integer :: repid,niso
        integer,pointer :: idiso(:)
    end type
    type(cattype),pointer :: ResoCat(:)
    logical,pointer :: ResoCatUse(:,:)

    CONTAINS

    SUBROUTINE CoreResIsoUpdate(idiso, niso, iz)
        !Update the lists of resonance isotopes
        IMPLICIT NONE
        INTEGER :: niso, idiso(niso), iz
        INTEGER :: i,  j, id, id2, xxx
        LOGICAL :: flag
        do j=1,niso
          id=idiso(j)
          id2=mapnucl(id)
          IF(.not.ldiso(id2)%lreso) CYCLE
          xxx = mod(id,1000)
          flag = .false.
          if (xxx.gt.500) then
            flag = .true.
            id = id - 500
          endif
          DO i = 1, nreshel
            if (id .eq. idres_cor(i,iz)) exit
          ENDDO
          if (i.gt.nreshel) then
            nRes_Cor(iz) = nRes_Cor(iz) + 1
            idres_cor(nRes_Cor(iz),iz) = id
            mapnuclRes(id,iz) = nRes_Cor(iz)
            mapnuclRes(id+500,iz) = nRes_Cor(iz)
          endif
        enddo
    END SUBROUTINE

END MODULE
