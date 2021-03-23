#include <defines.h>
MODULE PointXSRT_MOD
USE PARAM
IMPLICIT NONE
!
CHARACTER*256 :: pwxsdir ! directory file of point-wise XS .pdl file
INTEGER :: nfbyte
INTEGER :: neg, ipsst, ipsed
REAL    :: psm_emin, psm_emax, del1pg, invdel1pg, Esdmin, Esdmax, invEsdmax
REAL(4), ALLOCATABLE, DIMENSION(:) :: geg4, xsa4, xss4, xsf4, xst4
REAL   , ALLOCATABLE, DIMENSION(:) :: geg, gug
!
REAL, PARAMETER :: E0 = 1.E+7
!
INTEGER, ALLOCATABLE, DIMENSION(:) :: icgb, icge, ipw2mg
REAL, ALLOCATABLE, DIMENSION(:,:) :: delu_pw
!
! First-flight collision probability of isolated pins
REAL :: PGrid_max=5.E+2, PGrid_min=1.E-4, delXS, invdelXS
INTEGER :: nPGrid = 5000
REAL, ALLOCATABLE, DIMENSION(:) :: Pij_Grid
!
INTEGER :: nrmax = 1! maximum number of fuel rings
INTEGER :: npin_iso ! number of isolated pincell
TYPE PSM_IsoPin_TYPE ! Isolated pin information
  REAL, POINTER :: Pji_iso(:,:,:) ! for fast calculation, j <- i
END TYPE
TYPE(PSM_ISOPIN_TYPE), POINTER :: PSM_ISOPIN(:)
! Slowing down variable
!
INTEGER, POINTER  :: mapnuclPXS(:,:), IDPXS_COR(:,:), NPXS_COR(:)
!
!
! TYPES for PSM
TYPE PSMMacXs_Type
  !SEQUENCE
  INTEGER :: id = 0
  logical :: lalloc = .FALSE.
  logical :: lfuel
  INTEGER :: ng = 0, nfxr = 0
  INTEGER :: niso = 0
  Logical :: lIsoAlloc = .FALSE.
  LOGICAL :: lallocIsoDataSLD = .FALSE.
  INTEGER :: nrmax = 0, neg = 0
  !REAL :: mavgm(3), alpham(3), invsubalpham(3) ! average mass of moderator
  ! XsMac (ig, ifxr)
  REAL, POINTER :: xsmact(:,:)               !Total Xs
  REAL, POINTER :: xsmactr(:,:)              !Transport
  REAL, POINTER :: xsmaca(:,:)               !Absorption
  REAL, POINTER :: xsmacf(:,:)               !Fission
  REAL, POINTER :: xsmacnf(:,:)              !nu-Fission
  REAL, POINTER :: xsmackf(:,:)              !kappa-Fission
  REAL, POINTER :: xsmacs(:,:)               !Total Scattering Xs
  ! IsoXsMac (ig, iiso, ifxr)
  REAL, POINTER, DIMENSION(:,:,:) :: IsoMGXSMacA, IsoMGXsMacTr, IsoMGXsMacT !BYS edit 14/06/01
  REAL, POINTER, DIMENSION(:,:,:) :: IsoMGXSMacS, IsoMGXsMacf,  IsoMGXSRadCap, IsoMGXSMacNF, IsoMGMacNu
  ! IsoXSMac Slowing Down Data
  REAL, POINTER :: phip(:,:), phip_(:,:), phim(:), SSrcal(:,:,:), xss(:,:,:), xss0(:,:), xst(:,:), Qf(:), xsf(:,:), xsa(:,:)
  REAL, POINTER, DIMENSION(:) :: Pim, Pei_iso, Pmi
  LOGICAL, POINTER :: lfsrc(:,:)
  REAL, POINTER, DIMENSION(:,:) :: Pji
  REAL, POINTER :: favgm(:,:), alphaf(:,:), invsubalphaf(:,:),xspf(:,:),xspf_inval(:,:) ! average mass and inv(1-alpha) of fuel
END TYPE

  CONTAINS
  
! Subroutine to read PXS file
SUBROUTINE ReadPWXS(indev, pdlfile)
USE ioutil,         only : terminate, toupper, openfile
USE XSLIB_MOD,       ONLY : libdata, ldiso, nelthel, pwxs_type, pwxs, npwxs
IMPLICIT NONE
INTEGER :: indev
INTEGER :: ixs, i, j, k, it, ig, nid
CHARACTER*256 :: pdlfile 
REAL(4) :: tsp1, tsp2, tsp3
TYPE(pwxs_type), POINTER :: txs
TYPE(libdata),POINTER :: lib,ljb
INTEGER :: xxx, njd, ir, mlbnid

pwxsdir = pdlfile
close(indev)
call openfile(indev,.TRUE.,.TRUE.,.FALSE.,pdlfile)
READ(indev) npwxs
READ(indev) nfbyte
READ(indev) neg
READ(indev) psm_emin, psm_emax
ALLOCATE(geg(0:neg),gug(0:neg),geg4(neg))
READ(indev) geg(1:neg)
!geg(1:neg) = geg4
geg(0) = geg(1) + (geg(1)-geg(2))
IF (psm_emax .NE. geg(1)) psm_emax = geg(1)
IF (psm_emin .NE. geg(neg)) psm_emin = geg(neg)
!DEALLOCATE(geg4)
DO i = 0, neg
  gug(i) = LOG(E0/geg(i))
END DO
ALLOCATE(pwxs(npwxs))
!
DO i = 1, npwxs
  READ(indev) pwxs(i)%nid
  READ(indev) pwxs(i)%ntemp
  ALLOCATE(pwxs(i)%temp(pwxs(i)%ntemp))
  ALLOCATE(pwxs(i)%tempsq(pwxs(i)%ntemp))
  READ(indev) pwxs(i)%temp
  pwxs(i)%tempsq = SQRT(pwxs(i)%temp)
  READ(indev) pwxs(i)%lfis
  READ(indev) pwxs(i)%xspot
  READ(indev) pwxs(i)%awr
  READ(indev) pwxs(i)%nbyte
  READ(indev) pwxs(i)%location
END DO
!
DO j=1,npwxs
  DO i=1,nelthel
      mlbnid=ldiso(i)%nid
      xxx=mod(mlbnid,1000)
      IF (xxx.gt.500) mlbnid=mlbnid-500
      IF (pwxs(j)%nid.eq.mlbnid) THEN
        ldiso(i)%pxsid =  j
      END IF
  END DO
END DO
DO i = 1, npwxs
  ALLOCATE(pwxs(i)%xsa(neg,pwxs(i)%ntemp),pwxs(i)%xss(neg,pwxs(i)%ntemp),pwxs(i)%xst(neg,pwxs(i)%ntemp))
  READ(indev) pwxs(i)%xsa
  READ(indev) pwxs(i)%xss
  pwxs(i)%xst = pwxs(i)%xsa + pwxs(i)%xss
  IF (pwxs(i)%lfis) THEN
    ALLOCATE(pwxs(i)%xsf(neg,pwxs(i)%ntemp))
    READ(indev) pwxs(i)%xsf
  END IF
END DO

END SUBROUTINE

! Subroutine to set Pij table
SUBROUTINE CalcPSM_ISOPIN
USE PARAM
USE GEOM, ONLY : nCellTYPE0, Core
USE RAYS,             ONLY : RayInfo
USE Core_mod,         ONLY : CMInfo, FmInfo, THInfo, GroupInfo
USE PE_MOD,           ONLY : PE
USE CNTL,             ONLY : nTracerCntl
USE TYPEDEF, ONLY : Cell_TYPE, basicgeom
USE CP_mod, ONLY : CalcPijGam_
IMPLICIT NONE
TYPE(Cell_Type),POINTER :: CellInfo(:)
INTEGER :: icel, celstr, jcel, ig
LOGICAL :: lsame
INTEGER, ALLOCATABLE :: ipin_iso2cel(:)
REAL, POINTER, DIMENSION(:) :: vol, gam, quartersurfvol,sigtvol,sigrvol, invsigt, sigt, invsigtvol
REAL, POINTER, DIMENSION(:,:) :: Pij
REAL :: sumvol, quartersurf
INTEGER :: nr, i, j
!
CellInfo => Core%CellInfo
ALLOCATE(ipin_iso2cel(nCellTYPE0))
!
! Count Base geom Isolated cells
DO icel = 1, nCellType0
  IF (CellInfo(icel)%lfuel .AND. .NOT.CellInfo(icel)%geom%lrect) EXIT
END DO
!
IF (icel .GT. nCellTYPE0) RETURN ! for the case there is no fuel in all the cell 
celstr = icel
npin_iso = 1
CellInfo(icel)%icelPSM = 1
CellInfo(icel)%lPSMcel = .TRUE.
ipin_iso2cel(npin_iso) = icel
!
DO icel = celstr+1, nCellType0
  IF (.NOT.CellInfo(icel)%lfuel .OR. CellInfo(icel)%geom%lrect) CYCLE
  CellInfo(icel)%lPSMcel = .TRUE.
  !--start searching
  DO jcel = 1, npin_iso
    CALL ChkSameFuelStr(CellInfo(icel), CellInfo(ipin_iso2cel(jcel)), lsame)
    IF (lSame) EXIT
  END DO
  IF (lSame) THEN
    CellInfo(iCel)%icelPSM = jcel
  ELSE
    npin_iso = npin_iso + 1
    CellInfo(iCel)%icelPSM = npin_iso
    ipin_iso2cel(npin_iso) = icel
  ENDIF
ENDDO
!
! Set the first-flight collision probability of isolated pins
ALLOCATE(Pij_Grid(nPGrid))
!delXS = (PGrid_max-PGrid_min)/nPGrid
!invdelXS = 1._8/delXS
!Pij_grid = PGrid_min
!DO i = 2, nPGrid
!  Pij_Grid(i) = Pij_Grid(i-1) + delXS
!END DO
delXS = LOG(PGrid_max/PGrid_min)/nPGrid
invdelXS = 1._8/delXS
Pij_grid(1) = PGrid_min
DO i = 2, nPGrid
  Pij_Grid(i) = Pij_Grid(i-1) * EXP(delXS)
END DO
!
ALLOCATE(PSM_ISOPIN(npin_iso))
DO jcel = 1, npin_iso
  icel = ipin_iso2cel(jcel)
  nr = CellInfo(icel)%nfueldiv
  nrmax = MAX(nr, nrmax)
  ALLOCATE(PSM_ISOPIN(jcel)%Pji_iso(nr,nr,nPGrid))
  ALLOCATE(vol(nr),sigtvol(nr),sigt(nr),invsigt(nr),sigrvol(nr),invsigtvol(nr),gam(nr), quartersurfvol(NR))
  ALLOCATE(Pij(nr,nr))
  
  vol(1)=pi*CellInfo(icel)%Rad_CP(1)*CellInfo(icel)%Rad_CP(1)
  sumvol=vol(1)
  DO i=2,nr
      vol(i)=pi*CellInfo(icel)%Rad_CP(i)*CellInfo(icel)%Rad_CP(i)-sumvol
      sumvol=sumvol+vol(i)
  ENDDO 
  quartersurf=hpi*CellInfo(icel)%Rad_CP(nr)
  
  DO ig = 1, nPGrid  
    sigt = Pij_Grid(ig)
    invsigt=1._8/sigt
    sigtvol=vol*sigt
    invsigtvol=1._8/sigtvol
    sigrvol=sigtvol   ! zero scattering XS
    DO i=1,nr
        quartersurfvol(i)=vol(i)*quartersurf
    ENDDO    
    call CalcPijGam_(Pij,gam,CellInfo(icel)%Rad_CP(1:nr),sigt,sigtvol,nr)
    DO i = 1, nr
      Pij(:,i) = Pij(:,i) * invsigtvol ! convert to first-flight CP
    END DO
    DO i = 1, nr
      DO j = 1, nr
        PSM_ISOPIN(jcel)%Pji_iso(j,i,ig) = Pij(i,j)
      END DO
    END DO 
  END DO
  DEALLOCATE(vol,sigtvol,sigt,invsigt,sigrvol,invsigtvol,quartersurfvol, gam)
  DEALLOCATE(Pij)
END DO
!
!
DEALLOCATE(ipin_iso2cel)
!
!nTRACERCntl%lPSMISO = .TRUE.
!
!
!CALL CalcDancoff(Core, FmInfo%Fxr, RayInfo, THInfo, nTracerCntl, PE)
!CALL SetMGnPWInfo

END SUBROUTINE

! SUBROUTINE TO VERIGY ONLY IF BOTH THE CELL HAVE THE SAME FUEL RADIUS WITH FXR...
SUBROUTINE ChkSameFuelStr(iCellinfo, jCellinfo, lsame)
USE TYPEDEF, ONLY : Cell_TYPE, basicgeom
IMPLICIT NONE
TYPE(Cell_Type) :: iCellInfo, jCellInfo
LOGICAL :: lsame
TYPE(basicgeom) :: iGeom, jGeom

INTEGER :: i, j

iGeom=iCellInfo%geom
jGeom=jCellInfo%geom

lsame=.TRUE.
IF( igeom%ncircle .NE. jgeom%ncircle )THEN
    lsame=.FALSE.
    RETURN
ENDIF
IF ( iCellInfo%ibfuel .NE. jCellInfo%ibfuel) THEN
  lsame=.FALSE.
  RETURN
END IF
IF ( iCellInfo%iefuel .NE. jCellInfo%iefuel) THEN
  lsame=.FALSE.
  RETURN
END IF
IF ( iCellInfo%FuelRad0 .NE. jCellInfo%FuelRad0) THEN
  lsame=.FALSE.
  RETURN
END IF
IF ( iCellInfo%nfueldiv .NE. jCellInfo%nfueldiv) THEN
  lsame=.FALSE.
  RETURN
END IF

DO i=1, iCellInfo%nfueldiv
  IF( abs(iCellInfo%Rad_CP(i) - jCellInfo%Rad_CP(i)) .GT. epsm10 )THEN
      lsame=.FALSE.
      RETURN
  ENDIF
ENDDO

ENDSUBROUTINE

! SUBROUTINE to set the infomation of group condensing...
SUBROUTINE SetMGnPWInfo
USE Core_mod,   ONLY : GroupInfo
USE XSLIB_MOD,  ONLY : noghel, enbhel
IMPLICIT NONE
INTEGER :: ng , ig, ip
!
ng = GroupInfo%norg
ALLOCATE(icgb(GroupInfo%iresogrp1:GroupInfo%iresogrp2), icge(GroupInfo%iresogrp1:GroupInfo%iresogrp2))
ALLOCATE(ipw2mg(neg))
ALLOCATE(delu_pw(neg,GroupInfo%iresogrp1:GroupInfo%iresogrp2))
delu_pw = 0._8
!
del1pg = LOG(geg(1)/geg(2))
invdel1pg = 1/del1pg
!
ig = GroupInfo%iresogrp1
ip = 1
DO WHILE (geg(ip).GT.enbhel(ig))
  ip = ip + 1
END DO
!
icgb(ig) = ip
DO WHILE (geg(ip).GT.enbhel(ig+1))
  ipw2mg(ip) = ig
  ip = ip + 1
END DO
icge(ig) = ip - 1
!
DO ig = GroupInfo%iresogrp1+1, GroupInfo%iresogrp2
  icgb(ig) = ip
  DO WHILE ((enbhel(ig)-geg(ip))*(enbhel(ig+1)-geg(ip)) .LE. 0.)
    ipw2mg(ip) = ig
    ip = ip + 1
    IF (ip .GT. neg) GOTO 10
  END DO
  icge(ig) = ip - 1
END DO
!
!DO ig = GroupInfo%iresogrp1, GroupInfo%iresogrp2
!  WRITE(100, '(3I7)'), ig, icgb(ig), icge(ig)
!  DO ip = icgb(ig), icge(ig)
!    WRITE(95, '(I7, I4, ES14.6)') ip, ig, geg(ip)
!  END DO
!END DO
10 CONTINUE
IF (ig .LE. ng)   STOP 'PointXSRT_MOD.F90 ** SetMGnPWInfo ** The PWXS could not cover the MG range (Lower E)'
ipsst = icgb(GroupInfo%iresogrp1)
Esdmax = geg(ipsst)
invEsdmax = 1._8/Esdmax
ipsed = icge(GroupInfo%iresogrp2)
Esdmin = geg(ipsed)
!
DO ig = GroupInfo%iresogrp1, GroupInfo%iresogrp2
  delu_pw(icgb(ig), ig) = LOG(enbhel(ig)/geg(icgb(ig))) + 0.5 * del1pg
  delu_pw(icgb(ig)+1:icge(ig)-1, ig) = del1pg
  delu_pw(icge(ig), ig) = LOG(geg(icgE(ig))/enbhel(ig+1)) + 0.5 * del1pg
END DO
END SUBROUTINE

! SUBROUTINE to search index of Pij grid (ascending order of XS) for shadowing correction and Pij inter..
SUBROUTINE IdxWgtPGrid(XS, idx1, idx2, wgt1, wgt2)
IMPLICIT NONE
REAL :: XS
INTEGER :: IDX, idx1, idx2
REAL :: wgt1, wgt2
!IDX = INT((XS - PGrid_min)*invdelXS)
!idx1 = MAX(idx,1)
!idx2 = idx + 1
!wgt2 = (XS - Pij_Grid(idx1))*invdelXS
!wgt1 = 1._8-wgt2
IDX = INT(invdelXS*LOG(XS/PGrid_min))+1
idx1 = MAX(idx,1)
idx2 = idx1+1
wgt1 = LOG(Pij_Grid(idx2)/XS)*invdelXS
wgt2 = 1._8-wgt1
END SUBROUTINE

! SUBROUTINE to search index of point-wise energy grid (ascending order of E) 
SUBROUTINE IdxWgtPWEA(E, idx1, idx2, wgt1, wgt2)
REAL :: E
INTEGER :: IDX,idx1, idx2
REAL :: wgt1, wgt2

IDX = INT(invdel1pg*LOG(E/psm_emin))
idx1 = MAX(idx,1)
idx2 = idx1+1
wgt1 = LOG(geg(idx2)/E)*invdel1pg
wgt2 = 1._8-wgt1
END SUBROUTINE

! SUBROUTINE to search index of point-wise energy grid (descending order of E)
SUBROUTINE IdxWgtPWED(E, idx1, idx2, wgt1, wgt2)
REAL :: E
INTEGER :: IDX,idx1, idx2
REAL :: wgt1, wgt2
!IDX = INT(invdel1pg*LOG(psm_emax/E))+1
IDX = INT(invdel1pg*LOG(psm_emax/E))+1
idx1 = MAX(idx,1)
idx2 = idx1+1
wgt1 = LOG(E/geg(idx2))*invdel1pg
wgt2 = 1._8-wgt1
END SUBROUTINE

SUBROUTINE CalcShadowingCorrection(Core, Fxr, GroupInfo, THInfo, nTracerCntl, PE)
USE PARAM
USE FILES,          ONLY : IO8
USE IOUTIL,         ONLY : message
USE TYPEDEF,        ONLY : coreinfo_type,    Fxrinfo_type, ResVarPin_Type,  Cell_TYPE, Pin_TYPE,  GroupInfo_Type, PE_TYPE, THInfo_Type
USE CNTL,           ONLY : nTracerCntl_Type
USE TIMER,          ONLY : nTracer_dclock,       TimeChk 
USE OMP_LIB
USE TH_Mod,         ONLY : GetPinFuelTemp
#ifdef MPI_ENV
USE MPICOMM_MOD,    ONLY : MPI_SYNC,             REDUCE
#endif
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type),POINTER :: Fxr(:, :), myFxr
TYPE(THInfo_Type) :: THInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE
TYPE(ResVarPin_Type), POINTER :: ResVarPin(:,:),RP
TYPE(Pin_TYPE), POINTER :: PIN(:)
TYPE(Cell_TYPE), POINTER :: CellInfo(:)
INTEGER :: myzb, myze
logical :: master
real :: Tbeg, Tend, TempAvgsq
INTEGER :: iz, nxy, ipin, icel, nlocalfxr, fxridxst, ig
REAL :: D, A, alpha1, alpha2, beta1, beta2, Pe_lat, Pe_iso, XSe

myzb = PE%myzb; myze = PE%myze
master = PE%master

Tbeg = nTracer_Dclock(FALSE, FALSE)
nxy = Core%nxy
Pin => Core%Pin; CellInfo => Core%CellInfo
ResVarPin => Core%ResVarPin

WRITE(mesg,'(a,i3)') 'Calculating Pin-wise Shadowing Effect Correction Factor for PSM'
IF(master) CALL message(io8, TRUE, TRUE, mesg) 

DO iz = myzb, myze
  IF(.NOT. Core%lFuelPlane(iz)) CYCLE  
  DO ipin = 1, nxy
    if (.not.ResVarPin(ipin,iz)%lres) cycle ! Should AIC also be treated?
    !if (.not.ResVarPin(ipin,iz)%lfuel) cycle
    icel = Pin(ipin)%Cell(iz)
    RP => ResVarPin(ipin,iz)
    D = RP%Dancoff
!    D = 0.626261851999634
    XSE = RP%XSEsc
    ! Carlvik's two-term rational approx.
    A = D / (1-D)
    alpha1 = (5*A+6-SQRT(A**2+36*A+36))/(2*A+2)
    alpha2 = (5*A+6+SQRT(A**2+36*A+36))/(2*A+2)
    beta1 = ((4*A+6)/(A+1)-alpha1)/(alpha2-alpha1)
    beta2 = 1 - beta1
    DO ig = 1, nPGrid
      Pe_iso = 4*XSE/(Pij_Grid(ig)+2*XSE) - 3*XSE/(Pij_Grid(ig)+3*XSE)
      Pe_lat = beta1 * alpha1 * XSE / (Pij_Grid(ig)+alpha1*XSE) + beta2 * alpha2 * XSE / (Pij_Grid(ig)+alpha2*XSE)
      RP%eta(ig) = Pe_lat/Pe_iso
    END DO
    RP%icelPSM = CellInfo(iCel)%icelPSM
    RP%fuelvol = CellInfo(icel)%fuelvol
    RP%nonfuelvol = CellInfo(icel)%nonfuelvol
    RP%nfuelring = CellInfo(icel)%nfueldiv
!    RP%Pji => PSM_ISOPIN(CellInfo(iCel)%icelPSM)%Pji_iso ! Pointing Pin-wise Isolated First Flight Collision Prob. for PSM
  END DO ! ipin
END DO ! iz
  
!WRITE(mesg,'(a,i3)') 'Pointing Pin-wise Isolated First Flight Collision Prob. for PSM'
!IF(master) CALL message(io8, TRUE, TRUE, mesg) 
!DO iz = myzb, myze
!  IF(.NOT. Core%lFuelPlane(iz)) CYCLE  
!  DO ipin = 1, nxy
!    if (.not.ResVarPin(ipin,iz)%lres) cycle
!    icel = Pin(ipin)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
!    FxrIdxSt = Pin(ipin)%FxrIdxSt
!  END DO ! ipin
!END DO ! iz

END SUBROUTINE

SUBROUTINE CorePXSIsoUpdate(idiso, niso, iz)
!Update the lists of resonance isotopes
! PLANE 별로 RESONANT NUCLIDE 수를 세고 INDEX해주는 ROUTINE
USE XSLIB_MOD, ONLY : ldiso, mapnucl, npwxs
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
    DO i = 1, npwxs
      if (id .eq. IDPXS_COR(i,iz)) exit
    ENDDO
    if (i.gt.npwxs) then
      NPXS_COR(iz) = NPXS_COR(iz) + 1
      IDPXS_COR(NPXS_COR(iz),iz) = id
      mapnuclPXS(id,iz) = NPXS_COR(iz)
      mapnuclPXS(id+500,iz) = NPXS_COR(iz)
    endif
  enddo
END SUBROUTINE

SUBROUTINE UpdtCoreIsoInfoPSM(Core, Fxr, PE)
USE TYPEDEF,      ONLY : CoreInfo_Type,    FxrInfo_Type,   Pin_TYPE, ResVarPin_TYPE, Cell_TYPE, PE_Type   
USE XSLIB_MOD,    ONLY : nelthel, CoreResIsoUpdate
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type),POINTER :: Fxr(:, :), myFxr
TYPE(PE_TYPE) :: PE
TYPE(Pin_TYPE), POINTER :: PIN(:)
TYPE(ResVarPin_TYPE), POINTER :: ResVarPIN(:,:)
TYPE(Cell_TYPE), POINTER :: CellInfo(:)

INTEGER :: nxy, myzb, myze, xyb, xye, iz, ipin, FxrIdxSt, icel, nlocalFxr, j, ifxr 
INTEGER :: niso, nisotot, iso, nid, idx, idisotot(nelthel), jso
INTEGER, POINTER :: idiso(:)
REAL :: pnumtot(nelthel),vol,volsum
REAL, POINTER :: pnum(:)

Pin => Core%Pin; CellInfo => Core%CellInfo
ResVarPin => Core%ResVarPin
nxy = Core%nxy;  myzb = PE%myzb; myze = PE%myze

!--- CNJ Edit : Domain Decomposition + MPI
xyb = PE%myPinBeg; xye = PE%myPinEnd
IF (PE%RTMASTER) THEN
  xyb = 1; xye = nxy
ENDIF

DO iz = myzb, myze
  DO ipin = xyb, xye
      icel = Pin(ipin)%Cell(iz)
      nlocalFxr = CellInfo(icel)%nFxr
      FxrIdxSt = Pin(ipin)%FxrIdxSt
      
      idisotot = 0; pnumtot = 0._8; nisotot = 0; volsum = 0._8
      DO j =1, nLocalFxr
        ifxr = FxrIdxSt + j -1
        myFxr => Fxr(ifxr, iz)
        IF( .NOT. myFxr%lfuel ) CYCLE    
        !IF ( myFxr%lCLD ) CYCLE   
        niso = myFxr%niso
        pnum => myFxr%pnum; idiso => myFxr%idiso
        CALL CorePXSIsoUpdate(myFXR%idiso(1:niso), niso, iz)
        vol = myFxr%area
        DO iso = 1, niso
          nid = idiso(iso)
          DO idx = 1, nisotot
            if (nid.eq.idisotot(idx)) exit
          ENDDO
          if (idx.le.nisotot) then
            pnumtot(idx) =  pnumtot(idx) + pnum(iso) * vol
          else
            nisotot = nisotot + 1
            idisotot(nisotot) = nid
            pnumtot(nisotot) =  pnumtot(nisotot) + pnum(iso) * vol 
          endif
        ENDDO ! iso
        volsum = volsum + vol
      ENDDO ! j (fxr)
      DO iso = 1, nisotot
        pnumtot(iso) = pnumtot(iso) / volsum  
      ENDDO
      ResVarPin(ipin,iz)%niso = nisotot
      ResVarPin(ipin,iz)%idiso(1:nisotot)=idisotot(1:nisotot)
      ResVarPin(ipin,iz)%pnum(1:nisotot)=pnumtot(1:nisotot)           
      
  ENDDO
ENDDO
      
NULLIFY(Pin,ResVarPin,CellInfo,myFxr,idiso,pnum)
END SUBROUTINE

! Allocation and Deallocation of PSMMacXS_TYPE type
SUBROUTINE AllocPSMMacXs(XsMac)
IMPLICIT NONE
TYPE(PSMMacXs_Type) :: XsMac
INTEGER :: ng, nr
nr = XsMac%nfxr
ng = XsMac%ng
ALLOCATE(XsMac%xsmact(ng, nr))
ALLOCATE(XsMac%xsmactr(ng, nr))
ALLOCATE(XsMac%xsmaca(ng, nr))
ALLOCATE(XsMac%xsmacf(ng, nr))
ALLOCATE(XsMac%xsmacnf(ng, nr))
ALLOCATE(XsMac%xsmackf(ng, nr))
ALLOCATE(XsMac%xsmacs(ng, nr)) ! umm..
XsMac%lalloc = .TRUE.
END SUBROUTINE

SUBROUTINE FreePSMMacXs(XsMac)
IMPLICIT NONE
TYPE(PSMMacXs_Type) :: XsMac
INTEGER :: ng, nr
XsMac%nfxr = 0
XsMac%ng = 0
DEALLOCATE(XsMac%xsmact)
DEALLOCATE(XsMac%xsmactr)
DEALLOCATE(XsMac%xsmaca)
DEALLOCATE(XsMac%xsmacf)
DEALLOCATE(XsMac%xsmacnf)
DEALLOCATE(XsMac%xsmackf)
DEALLOCATE(XsMac%xsmacs) ! umm..
XsMac%lalloc = .FALSE.
END SUBROUTINE

SUBROUTINE AllocPSMMacXsIso(XsIsoMac, ng, niso, nr)
IMPLICIT NONE
INTEGER :: ng, niso, nr
TYPE(PSMMacXs_Type) :: XsIsoMac
IF (XsIsoMac%lIsoalloc) CALL FreePSMMacXsIso(XsIsoMac)
ALLOCATE(XsIsoMac%IsoMGXsMacA(niso,ng,nr))
ALLOCATE(XsIsoMac%IsoMGXsMacf(niso,ng,nr))
ALLOCATE(XsIsoMac%IsoMGXsMacnf(niso,ng,nr))
ALLOCATE(XsIsoMac%IsoMGMacnu(niso,ng,nr))
ALLOCATE(XsIsoMac%IsoMGXsMacS(niso,ng,nr))
ALLOCATE(XsIsoMac%IsoMGXsMacTr(niso,ng,nr))
ALLOCATE(XsIsoMac%IsoMGXsMacT(niso,ng,nr)) !BYS edit
ALLOCATE(XsIsoMac%IsoMGXsRadCap(niso,ng,nr))
XsIsoMac%lIsoalloc = .TRUE.
XsIsoMac%niso = niso
END SUBROUTINE

SUBROUTINE FreePSMMacXsIso(XsIsoMac)
IMPLICIT NONE
TYPE(PSMMacXs_Type) :: XsIsoMac
DEALLOCATE(XsIsoMac%IsoMGXsMacA)
DEALLOCATE(XsIsoMac%IsoMGXsMacf)
DEALLOCATE(XsIsoMac%IsoMGXsMacnf)
DEALLOCATE(XsIsoMac%IsoMGMacnu)
DEALLOCATE(XsIsoMac%IsoMGXsMacS)
DEALLOCATE(XsIsoMac%IsoMGXsMacTr)
DEALLOCATE(XsIsoMac%IsoMGXsMacT) !BYS edit
DEALLOCATE(XsIsoMac%IsoMGXsRadCap)
XsIsoMac%lIsoalloc = .FALSE.
XsIsoMac%niso = 0
END SUBROUTINE

SUBROUTINE AllocPSMSLDdata(XsMac)
IMPLICIT NONE
TYPE(PSMMacXs_Type) :: XsMac
INTEGER :: ng, nr
ng = XsMac%neg;
nr = XsMac%nrmax;
ALLOCATE(XsMac%phip(nr,0:ng))
ALLOCATE(XsMac%phip_(0:ng,nr))
ALLOCATE(XsMac%phim(0:ng))
!ALLOCATE(XsMac%SSrc(ng*3,0:nr))
ALLOCATE(XsMac%xst(nr,ng))
ALLOCATE(XsMac%xsa(ng,nr))
ALLOCATE(XsMac%xsf(ng,nr))
ALLOCATE(XsMac%xss(ng,3,nr))
ALLOCATE(XsMac%xss0(ng,nr))
ALLOCATE(XsMac%SSRCAL(ng,3,0:nr))
ALLOCATE(XsMac%favgm(3,nr))
ALLOCATE(XsMac%lfsrc(3,nr))
ALLOCATE(XsMac%alphaf(3,nr))
ALLOCATE(XsMac%invsubalphaf(3,nr))
ALLOCATE(XsMac%XSPF(3,nr))
ALLOCATE(XsMac%xspf_inval(3,nr))
ALLOCATE(XsMac%Qf(nr))
ALLOCATE(XsMac%Pim(nr),XsMac%Pei_iso(nr),XsMac%Pmi(nr),XsMac%Pji(nr,nr))
XsMac%lallocIsoDataSLD = .TRUE.
END SUBROUTINE

SUBROUTINE FreePSMSLDdata(XsMac)
IMPLICIT NONE
TYPE(PSMMacXs_Type) :: XsMac
INTEGER :: ng, nr
DEALLOCATE(XsMac%phip)
DEALLOCATE(XsMac%phip_)
DEALLOCATE(XsMac%phim)
DEALLOCATE(XsMac%xst)
DEALLOCATE(XsMac%xsa)
DEALLOCATE(XsMac%xsf)
DEALLOCATE(XsMac%xss)
DEALLOCATE(XsMac%SSRCAL)
DEALLOCATE(XsMac%favgm)
DEALLOCATE(XsMac%lfsrc)
DEALLOCATE(XsMac%alphaf)
DEALLOCATE(XsMac%invsubalphaf)
DEALLOCATE(XsMac%XSPF)
DEALLOCATE(XsMac%xspf_inval)
DEALLOCATE(XsMac%Qf)
DEALLOCATE(XsMac%Pim,XsMac%Pei_iso,XsMac%Pmi,XsMac%Pji)
XsMac%lallocIsoDataSLD = .FALSE.
END SUBROUTINE

!SUBROUTINE TempIntIdxWgtPSM(pxs, temp, it1, it2, wgt1, wgt2)
!USE XSLIB_MOD, ONLY : pwxs_TYPE
!IMPLICIT NONE
!INTEGER :: ntemp, it, it1, it2
!REAL :: wgt1, wgt2
!REAL :: temp
!type(pwxs_TYPE) :: pxs
!ntemp = pxs%ntemp
!DO it = 1, ntemp-1
!  IF ((pxs%temp(it)-temp)*(pxs%temp(it+1)-temp) .LE. 0._8) EXIT
!END DO
!IF (it .GE. ntemp) it = ntemp-1
!it1 = it
!it2 = it + 1
!wgt2 = (temp - pxs%temp(it1))/(pxs%temp(it2)-pxs%temp(it1))
!wgt1 = 1._8-wgt2
!END SUBROUTINE

SUBROUTINE TempIntIdxWgtPSM(pxs, temp, it1, it2, wgt1, wgt2)
USE XSLIB_MOD, ONLY : pwxs_TYPE
IMPLICIT NONE
INTEGER :: ntemp, it, it1, it2
REAL :: wgt1, wgt2
REAL :: temp, tempsq
type(pwxs_TYPE) :: pxs
ntemp = pxs%ntemp
tempsq = SQRT(temp)
DO it = 1, ntemp-1
  IF ((pxs%tempsq(it)-tempsq)*(pxs%tempsq(it+1)-tempsq) .LE. 0._8) EXIT
END DO
IF (it .GE. ntemp) it = ntemp-1
it1 = it
it2 = it + 1
wgt2 = (tempsq - pxs%tempsq(it1))/(pxs%tempsq(it2)-pxs%tempsq(it1))
wgt1 = 1._8-wgt2
END SUBROUTINE

! Calculate Effective XS of FXRs in the Pin with PSM
SUBROUTINE EffMacXsPSM(XsMac, Fxr, mypin, tempref, nring, iresogrp1, iresogrp2, ng, lIsoXsOut, iPSMcel)
USE TYPEDEF,   ONLY : FxrInfo_Type, Pin_Type, ResVarPin_Type
USE CNTL,      ONLY : nTracerCntl
USE XSLIB_MOD, ONLY : libdata, mapnucl, ldiso, pwxs_type, pwxs, nelthel
USE OMP_LIB
!USE TIMER,            ONLY : nTracer_dclock, TimeChk
#ifdef __PGI
USE IEEE_ARITHMETIC   !--- CNJ Edit : F2003 Standard
#endif
IMPLICIT NONE
TYPE(PSMMacXs_Type) :: XsMac     !Microscopic XS
TYPE(ResVarPin_Type) :: mypin
INTEGER :: nring
TYPE(FxrInfo_Type), TARGET :: Fxr(nring)
!
INTEGER,INTENT(IN) :: ng, iPSMcel, iresogrp1, iresogrp2
REAL,INTENT(IN) :: tempref
LOGICAL,INTENT(IN) :: lIsoXsOut
!
INTEGER :: iso, jso, id, jd, idres, jdres, repnid, repid, icat, idxsig(2,2), ir, kso, idxrat(2), kdres
INTEGER :: nlv, nlvflx, ilv, it, idxtemp(2), nThread
INTEGER :: ig, ifxr, i, j, niso, idpxs, it1, it2, idx1, idx2, img

REAL :: temp, TempAvgsq, Tempsq
REAL, POINTER :: Pji_iso(:,:,:)
REAL, POINTER :: phip(:,:), phip_(:,:), xst(:,:), Qf(:), phim(:), xss(:,:,:), xsa(:,:), xsf(:,:), SSRCAL(:,:,:), xss0(:,:)
REAL, POINTER, DIMENSION(:,:) :: xsmaca,xsmacf,xsmacnf,xsmacKf,xsmacs,xsmact
REAL, POINTER, DIMENSION(:,:,:) :: IsoMGXSMacA, IsoMGXsMacS, IsoMGXsMacT, IsoMGXsMacF
REAL :: xspm, fringvol, fuelvol, nonfuelvol, xstf, Qm, invnonfuelsrc
TYPE(FxrInfo_Type), POINTER :: myFXR
!
REAL,POINTER,DIMENSION(:)   :: pnum
INTEGER,POINTER :: idiso(:)
TYPE(libdata), POINTER :: isodata, jsodata, repisodata
TYPE(pwxs_type), pointer :: pxs
REAL :: Pmm, eta, Pjicorr,invxspm,ndsum(3)
REAL, POINTER, DIMENSION(:) :: Pim, Pei_iso, Pmi
REAL, POINTER, DIMENSION(:,:) :: Pji
REAL, POINTER, DIMENSION(:,:) :: favgm, alphaf, invsubalphaf, xspf, xspf_inval
REAL ::  siglpMcat(3), sigpm_inval(3), mavgm(3), alpham(3), invsubalpham(3), ND, wgt1, wgt2, xssp(3), mgflx,delphi, mgregflx
LOGICAL,POINTER :: lfsrc(:,:)
LOGICAL :: lmsrc(3)
!REAL :: TimeBeg(15), TimeEnd(15), ElaspedTime(15)
!
IF (.NOT.myPin%lfuel) RETURN
!
!ElaspedTime = 0._8
!TimeBeg(1) = nTracer_dclock(FALSE,FALSE)

!IF(XsMac%lallocIsoDataSLD) CALL FreePSMSLDdata(XsMac)
IF (.NOT.XsMac%lallocIsoDataSLD) THEN
  XsMac%neg = neg
  XsMac%nrmax = nrmax
  CALL AllocPSMSLDdata(XsMac)
END IF
!
!IF(XsMac%lalloc) CALL FreePSMMacXs(XsMac)
IF(.NOT. XsMac%lalloc) THEN
  XsMac%ng = ng
  XsMac%nfxr = nrmax
  CALL AllocPSMMacXs(XsMac)
ENDIF  
!
!IF(XsMac%lIsoAlloc) CALL FreePSMMacXsIso(XsMac)
IF( lIsoXsOut .AND. ( .NOT. XsMac%lIsoAlloc) ) THEN
  CALL AllocPSMMacXsIso(XsMac, ng, nelthel, nrmax) ! 19/09/24 using niso considering increasable niso |--Edit by JSU
ENDIF
!
Pji_iso => PSM_ISOPIN(iPSMcel)%Pji_iso
!
xsmact      => XsMac%xsmact
xsmaca      => XsMac%xsmaca
xsmacf      => XsMac%xsmacf
xsmacnf     => XsMac%xsMacnF
xsmackf     => XsMac%xsmackf
xsmacs      => XsMac%xsmacs
!
Pim         => XsMac%Pim    
Pei_iso     => XsMac%Pei_iso
Pmi         => XsMac%Pmi
Pji         => XsMac%Pji
favgm       => XsMac%favgm        
alphaf      => XsMac%alphaf       
invsubalphaf=> XsMac%invsubalphaf 
lfsrc       => XsMac%lfSRC
xspf        => XsMac%xspf
xspf_inval  => XsMac%xspf_inval
!
IsoMGXsMacA => XsMac%IsoMGXsMacA
IsoMGXsMacS => XsMac%IsoMGXsMacS
IsoMGXsMacF => XsMac%IsoMGXsMacF
phip        => XsMac%phip
phip_       => XsMac%phip_
phim        => XsMac%phim
xst         => XsMac%xst
xsa         => XsMac%xsa
xsf         => XsMac%xsf
xss         => XsMac%xss
xss0        => XsMac%xss0
SSRCAL      => XsMac%SSRCAL
Qf          => XsMac%Qf
!
fuelvol      = mypin%fuelvol
nonfuelvol   = mypin%nonfuelvol
siglpMcat    = myPin%siglpMcat   
mavgm        = myPin%mavgm       
alpham       = myPin%alpham      
invsubalpham = myPin%invsubalpham
sigpm_inval  = mypin%sigpm_inval
lmsrc        = mypin%lmsrc
!
fringvol = fuelvol/nring
xspm     = mypin%siglpM
DO i = 1, 3
  xspm = xspm - sigpm_inval(i) * 0.5_8 * del1pg
END DO
invxspm = 1._8/ xspm
invnonfuelsrc = 1._8/(xspm * nonfuelvol)
!invnonfuelsrc = 1._8/(nonfuelvol)
!
TempAvgsq = tempref 
!#define PSM_debug_noABS
!
! Set point-wise XS of each FXR (xst in fuels -> temperature difference/ND consideration)
xst = 0._8; xspf = 0._8
xss = 0._8; xsf = 0._8; xsa = 0._8
favgm = 0._8
lfsrc = .TRUE.
DO ir = 1, nring
  ifxr = nring - ir + 1
  myFxr => FXR(ifxr)
  niso = myFXR%niso
  pnum  => myFXR%pnum
  idiso => myFXR%idiso
  ndsum = 0._8
  DO iso = 1, niso
    id = mapnucl(idiso(iso));    isodata => ldiso(id)
    idpxs = isodata%pxsid
    ND = pnum(iso)
    favgm(isodata%icatf,ir) = favgm(isodata%icatf,ir) + pnum(iso) * isodata%aw
    ndsum(isodata%icatf) = ndsum(isodata%icatf) + pnum(iso)
    xspf(isodata%icatf,ir) = xspf(isodata%icatf,ir) + ND*isodata%sigp
    IF (idpxs .EQ. 0) THEN
      DO ig = 1, ipsed
        !xst(ir,ig) = xst(ir,ig) + ND*isodata%sigp
        xss(ig,isodata%icatf,ir) = xss(ig,isodata%icatf,ir) + ND*isodata%sigp
      END DO
    ELSE
      pxs => pwxs(idpxs)
      CALL TempIntIdxWgtPSM(pxs, myFXR%temp, it1, it2, wgt1, wgt2)
      DO ig = 1, ipsed
        !xst(ir,ig) = xst(ir,ig) + ND*(pxs%xst(ig,it1)*wgt1+pxs%xst(ig,it2)*wgt2)
        xsa(ig,ir) = xsa(ig,ir) + ND*(pxs%xsa(ig,it1)*wgt1+pxs%xsa(ig,it2)*wgt2)
        xss(ig,isodata%icatf,ir) = xss(ig,isodata%icatf,ir) + ND*(pxs%xss(ig,it1)*wgt1+pxs%xss(ig,it2)*wgt2)
        IF (pxs%lfis) xsf(ig,ir) = xsf(ig,ir) + ND*(pxs%xsf(ig,it1)*wgt1+pxs%xsf(ig,it2)*wgt2)
      END DO
    END IF
  END DO ! iso
  DO ig = 1, ipsed
    xst(ir,ig) = xsa(ig,ir)
    xss0(ig,ir) = 0._8
    DO i = 1, 3
      xss0(ig,ir) = xss0(ig,ir) + xss(ig,i,ir)
      xst(ir,ig) = xst(ir,ig) + xss(ig,i,ir)
    END DO
  END DO
  DO i = 1, 3
    IF (ndsum(i).LT.1.E-20_8) THEN
      lfsrc(i,ir) = .FALSE.
      CYCLE
    END IF
    favgm(i,ir) = favgm(i,ir)/ndsum(i)
    alphaf(i,ir) = ((favgm(i,ir)-1._8)/(favgm(i,ir)+1._8))**2._8
    invsubalphaf(i,ir) = 1._8/(1._8-alphaf(i,ir))
    xspf_inval(i,ir) = xspf(i,ir) * invsubalphaf(i,ir)
    xss(:,i,ir) = xss(:,i,ir) * invsubalphaf(i,ir) * 0.5_8 * del1pg
  END DO
  DO ig = 1, ipsed
    DO i=1,3
      xst(ir,ig) = xst(ir,ig) - xss(ig,i,ir)!*invsubalphaf(i,ir)*0.5*del1pg
    END DO
  END DO
END DO ! ir
!TimeEnd(1) = nTracer_dclock(FALSE,FALSE)
! Solve Slowingdown calculation
!TimeBeg(2) = nTracer_dclock(FALSE,FALSE)
phip = 0._8
phip_ = 0._8
phim = 0._8
SSrcal = 0._8
DO i = 1, nring
  phip(i,0:ipsst-1)  = 1._8/geg(0:ipsst-1)
  phip_(0:ipsst-1,i) = 1._8/geg(0:ipsst-1)
END DO
phim(0:ipsst-1) = 1._8/geg(0:ipsst-1)
! source at (ipsst-1)-th point
Qf = 0._8
DO ir = 1, nring
  DO i = 1, 3
    IF (.NOT.lfsrc(i,ir)) CYCLE
    Qf(ir) = Qf(ir) + xspf(i,ir)*(1._8-del1pg*0.5_8*invsubalphaf(i,ir))!/geg(ipsst-1)
  END DO
END DO
Qm = 0._8
DO i = 1, 3
  Qm = Qm + siglpMcat(i)*(1._8-del1pg*0.5_8*invsubalpham(i))!/geg(ipsst-1)
END DO
xssp = sigpm_inval * 0.5_8 * del1pg
!
DO ir = 1, nring
  DO i = 1, 3
    xspf(i,ir) = xspf(i,ir)*(1-invsubalphaf(i,ir) * 0.5_8 * del1pg)
  END DO
END DO

IG = ipsst-1
!xstf = SUM(xst(:,ig)*phip(:,ig-1))/SUM(phip(:,ig-1))
xstf = SUM(xspf(:,1))
!WRITE(86, '(I7, 20ES14.6)') ig, geg(ig), xstf, (Phip(i,ig), i = 1, nring), Phim(ig)
!WRITE(87, '(I7, 20ES14.6)') ig, geg(IG), xss0(1,1), (Qf(i), i = 1, nring), Qm
!TimeEnd(2) = nTracer_dclock(FALSE,FALSE)
!TimeBeg(3) = nTracer_dclock(FALSE,FALSE)
!TimeBeg(13) = nTracer_dclock(FALSE,FALSE)
#ifdef PSM_debug_noABS
DO ir = 1, nring
  DO i =1, 3
    SSRCAL(ipsst-1,i,ir) = phip_(ipsst-1,ir)*xspf_inval(i,ir)*0.5_8*del1pg
  END DO
END DO
#else
DO ir = 1, nring
  DO i =1, 3
    SSRCAL(ipsst-1,i,ir) = phip_(ipsst-1,ir)*xss(ipsst-1,i,ir)
  END DO
END DO
#endif
DO i =1, 3
  SSRCAL(ipsst-1,i,0) = phim(ipsst-1)*xssp(i)
END DO
!TimeEnd(13) = nTracer_dclock(FALSE,FALSE)
!ElaspedTime(13) = TimeEnd(13) - TimeBeg(13)
DO ig = ipsst, ipsed
  !TimeBeg(5) = nTracer_dclock(FALSE,FALSE)
#ifdef PSM_debug_noABS
  xstf = SUM(xspf(:,1))
#else
  xstf = SUM(xst(:,ig)*phip(:,ig-1))/SUM(phip(:,ig-1))
#endif
  ! Src update (Qm and Q)
  CALL UpdtPSMSrc(geg(ig),geg(ig-1),IG,Qf,Qm,SSRcal,xspf_inval,alphaf,            &
                  sigpm_inval,alpham,nring,lfsrc,lmsrc)!,ElaspedTime(6:10))
  !TimeEnd(5) = nTracer_dclock(FALSE,FALSE)
  !ElaspedTime(5) = ElaspedTime(5) + TimeEnd(5) - TimeBeg(5)
!  WRITE(87, '(I7, 30ES14.6)') ig, geg(ig), xss0(ig,1), (Qf(i), i = 1, nring), Qm
  ! Set Collision Probabilities ******************************************************
  ! Find Shadowing correction factor according to xstf
  !TimeBeg(11) = nTracer_dclock(FALSE,FALSE)
  CALL IdxWgtPGrid(xstf, idx1, idx2, wgt1, wgt2)
  eta = myPin%eta(idx1)*wgt1 + myPin%eta(idx2)*wgt2
  Pmm = 1._8
  DO i = 1, nring
    Pei_iso(i) = 1._8
    DO j = 1, nring
      ! Find Pij_iso and Pei_iso according to xstf, note that this Pji is FFCP i -> j
      Pji(j,i) = Pji_iso(j,i,idx1)*wgt1 + Pji_iso(j,i,idx2)*wgt2
      Pei_iso(i) = Pei_iso(i) - Pji(j,i)
    END DO
    ! Find Pim (or Pei) and correction factor of Pij
    Pim(i) = Pei_iso(i)*eta
    Pjicorr = (1._8-Pim(i))/(1._8-Pei_iso(i))
    ! Find Pij
#ifdef PSM_debug_noABS
    DO j = 1, nring
      Pji(j,i) = Pji(j,i)*Pjicorr/xstf
    END DO
#else
    DO j = 1, nring
      Pji(j,i) = Pji(j,i)*Pjicorr/xst(j,ig)
    END DO
#endif
    ! Find PM,i and Pmm
    Pmi(i) = xstf*fringvol*Pim(i)*invnonfuelsrc
    Pmm = Pmm - Pmi(i)
#ifdef PSM_debug_noABS
    Pmi(i) = Pmi(i)/xstf
#else
    Pmi(i) = Pmi(i)/xst(i,ig)
#endif
  END DO ! i (ring idx)
!  TimeEnd(11) = nTracer_dclock(FALSE,FALSE)
!  ElaspedTime(11) = ElaspedTime(11) + TimeEnd(11) - TimeBeg(11)
  ! Calc flux of each region ********************************************************
!  TimeBeg(12) = nTracer_dclock(FALSE,FALSE)
  phim(ig) = Pmm*Qm*invxspm
  DO i = 1, nring
    phip(i,ig) = PiM(i)*Qm*invxspm
    DO j = 1, nring
      phip(i,ig) = phip(i,ig) + Pji(j,i)*Qf(j)
    END DO
    phip_(ig,i) = phip(i,ig)
    phim(ig) = phim(ig) + Pmi(i)*Qf(i)
  END DO ! i (ring idx)
  WRITE(86, '(I7, 20ES14.6)') ig, geg(ig), xstf, (Phip(i,ig), i = 1, nring), Phim(ig)
  phip(:,ig) = phip(:,ig)/geg(ig)
  phim(ig) = phim(ig)/geg(ig)
  DO ir = 1, nring
    phip_(ig,ir) = phip(ir,ig)
  END DO
  !WRITE(86, '(I7, 20ES14.6)') ig, geg(ig), xstf, (Phip(i,ig)*geg(ig), i = 1, nring), Phim(ig)*geg(ig)
#ifdef PSM_debug_noABS
  phim(ig) = 1._8/geg(ig)
  phip(:,ig) = 1._8/geg(ig)
  DO i = 1, nring
    phip_(ig,i) = 1._8/geg(ig)
  END DO
#endif
#ifdef PSM_debug_noABS
  DO ir = 1, nring
    DO i =1, 3
      SSRcal(ig,i,ir) = phip_(ig,ir)*xspf_inval(i,ir)*0.5_8*del1pg
    END DO
  END DO
#else
  DO ir = 1, nring
    DO i =1, 3
      SSRcal(ig,i,ir) = phip_(ig,ir)*xss(ig,i,ir)
    END DO
  END DO
#endif
  DO i =1, 3
    SSRcal(ig,i,0) = phim(ig)*xssp(i)
  END DO
!  TimeEnd(12) = nTracer_dclock(FALSE,FALSE)
!  ElaspedTime(12) = ElaspedTime(12) + TimeEnd(12) - TimeBeg(12)
END DO ! IG
!TimeEnd(3) = nTracer_dclock(FALSE,FALSE)
!TimeBeg(4) = nTracer_dclock(FALSE,FALSE)
!
xsmaca = 0._8
xsmacf = 0._8
xsmacs = 0._8
IsoMGXSMacA = 0._8
IsoMGXSMacS = 0._8
IsoMGXSMacF = 0._8
! Condensing the MGXS and flux     
DO ir = 1, nring
  ifxr = nring - ir + 1
  myFxr => FXR(ifxr);  niso = myFXR%niso
  pnum  => myFXR%pnum;  idiso => myFXR%idiso
  DO img = iresogrp1, iresogrp2
    mgregflx = 0._8
    DO ig = icgb(img), icge(img)
      delphi = delu_pw(ig, img)*phip(ir,ig)*geg(ig)
      mgregflx = mgregflx + delphi
      xsmaca(img, ir) = xsmaca(img, ir) + xsa(ig,ir)*delphi
      xsmacs(img, ir) = xsmacs(img, ir) + xss0(ig,ir)*delphi
      xsmacf(img, ir) = xsmacf(img, ir) + xsf(ig,ir)*delphi
    END DO
    xsmaca(img, ir) = xsmaca(img, ir) / mgregflx
    xsmacs(img, ir) = xsmacs(img, ir) / mgregflx
    xsmacf(img, ir) = xsmacf(img, ir) / mgregflx    
  END DO
  DO iso = 1, niso
    id = mapnucl(idiso(iso));    isodata => ldiso(id)
    idpxs = isodata%pxsid;     ND = pnum(iso)
    IF (idpxs .EQ. 0) CYCLE
    pxs => pwxs(idpxs)
    CALL TempIntIdxWgtPSM(pxs, myFXR%temp, it1, it2, wgt1, wgt2)
    IF (pxs%lfis) THEN
      DO img = iresogrp1, iresogrp2
        mgflx = 0._8
        DO ig = icgb(img), icge(img)
          delphi = delu_pw(ig, img)*phip(ir,ig)*geg(ig)
          mgflx = mgflx + delphi
          IsoMGXSMacA(iso,img,ir) = IsoMGXSMacA(iso,img,ir) + (pxs%xsa(ig,it1)*wgt1+pxs%xsa(ig,it2)*wgt2)*delphi*ND
          IsoMGXSMacS(iso,img,ir) = IsoMGXSMacS(iso,img,ir) + (pxs%xss(ig,it1)*wgt1+pxs%xss(ig,it2)*wgt2)*delphi*ND
          IF (pxs%lfis) IsoMGXSMacF(iso,img,ir) = IsoMGXSMacF(iso,img,ir) + (pxs%xsf(ig,it1)*wgt1+pxs%xsf(ig,it2)*wgt2)*delphi*ND
        END DO ! ig ENERGY POINT
        IsoMGXSMacA(iso,img,ir) = IsoMGXSMacA(iso,img,ir) / mgflx
        IsoMGXSMacS(iso,img,ir) = IsoMGXSMacS(iso,img,ir) / mgflx
        IsoMGXSMacF(iso,img,ir) = IsoMGXSMacF(iso,img,ir) / mgflx
      END DO ! img MULTIGROUP
    ELSE
      DO img = iresogrp1, iresogrp2
        mgflx = 0._8
        DO ig = icgb(img), icge(img)
          delphi = delu_pw(ig, img)*phip(i,ig)*geg(ig)
          mgflx = mgflx + delphi
          IsoMGXSMacA(iso,img,ir) = IsoMGXSMacA(iso,img,ir) + (pxs%xsa(ig,it1)*wgt1+pxs%xsa(ig,it2)*wgt2)*delphi*ND
          IsoMGXSMacS(iso,img,ir) = IsoMGXSMacS(iso,img,ir) + (pxs%xss(ig,it1)*wgt1+pxs%xss(ig,it2)*wgt2)*delphi*ND
        END DO ! ig ENERGY POINT
        IsoMGXSMacA(iso,img,ir) = IsoMGXSMacA(iso,img,ir) / mgflx
        IsoMGXSMacS(iso,img,ir) = IsoMGXSMacS(iso,img,ir) / mgflx
      END DO ! img MULTIGROUP
    END IF
!    WRITE(88, '(I3, I7, 50ES14.6)') ir, idiso(iso), (IsoMGXSMacA(iso,img,ir), img=iresogrp1, iresogrp2), (IsoMGXSMacS(iso,img,ir), img=iresogrp1, iresogrp2), (IsoMGXSMacF(iso,img,ir), img=iresogrp1, iresogrp2)
  END DO ! iso NUCLIDE
END DO ! ir RING
xsmacnf = xsmacf
xsmackf = xsmacf
!TimeEnd(4) = nTracer_dclock(FALSE,FALSE)
!DO i = 1, 4
!  ElaspedTime(i) = TimeEnd(i) - TimeBeg(i)
!END DO
!TimeChk%PSMTime = TimeChk%PSMTime + (TimeEnd(4) - TimeBeg(1))
!PRINT '(10ES14.6)', (ElaspedTime(i),i=1,7),ElaspedTime(11),ElaspedTime(12),ElaspedTime(13)

END SUBROUTINE

SUBROUTINE UpdtPSMSrc(E,Eold,ig,Qf,Qm,SSRcal,sigpf,alphaf,sigpm, &
                      alpham,nring,lfsrc,lmsrc)!,timemeasure)
!USE TIMER,            ONLY : nTracer_dclock
IMPLICIT NONE
REAL,POINTER :: sigpf(:,:),alphaf(:,:),SSRcal(:,:,:)
LOGICAL :: lfsrc(3,nring), lmsrc(3)
REAL :: sigpm(3),alpham(3),Qf(nring), Qm!, sigpm_inval(3)! sigpm is sigp/(1-a) form
REAL :: E, Eold, fEold, fE
REAL :: alpha(3), S, Stmp!, timemeasure(5)
INTEGER :: nring
INTEGER :: idxold1, idxold2, idxnew1, idxnew2, idx1, idx2, ir, i, ig, igfrom,j
REAL    :: wgtold1, wgtold2, wgtnew1, wgtnew2, srcold, srcnew, delutemp, addsrc, subsrc, time1, time2
!ALLOCATE(sSrc(ipsst-1:ig))
Qf = Qf*E/Eold
Qm = Qm*E/Eold
!time1 = nTracer_dclock(FALSE,FALSE)
DO ir = 1, nring
  addsrc = 0._8
  subsrc = 0._8
  DO i = 1,3
    IF (.NOT.lfsrc(i,ir)) CYCLE
    fE = E/alphaf(i,IR)
    fEold = Eold/alphaf(i,ir)
    idx1 = ipsst;idx2 = ipsst-1
    addsrc = addsrc + SSrcal(ig-1,i,ir)*2
    ! Source between E/alpha ~ Eold/alpha (need to be subtracted)
    IF (fE.LE.Esdmax) THEN
      CALL IdxWgtPWED(fE, idxnew1, idxnew2, wgtnew1, wgtnew2)
      idx2 = idxnew1
      srcnew = SSrcal(idxnew1,i,ir)*wgtnew1+SSrcal(idxnew2,i,ir)*wgtnew2
      delutemp = LOG(geg(idx2)/fE)/del1pg
      subsrc = subsrc + (srcnew + SSrcal(idx2,i,ir))*delutemp ! lethargy
      IF (fEold.LE.Esdmax) THEN
        CALL IdxWgtPWED(fEold, idxold1, idxold2, wgtold1, wgtold2)
        idx1 = idxold2
        srcold = SSrcal(idxold1,i,ir)*wgtold1+SSrcal(idxold2,i,ir)*wgtold2
        delutemp = LOG(fEold/geg(idx1))/del1pg
        subsrc = subsrc + (srcold + SSrcal(idx1,i,ir))*delutemp ! lethargy
      ELSE
        subsrc = subsrc + sigpf(i,ir)*(1._8/Esdmax - 1._8/fEold) ! lethargy
      END IF
      !
      DO igfrom= idx1, idx2-1
        subsrc = subsrc + (SSrcal(igfrom,i,ir)+SSrcal(igfrom+1,i,ir)) ! lethargy
      END DO
    ELSE
      subsrc = subsrc + alphaf(i,ir)*sigpf(i,ir)*(1._8/E - 1._8/Eold) ! lethargy
    END IF
  END DO ! i (categorized scatterer)
  Qf(ir) = Qf(ir) + (addsrc - subsrc)*E
END DO ! nring
!time2 = nTracer_dclock(FALSE,FALSE)
!timemeasure(1) = timemeasure(1) + time2-time1
!time1 = nTracer_dclock(FALSE,FALSE)
! moderator source
addsrc = 0._8
subsrc = 0._8
DO i = 1,3
  IF (.NOT.lmsrc(i)) CYCLE
  fE = E/alpham(i)
  fEold = Eold/alpham(i)
  idx1 = ipsst;idx2 = ipsst-1
  addsrc = addsrc + SSrcal(ig-1,i,0)*2
  ! Source between E/alpha ~ Eold/alpha (need to be subtracted)
  IF (fE.LE.Esdmax) THEN
    CALL IdxWgtPWED(fE, idxnew1, idxnew2, wgtnew1, wgtnew2)
    idx2 = idxnew1
    srcnew = SSrcal(idxnew1,i,0)*wgtnew1+SSrcal(idxnew2,i,0)*wgtnew2
    delutemp = LOG(geg(idx2)/fE)/del1pg
    subsrc = subsrc + (srcnew + SSrcal(idx2,i,0))*delutemp ! lethargy
    !
    IF (fEold.LE.Esdmax) THEN
      CALL IdxWgtPWED(fEold, idxold1, idxold2, wgtold1, wgtold2)
      idx1 = idxold2
      srcold = SSrcal(idxold1,i,0)*wgtold1+SSrcal(idxold2,i,0)*wgtold2
      delutemp = LOG(fEold/geg(idx1))/del1pg
      subsrc = subsrc + (srcold + SSrcal(idx1,i,0))*delutemp ! lethargy
    ELSE
      subsrc = subsrc + sigpm(i)*(1._8/Esdmax - 1._8/fEold) ! lethargy
    END IF
    !
    DO igfrom= idx1, idx2-1
      subsrc = subsrc + (SSrcal(igfrom,i,0)+SSrcal(igfrom+1,i,0)) ! lethargy
    END DO
  ELSE
    subsrc = subsrc + alpham(i)*sigpm(i)*(1._8/E - 1._8/Eold) ! lethargy
  END IF
END DO ! i (categorized scatterer)
Qm = Qm + (addsrc - subsrc)*E
!time2 = nTracer_dclock(FALSE,FALSE)
!timemeasure(2) = timemeasure(2) + time2-time1

END SUBROUTINE
!
!! Subroutine to generate effective XS with PSM
!SUBROUTINE PSMEffXSGen(Core, Fxr, THInfo, eigv, GroupInfo, nTracerCntl, PE)
!!Effective Xs Generation
!USE PARAM
!USE TYPEDEF,        ONLY : CoreInfo_Type,       FxrInfo_Type,          GroupInfo_Type,     &
!                           PE_Type,             THInfo_Type,           pin_type,           &
!                           ResVarpin_type,      cell_type,             XsMac_Type
!USE CNTL,           ONLY : nTracerCntl_Type
!USE MacXsLib_mod,   ONLY : MacXsBase
!USE BasicOperation, ONLY : CP_VA
!USE IOUTIL,         ONLY : message
!USE FILES,          ONLY : io8
!USE TH_Mod,         ONLY : GetPinFuelTemp
!USE XSLIB_MOD,      ONLY : libdata, ldiso, mapnucl
!USE MPIComm_Mod,    ONLY : MPI_SYNC
!USE OMP_LIB
!USE TIMER,            ONLY : nTracer_dclock, TimeChk
!!
!IMPLICIT NONE
!!
!TYPE(CoreInfo_Type) :: Core
!TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
!TYPE(THInfo_Type) :: THInfo
!REAL :: eigv
!TYPE(GroupInfo_TYPE) :: GroupInfo
!TYPE(nTracerCntl_Type) :: nTracerCntl
!TYPE(PE_Type) :: PE
!!
!TYPE(Pin_Type), POINTER :: PIN(:)
!TYPE(ResVarPin_Type), POINTER :: ResVarPin(:,:),mypin
!TYPE(Cell_Type), POINTER :: CellInfo(:), myCell
!TYPE(FxrInfo_Type), POINTER :: myFxr
!TYPE(XsMac_Type), SAVE :: XsMac(nThreadMax)
!TYPE(PSMMacXs_Type), SAVE :: PSMXsMac(nThreadMax)
!!
!INTEGER :: nxy, myzb, myze, iz, ipin, xyb, xye, FxrIdxSt, nlocalFxr, j, ifxr, icel  !--- (xyb, xye)CNJ Edit : Domain Decomposition + MPI
!INTEGER :: ng, ig, iresoGrp1, iresoGrp2, niso, iso, tid
!REAL :: PinFuelTempAvgsq, XsMacAold, XsMacFold, XsMacNFold, XsMackFold, XsMacSold, XsMacStrold,ndsumcat(3),volsum
!REAL :: isoXsMacAold(500),isoXsMacFold(500),isoXsMacNFold(500),isoXsMackFold(500),isoXsMacSold(500),isoXsMacS1old(500),isoXsMacSSold(500),isoxsmaccapold(500)
!LOGICAL :: lxslib, lAIC, lmcxs
!INTEGER :: nfueldiv, idxfuelsrt, iPSMcel, id, i, ifuelring, idpxs
!TYPE(LIBDATA), POINTER :: isodata
!REAL :: TimeEnd, TimeBeg 
!INTEGER, POINTER :: idiso(:)
!!
!IF (.not.nTracerCntl%lrestrmt) RETURN
!lxslib  = nTracerCntl%lxslib
!IF(.NOT. lxslib) RETURN
!!
!WRITE(mesg, '(A)') 'Update PSM Effective XSs...'
!IF(PE%Master) CALL message(io8, TRUE, TRUE, mesg)
!TimeBeg = nTracer_dclock(FALSE,FALSE)
!Pin => Core%Pin; CellInfo => Core%CellInfo
!ResVarPin => Core%ResVarPin;
!nxy = Core%nxy;  myzb = PE%myzb; myze = PE%myze
!
!!--- CNJ Edit : Domain Decomposition + MPI
!xyb = PE%myPinBeg; xye = PE%myPinEnd
!IF (PE%RTMASTER) THEN
!  xyb = 1; xye = nxy
!ENDIF
!
!ng = GroupInfo%ng
!iresoGrp1 = GroupInfo%nofg + 1; iresoGrp2 = GroupInfo%nofg + GroupInfo%norg
!
!CALL UpdtCoreIsoInfoPSM(Core, Fxr, PE)
!!
!CALL omp_set_num_threads(PE%nThread)
!tid = 1
!DO iz = myzb, myze
!  IF (.NOT. Core%lFuelPlane(iz)) CYCLE
!  !$OMP PARALLEL DEFAULT(SHARED)      &
!  !$OMP PRIVATE(tid, ipin, mypin, FxrIdxSt, icel, myCell, iPSMcel, nLocalFxr, PinFuelTempAvgsq, ndsumcat, volsum, j, ifxr, myfxr, niso, &
!  !$OMP         iso, id, isodata, nfueldiv, idxfuelsrt, ig, idiso, ifuelring, idpxs, &
!  !$OMP         isoXsMacAold,isoXsMacFold,isoXsMacNFold,isoXsMackFold,isoXsMacSold,isoXsMacS1old,isoXsMacSSold,isoxsmaccapold)
!  !$ tid = omp_get_thread_num()+1
!  !$OMP DO
!  DO ipin = xyb, xye
!    mypin => ResVarPin(ipin,iz)
!    FxrIdxSt = Pin(ipin)%FxrIdxSt
!    icel = Pin(ipin)%Cell(iz)
!    myCell => CellInfo(icel)
!    IF (.NOT.myCell%lPSMcel) CYCLE
!    iPSMcel = myCell%icelPSM
!    nlocalFxr = CellInfo(icel)%nFxr
!    PinFuelTempAvgsq = dsqrt(GetPinFuelTemp(Core, Fxr, iz, ipin))
!    ! Set the Potential, total XS o each FXR  in case the moderator condition changes -> after debug sink to TH module
!    mypin%siglpM = 0._8
!    mypin%siglpMcat = 0._8
!    mypin%mavgm = 0._8
!    ndsumcat = 0._8
!    volsum = 0._8
!    DO j = 1, nLocalFxr
!      ifxr = FxrIdxSt + j -1
!      myFxr => Fxr(ifxr, iz)
!      IF( myFxr%lfuel ) EXIT
!      niso = myFxr%niso
!      volsum = volsum + myCell%fxrvol(nLocalFxr-j+1)
!      DO iso = 1, niso
!        id = mapnucl(myFxr%idiso(iso));
!        isodata => ldiso(id)
!        ndsumcat(isodata%icatm) = ndsumcat(isodata%icatm) + myFXR%pnum(iso) * myCell%fxrvol(nLocalFxr-j+1)
!        mypin%mavgm(isodata%icatm) = mypin%mavgm(isodata%icatm) + myFXR%pnum(iso) * isodata%aw * myCell%fxrvol(nLocalFxr-j+1)
!        IF (isodata%sigp.EQ. 0._8) THEN
!          mypin%siglpM = mypin%siglpM + myFXR%pnum(iso) * isodata%lamsigp1G * myCell%fxrvol(nLocalFxr-j+1)
!          mypin%siglpMcat(isodata%icatm) = mypin%siglpMcat(isodata%icatm) + myFXR%pnum(iso) * isodata%lamsigp1G * myCell%fxrvol(nLocalFxr-j+1)
!        ELSE
!          mypin%siglpM = mypin%siglpM + myFXR%pnum(iso) * isodata%sigp * myCell%fxrvol(nLocalFxr-j+1)
!          mypin%siglpMcat(isodata%icatm) = mypin%siglpMcat(isodata%icatm) + myFXR%pnum(iso) * isodata%sigp * myCell%fxrvol(nLocalFxr-j+1)
!        END IF
!      END DO
!    END DO
!    mypin%siglpM = mypin%siglpM / volsum
!    myPin%lmsrc = .TRUE.
!    DO i = 1, 3
!      IF (ndsumcat(i) .LT. 1E-20) THEN
!        myPin%lmsrc(i) = .FALSE.
!        CYCLE
!      END IF
!      mypin%siglpMcat(i) = mypin%siglpMcat(i) / volsum
!      mypin%mavgm(i) = mypin%mavgm(i) / ndsumcat(i)
!      mypin%alpham(i) = ((mypin%mavgm(i)-1._8)/(mypin%mavgm(i)+1._8))**2
!      mypin%invsubalpham(i) = 1._8/(1._8-mypin%alpham(i))
!      mypin%sigpm_inval(i) = mypin%siglpMcat(i)*mypin%invsubalpham(i)
!    END DO
!    !
!    nfueldiv = myCell%nfueldiv
!    idxfuelsrt = FxrIdxSt+myCell%nnonfuel
!    IF (.NOT.myPin%lfuel .OR. myCell%lrect) CYCLE
!    CALL EffMacXsPSM(PSMXsMac(tid), Fxr(idxfuelsrt:FxrIdxSt+nLocalFXR-1, iz), mypin, PinFuelTempAvgsq, nfueldiv, iresoGrp1, iresoGrp2, ng, .TRUE.,iPSMcel)
!    
!    DO j = 1, nLocalFxr
!      ifxr = FxrIdxSt + j -1
!      myFxr => Fxr(ifxr, iz)
!      IF (idxfuelsrt.GT.ifxr) CYCLE ! beyond the fuel ring !.OR.ifxr.GT.FxrIdxSt+nLocalFXR-1
!      IF( .NOT. myFxr%lfuel ) CYCLE
!      niso = myFxr%niso
!      myfxr%idiso_pastpsm(1:niso) = myfxr%idiso(1:niso)
!      do ig = iResoGrp1, iResoGrp2
!        do iso = 1, niso
!          myFXR%fresoAIso(iso,ig) = 1._8
!          myFXR%fresoFIso(iso,ig) = 1._8
!        enddo
!      enddo
!      IF (nTRACERCntl%lgamma) THEN
!        do ig = iResoGrp1, iResoGrp2
!            do iso = 1, niso
!              myFxr%fresocapIso(iso, ig) = 1._8
!            enddo
!        enddo
!      END IF
!      IF (nTracerCntl%lRST) THEN
!        do ig = iResoGrp1, iResoGrp2
!          do iso = 1, niso
!            myFXR%fresoSIso(iso,ig) = 1._8
!            myFXR%fresoSSIso(iso,ig) = 1._8
!            myFXR%fresoS1Iso(iso,ig) = 1._8
!          enddo
!        enddo
!      ENDIF
!      idiso =>  myFXR%idiso
!      ifuelring = nfueldiv - (ifxr - idxfuelsrt)
!      CALL MacXsBase(XsMac(tid), myFxr, iResoGrp1, iResoGrp2, ng, eigv, FALSE, TRUE)
!      !nuclides without pxs information (fresoXiso = 1._8)
!      DO iso = 1, niso
!        id = mapnucl(idiso(iso));    isodata => ldiso(id)
!        idpxs = isodata%pxsid
!        IF (idpxs .NE. 0) CYCLE
!        !IF (idpxs .EQ. 0 .OR. idiso(iso) .LT. 90000) THEN
!        DO ig = iResoGrp1, iResoGrp2
!          PSMXsMac(tid)%IsoMGXsMacA(iso,ig,ifuelring) = XsMac(tid)%IsoXsMacA(iso,ig)
!          IF (isoXsMacFold(iso) .GT. epsm8) THEN
!            PSMXsMac(tid)%IsoMGXsMacF(iso,ig,ifuelring) = XsMac(tid)%IsoXsMacF(iso,ig)
!          END IF
!        END DO ! IG
!        !END IF
!      END DO ! ISO
!      ! Calculate freso
!      DO ig = iResoGrp1, iResoGrp2
!        CALL CP_VA(isoXsMacFold(1:niso),XsMac(tid)%IsoXsMacF(1:niso,ig),niso)
!        CALL CP_VA(isoXsMacNFold(1:niso),XsMac(tid)%IsoXsMacNF(1:niso,ig),niso)
!        CALL CP_VA(isoXsMackFold(1:niso),XsMac(tid)%IsoXsMackF(1:niso,ig),niso)
!        CALL CP_VA(isoXsMacAold(1:niso),XsMac(tid)%IsoXsMacA(1:niso,ig),niso)
!        IF (nTracerCntl%lRST) THEN
!          CALL CP_VA(isoXsMacSold(1:niso),XsMac(tid)%IsoXsMacS0(1:niso,ig),niso)
!          CALL CP_VA(isoXsMacS1old(1:niso),XsMac(tid)%IsoXsMacS1(1:niso,ig),niso)
!          CALL CP_VA(isoXsMacSSold(1:niso),XsMac(tid)%IsoXsMacSS(1:niso,ig),niso)
!        ENDIF
!        CALL CP_VA(IsoXsMacCapOld(1:niso),XsMac(tid)%IsoXsRadCap(1:niso,ig),niso)     !-- JSU EDIT 20170727
!        !
!        DO iso = 1, niso
!          myFXR%fresoAIso(iso,ig) = PSMXsMac(tid)%IsoMGXsMacA(iso,ig,ifuelring) / isoXsMacAold(iso)
!          IF (isoXsMacFold(iso) .GT. epsm8)   myFXR%fresoFIso(iso,ig) = PSMXsMac(tid)%IsoMGXsMacF(iso,ig,ifuelring) / isoXsMacFold(iso)
!        END DO
!        
!        IF (nTRACERCntl%lGamma) THEN
!          IF(isoXsMacFold(iso) .gt. epsm8) THEN ! Radioactive Capture resonance fraction...
!            DO iso = 1, niso
!              myFxr%fresocapIso(iso, ig) = (PSMXsMac(tid)%IsoMGXsMacA(iso,ig,ifuelring) - PSMXsMac(tid)%IsoMGXsMacF(iso,ig,ifuelring)) / (isoXsMacAold(iso) - isoXsMacFold(iso))
!            END DO
!          ELSE
!            DO iso = 1, niso
!              myFxr%fresocapIso(iso, ig) = PSMXsMac(tid)%IsoMGXsMacA(iso,ig,ifuelring) / isoXsMacAold(iso)
!            END DO
!          END IF
!        END IF
!        
!        IF (nTracerCntl%lRST) THEN
!          DO iso = 1, niso
!            if (abs(isoXsMacSold(iso)).GT. epsm8) myFXR%fresoSIso(iso,ig) = PSMXsMac(tid)%IsoMGXsMacS(iso,ig,ifuelring)/isoXsMacSold(iso)
!          ENDDO
!        ENDIF
!
!      ENDDO ! ig
!    ENDDO ! fxr
!  ENDDO ! ipin
!!  STOP 'DEBUGGING ** SUBGROUPEFFXS LINE 231'
!  !$OMP END DO
!  !$OMP END PARALLEL
!ENDDO ! iz
!TimeEnd = nTracer_dclock(FALSE,FALSE)
!TimeChk%PSMTime = TimeChk%PSMTime + TimeEnd - TimeBeg 
!CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
!
!NULLIFY(PIN, ResVarPin, CellINfo, myFxr, myPin, myCell)
!END SUBROUTINE
!
!SUBROUTINE GetfresoFXR(Core, Fxr, THInfo, eigv, GroupInfo, nTracerCntl, PE)
!! Subroutine to get freso of each FXR's macroscopic XSs based on fresoISO
!USE PARAM
!USE PointXSRT_MOD
!USE TYPEDEF,        ONLY : CoreInfo_Type,       FxrInfo_Type,          GroupInfo_Type,     &
!                           PE_Type,             THInfo_Type,           pin_type,           &
!                           ResVarpin_type,      cell_type,             XsMac_Type
!USE CNTL,           ONLY : nTracerCntl_Type
!!USE MacXsLib_mod,   ONLY : MacXsBase
!USE BasicOperation, ONLY : CP_VA
!USE IOUTIL,         ONLY : message
!USE FILES,          ONLY : io8
!USE TH_Mod,         ONLY : GetPinFuelTemp
!USE XSLIB_MOD,      ONLY : libdata, ldiso, mapnucl
!USE MPIComm_Mod,    ONLY : MPI_SYNC
!USE OMP_LIB
!USE TIMER,            ONLY : nTracer_dclock, TimeChk
!!
!IMPLICIT NONE
!!
!TYPE(CoreInfo_Type) :: Core
!TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
!TYPE(THInfo_Type) :: THInfo
!REAL :: eigv
!TYPE(GroupInfo_TYPE) :: GroupInfo
!TYPE(nTracerCntl_Type) :: nTracerCntl
!TYPE(PE_Type) :: PE
!!
!TYPE(Pin_Type), POINTER :: PIN(:)
!TYPE(ResVarPin_Type), POINTER :: ResVarPin(:,:),mypin
!TYPE(Cell_Type), POINTER :: CellInfo(:), myCell
!TYPE(FxrInfo_Type), POINTER :: myFxr
!TYPE(XsMac_Type), SAVE :: XsMac(nThreadMax)
!TYPE(PSMMacXs_Type), SAVE :: PSMXsMac(nThreadMax)
!!
!INTEGER :: nxy, myzb, myze, iz, ipin, xyb, xye, FxrIdxSt, nlocalFxr, j, ifxr, icel  !--- (xyb, xye)CNJ Edit : Domain Decomposition + MPI
!INTEGER :: ng, ig, iresoGrp1, iresoGrp2, niso, iso, tid
!REAL :: PinFuelTempAvgsq, XsMacAold, XsMacFold, XsMacNFold, XsMackFold, XsMacSold, XsMacStrold
!REAL :: XsMacAnew, XsMacFnew, XsMacNFnew, XsMackFnew, XsMacSnew,XsMacStrnew
!REAL :: isoXsMacAold(500),isoXsMacFold(500),isoXsMacNFold(500),isoXsMackFold(500),isoXsMacSold(500),isoXsMacS1old(500),isoXsMacSSold(500),isoxsmaccapold(500)
!LOGICAL :: lxslib, lAIC, lmcxs
!INTEGER :: nfueldiv, iPSMcel, id, i, ifuelring, idpxs, idxfuelsrt
!TYPE(LIBDATA), POINTER :: isodata
!REAL :: TimeEnd, TimeBeg 
!INTEGER, POINTER :: idiso(:), idiso_psm(:)
!INTEGER :: idxiso_psm(1000)
!INTEGER :: jso
!!
!IF (.not.nTracerCntl%lrestrmt) RETURN
!lxslib  = nTracerCntl%lxslib
!IF(.NOT. lxslib) RETURN
!!
!WRITE(mesg, '(A)') 'Update FXR freso by PSM...'
!IF(PE%Master) CALL message(io8, TRUE, TRUE, mesg)
!TimeBeg = nTracer_dclock(FALSE,FALSE)
!Pin => Core%Pin; CellInfo => Core%CellInfo
!ResVarPin => Core%ResVarPin;
!nxy = Core%nxy;  myzb = PE%myzb; myze = PE%myze
!
!!--- CNJ Edit : Domain Decomposition + MPI
!xyb = PE%myPinBeg; xye = PE%myPinEnd
!IF (PE%RTMASTER) THEN
!  xyb = 1; xye = nxy
!ENDIF
!
!ng = GroupInfo%ng
!iresoGrp1 = GroupInfo%nofg + 1; iresoGrp2 = GroupInfo%nofg + GroupInfo%norg
!
!CALL UpdtCoreIsoInfoPSM(Core, Fxr, PE)
!!
!CALL omp_set_num_threads(PE%nThread)
!tid = 1
!DO iz = myzb, myze
!  IF (.NOT. Core%lFuelPlane(iz)) CYCLE
!  !$OMP PARALLEL DEFAULT(SHARED)      &
!  !$OMP PRIVATE(tid, ipin, mypin, FxrIdxSt, icel, myCell, nLocalFxr, j, ifxr, myfxr, niso, &
!  !$OMP         iso, id, isodata, idxfuelsrt, ig, idiso, idpxs, XsMacAold, XsMacFold, XsMacNFold, XsMackFold, XsMacSold, &
!  !$OMP         isoXsMacAold,isoXsMacFold,isoXsMacNFold,isoXsMackFold,isoXsMacSold, &
!  !$OMP         XsMacAnew, XsMacFnew, XsMacNFnew, XsMackFnew, XsMacSnew, idxiso_psm, idiso_psm, jso)
!  !$ tid = omp_get_thread_num()+1
!  !$OMP DO
!  DO ipin = xyb, xye
!    mypin => ResVarPin(ipin,iz)
!    FxrIdxSt = Pin(ipin)%FxrIdxSt
!    icel = Pin(ipin)%Cell(iz)
!    myCell => CellInfo(icel)
!    IF (.NOT.myCell%lPSMcel) CYCLE
!    nlocalFxr = CellInfo(icel)%nFxr
!    idxfuelsrt = FxrIdxSt+myCell%nnonfuel
!    DO j = 1, nLocalFxr
!      ifxr = FxrIdxSt + j -1
!      myFxr => Fxr(ifxr, iz)
!      IF( .NOT. myFxr%lfuel ) CYCLE
!      IF (idxfuelsrt.GT.ifxr) CYCLE ! beyond the fuel ring !.OR.ifxr.GT.FxrIdxSt+nLocalFXR-1
!      niso = myFxr%niso
!      ! initialize the macroscopic freso of each FXR in pin
!      do ig = iResoGrp1, iResoGrp2
!        myFxr%fresoa(ig) = 1._8
!        myFxr%fresof(ig) = 1._8
!        myFxr%fresoNf(ig) = 1._8
!        myFxr%fresokf(ig) = 1._8
!      enddo
!      IF (nTracerCntl%lRST) THEN
!        do ig = iResoGrp1, iResoGrp2
!          myFxr%fresos(ig) = 1._8
!        enddo
!      ENDIF
!      idiso =>  myFXR%idiso
!      idiso_psm => myFXR%idiso_pastpsm
!      DO iso = 1, niso
!        DO jso = 1, GroupInfo%ntiso
!          IF (idiso_psm(jso).EQ.idiso(iso)) EXIT
!        END DO
!        idxiso_psm(iso) = jso
!      END DO
!      CALL MacXsBase(XsMac(tid), myFxr, iResoGrp1, iResoGrp2, ng, eigv, FALSE, TRUE)
!      ! Calculate freso
!      DO ig = iResoGrp1, iResoGrp2
!        XsMacAold = XsMac(tid)%XsMacA(ig); XsMacFold = XsMac(tid)%XsMacF(ig);
!        XsMacNFold = XsMac(tid)%XsMacNF(ig); XsMackFold = XsMac(tid)%XsMackF(ig)
!        XsMacSold = XsMac(tid)%XsMacS(ig)
!        XsMacAnew   = XsMacAold; XsMacFnew   = XsMacFold; 
!        XsMacNFnew  = XsMacNFold; XsMackFnew  = XsMackFold; 
!        XsMacSnew   = XsMacSold;
!        CALL CP_VA(isoXsMacFold(1:niso),XsMac(tid)%IsoXsMacF(1:niso,ig),niso)
!        CALL CP_VA(isoXsMacNFold(1:niso),XsMac(tid)%IsoXsMacNF(1:niso,ig),niso)
!        CALL CP_VA(isoXsMackFold(1:niso),XsMac(tid)%IsoXsMackF(1:niso,ig),niso)
!        CALL CP_VA(isoXsMacAold(1:niso),XsMac(tid)%IsoXsMacA(1:niso,ig),niso)
!        IF (nTracerCntl%lRST) THEN
!          CALL CP_VA(isoXsMacSold(1:niso),XsMac(tid)%IsoXsMacS0(1:niso,ig),niso)
!        ENDIF
!        !
!        DO iso = 1, niso
!          id = mapnucl(idiso(iso));    isodata => ldiso(id)
!          idpxs = isodata%pxsid
!          IF (idpxs .EQ. 0) CYCLE
!          XsMacAnew   = XsMacAnew   + (myFxr%fresoAISO(idxiso_psm(iso),ig) - 1._8)*isoXsMacAold(iso)
!          XsMacFnew   = XsMacFnew   + (myFxr%fresoFISO(idxiso_psm(iso),ig) - 1._8)*isoXsMacFold(iso)
!          XsMacNFnew  = XsMacNFnew  + (myFxr%fresoFISO(idxiso_psm(iso),ig) - 1._8)*isoXsMacNFold(iso)
!          XsMackFnew  = XsMackFnew  + (myFxr%fresoFISO(idxiso_psm(iso),ig) - 1._8)*isoXsMacKFold(iso)
!        END DO
!        IF (nTracerCntl%lRST) THEN
!          DO iso = 1, niso
!            XsMacSnew   = XsMacSnew   + (myFxr%fresoSISO(idxiso_psm(iso),ig) - 1._8)*isoXsMacSold(iso)
!          END DO
!        END IF
!        !
!        myFxr%fresoa(ig) = XsMacAnew / XsMacAold
!        IF(XsMacFold .gt. epsm8) THEN
!          myFxr%fresof(ig) = XsMacFnew / XsMacFold
!          myFxr%fresonf(ig) = XsMacNFnew / XsMacNFold
!          myFxr%fresokf(ig) = XsMackFnew / XsMacKFold
!        END IF
!        IF (nTracerCntl%lRST) myFxr%fresos(ig) = XsMacSnew / XsMacSold
!      ENDDO ! ig
!    ENDDO ! fxr
!  ENDDO ! ipin
!!  STOP 'DEBUGGING ** SUBGROUPEFFXS LINE 231'
!  !$OMP END DO
!  !$OMP END PARALLEL
!ENDDO ! iz
!TimeEnd = nTracer_dclock(FALSE,FALSE)
!TimeChk%PSMfresoTime = TimeChk%PSMfresoTime + TimeEnd - TimeBeg 
!CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
!!DEALLOCATE(idxiso_psm)
!NULLIFY(PIN, ResVarPin, CellINfo, myFxr, myPin, myCell)
!END SUBROUTINE
END MODULE
  
