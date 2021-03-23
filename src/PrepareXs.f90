#include <defines.h>
SUBROUTINE PrepFxr(lXsLib, lfxr)
!PrePare
USE PARAM
USE TYPEDEF,         ONLY : coreinfo_type,  Pin_Type,    Cell_Type,             &
                            Mixture_Type
USE GEOM,            ONLY : Core,           Pin,         CellInfo,              &
                            CellInfo,       ng
USE CORE_mod,        ONLY : Fxr,            GroupInfo,   FmInfo      
USE Boron_mod,       ONLY : SetBoronCoolant
USE XeDyn_Mod,       ONLY : InitXeDynInfo
USE Material_mod,    ONLY : Mixture
USE MacXsLib_Mod,    ONLY : SetCriticalNumberDensity
USE Depl_mod,        ONLY : GetCoreHmMass
USE CNTL,            ONLY : nTracerCntl    
USE PE_MOD,          ONLY : PE
USE BasicOperation,  ONLY : CP_VA
USE ALLOCS
USE BenchXs,         ONLY : MacXsBen
USE XSLIB_MOD,       ONLY : ldiso, mapnucl, CoreResIsoUpdate, nreshel
USE PointXSRT_MOD,         ONLY : SetMGnPWInfo, nPGrid, CorePXSIsoUpdate, CalcPSM_ISOPIN
IMPLICIT NONE

TYPE(Mixture_Type), POINTER :: Mix       
INTEGER :: i, j, k, ii, id
INTEGER :: iz, ipin, icel, ireg, ifsr, ifxr, ixsreg, ifsrbeg, ig, imix, imix1, imixp
INTEGER :: iasytype
INTEGER :: myzb, myze      
INTEGER :: nCoreFxr, nCoreFsr, nxy, FsrIdxSt, FxrIdxSt, nFsrInFxr
REAL :: AREA
LOGICAL :: lXslib, lres, lXeDyn, lfxr, lflag
LOGICAL :: lfuelpin, lfuelplane, lfuelasy, lcladplane

!XS Library Variables
INTEGER :: niso, ntiso, nchi

nCoreFxr = Core%nCoreFxr; nxy = Core%nxy
myzb = PE%myzb; myze = PE%myze

CALL SetGroupInfo(GroupInfo, lXsLib)

! #ifdef MPI_ENV
! !Updae the lists of the resonance isotopes in the Active core
! IF(.NOT. PE%RTMASTER) THEN
!   IF(lXsLib) CALL MakeResIsoList()
!   RETURN
! ENDIF
! #endif

lXeDyn = nTracerCntl%lXeDyn

IF(lXsLib) CALL SetCriticalNumberDensity(300._8)

IF(lXsLib) lres = nTracerCntl%lrestrmt
AllocatE(Fxr(nCoreFxr, myzb:myze))
DO iz = myzb, myze
  lfuelplane=.FALSE.; lcladplane=.FALSE.
  DO ipin = 1, nxy
    lFuelPin=.FALSE.
    iasytype=Pin(ipin)%asytype
    FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
    icel = Pin(ipin)%Cell(iz)
    lflag = .true.
    DO j = 1, CellInfo(icel)%nFxr
      ixsreg = FxrIdxSt + j - 1
      nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
      ifsrbeg = CellInfo(icel)%MapFxr2FsrIdx(1, j)
      !Area Calculation
      area = 0
      DO i = 1, nFsrInFxr
        area = area + CellInfo(icel)%vol(CellInfo(icel)%MapFxr2FsrIdx(i, j))
      ENDDO
      Fxr(ixsreg, iz)%area = area
      Fxr(ixsreg, iz)%ipin = ipin
      
      imix = CellInfo(icel)%ireg(ifsrbeg)
      if(j.eq.1) imix1=imix
      IF( lfxr )THEN
          Fxr(ixsreg, iz)%imix = ixsreg !celtyp(it)%ireg(ir1)
          Fxr(ixsreg, iz)%lFuel = MacXsBen(ixsreg)%lfuel
          IF(MacXsBen(ixsreg)%lfuel)THEN
              lFuelPin=.TRUE.      
              lFuelPlane=.TRUE.
              IF(.NOT. Core%AsyInfo(iasytype)%lfuel) Core%AsyInfo(iasytype)%lfuel=.TRUE.
          ENDIF          
          !Fxr(ixsreg, iz)%lFuel
      ELSE
          Fxr(ixsreg, iz)%imix = imix !celtyp(it)%ireg(ir1)
          Fxr(ixsreg, iz)%lFuel = CellInfo(icel)%lFuel
      ENDIF
      Fxr(ixsreg, iz)%nFsrInFxr = nFsrInFxr
      Fxr(ixsreg, iz)%FsrIdxSt = FsrIdxSt + ifsrbeg - 1
      IF(.NOT. lXsLib) THEN
        IF(nTracerCntl%libtyp .EQ. 11) THEN 
          Fxr(ixsreg, iz)%temp = 306.6 + CKELVIN
          Fxr(ixsreg, iz)%Doptemp = 618.3 + CKELVIN
          Fxr(ixsreg, iz)%rho = 0.7125
        END IF
        CYCLE
      END IF
      !Only for XSlibrary 
      MIX => Mixture(imix)
      Fxr(ixsreg, iz)%lfuel = Mix%lfuel; Fxr(ixsreg, iz)%ldepl = Mix%ldepl
      Fxr(ixsreg, iz)%lCLD = Mix%lCLD
      Fxr(ixsreg, iz)%lAIC = Mix%lAIC
      IF(lres .and. Mixture(imix)%lres) then
          IF (CellInfo(icel)%Geom%lcircle) Fxr(ixsreg, iz)%lres = TRUE  !! Only for circular geometry
      ENDIF
      if (imix.eq.imix1.and.lflag) then
          Fxr(ixsreg, iz)%lRes = .FALSE.
      else
          lflag=.false.
      endif
      Fxr(ixsreg, iz)%lh2o = Mix%lh2o;  Fxr(ixsreg, iz)%temp = Mix%temp
        !Allocation Isotope list and thems compostion list
      niso = Mix%niso; nchi = GroupInfo%nchi; ntiso = GroupInfo%ntiso
      IF(Fxr(ixsreg, iz)%ldepl) THEN
        CALL Dmalloc(Fxr(ixsreg, iz)%idiso, ntiso)
        CALL Dmalloc(Fxr(ixsreg, iz)%pnum, ntiso)
        IF (Fxr(ixsreg, iz)%lres .OR. Fxr(ixsreg, iz)%lCrRes) THEN
          CALL Dmalloc(Fxr(ixsreg, iz)%idiso_pastpsm, ntiso)
          CALL Dmalloc(Fxr(ixsreg, iz)%idx_Res, ntiso)
          CALL Dmalloc(Fxr(ixsreg, iz)%idiso_Res, nreshel)
          CALL Dmalloc(Fxr(ixsreg, iz)%pnum_Res, nreshel)
        END IF
        !CALL Dmalloc(Fxr(ixsreg, iz)%pnumrat, ntiso, ntiso)
        CALL Dmalloc(Fxr(ixsreg, iz)%chi, nchi)
        Fxr(ixsreg, iz)%ndim = ntiso
      ELSEIF(Fxr(ixsreg, iz)%lh2o .AND. nTracerCntl%lInitBoron) THEN
        CALL Dmalloc(Fxr(ixsreg, iz)%idiso, niso + 4)
        CALL Dmalloc(Fxr(ixsreg, iz)%pnum, niso + 4)
        IF (Fxr(ixsreg, iz)%lres) THEN
          CALL Dmalloc(Fxr(ixsreg, iz)%idiso_pastpsm, niso+4)
          CALL Dmalloc(Fxr(ixsreg, iz)%idx_Res, niso+4)
          CALL Dmalloc(Fxr(ixsreg, iz)%idiso_Res, nreshel)
          CALL Dmalloc(Fxr(ixsreg, iz)%pnum_Res, niso+4)
        END IF
        Fxr(ixsreg, iz)%ndim = niso + 4
      ELSEIF(Fxr(ixsreg, iz)%lh2o) THEN
        CALL Dmalloc(Fxr(ixsreg, iz)%idiso, niso + 4)
        CALL Dmalloc(Fxr(ixsreg, iz)%pnum, niso + 4)
        IF (Fxr(ixsreg, iz)%lres) THEN
          CALL Dmalloc(Fxr(ixsreg, iz)%idiso_pastpsm, niso+4)
          CALL Dmalloc(Fxr(ixsreg, iz)%idx_Res, niso+4)
          CALL Dmalloc(Fxr(ixsreg, iz)%idiso_Res, nreshel)
          CALL Dmalloc(Fxr(ixsreg, iz)%pnum_Res, niso+4)
        END IF
        Fxr(ixsreg, iz)%ndim = niso + 4
      ELSE
        CALL Dmalloc(Fxr(ixsreg, iz)%idiso, niso)
        CALL Dmalloc(Fxr(ixsreg, iz)%pnum, niso)
        IF (Fxr(ixsreg, iz)%lres) THEN
          CALL Dmalloc(Fxr(ixsreg, iz)%idiso_pastpsm, niso)
          CALL Dmalloc(Fxr(ixsreg, iz)%idx_Res, niso)
          CALL Dmalloc(Fxr(ixsreg, iz)%idiso_Res, nreshel)
          CALL Dmalloc(Fxr(ixsreg, iz)%pnum_Res, niso)
        END IF
        !CALL Dmalloc(Fxr(ixsreg, iz)%pnumrat, niso, niso)
        Fxr(ixsreg, iz)%ndim = niso
      ENDIF
      Fxr(ixsreg, iz)%niso = niso; Fxr(ixsreg, iz)%niso_depl = niso
      CALL CP_VA(Fxr(ixsreg,iz)%idiso(1:niso),mix%idiso(1:niso),niso)
      CALL CP_VA(Fxr(ixsreg,iz)%pnum(1:niso),mix%pnum(1:niso),niso)

      CALL FIndGdFxr(Fxr(ixsreg, iz))
      IF(Fxr(ixsreg, iz)%lDepl) CALL GetFxrHmMass(Fxr(ixsreg,iz), Core%hz(iz))
      if (imix.eq.imix1.and.lflag) then
          Fxr(ixsreg, iz)%lRes = .FALSE.
      else
          lflag=.false.
      endif
      IF (Fxr(ixsreg, iz)%lRes .AND. Fxr(ixsreg, iz)%lCLD) lcladplane=.TRUE.
#ifndef MPI_ENV
      !Update the lists of the resonance isotopes in the Active core
      IF(Fxr(ixsreg, iz)%lRes) CALL CoreResIsoUpdate(Fxr(ixsreg,iz)%idiso(1:niso), niso, iz)
#endif
      NULLIFY(MIX)
    ENDDO    ! IFXR in Cell
#ifndef MPI_ENV
    IF (nTRACERCntl%lPointRT) THEN
      DO j = 1, CellInfo(icel)%nFxr
        ixsreg = FxrIdxSt + j - 1
        nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
        ifsrbeg = CellInfo(icel)%MapFxr2FsrIdx(1, j)
        niso = Fxr(ixsreg,iz)%niso
        CALL CorePXSIsoUpdate(Fxr(ixsreg,iz)%idiso(1:niso), niso, iz)
      END DO  
    END IF
#endif
    !Pin(ipin)%lFuel=lFuelPin
    IF( lFXR )    Pin(ipin)%lFuel=lFuelPin
    ! Define Cladding region
    IF ((.not.lres).OR.(.NOT. lXsLib)) CYCLE
    lflag=.false.
    DO j = CellInfo(icel)%nFxr, 1, -1
      ixsreg = FxrIdxSt + j - 1
      ifsrbeg = CellInfo(icel)%MapFxr2FsrIdx(1, j)

      imix = CellInfo(icel)%ireg(ifsrbeg)
      if (j.eq.CellInfo(icel)%nFxr) then
          imixp = imix
          cycle
      endif
      if (imix.ne.imixp) then
        lflag=.true.
      endif
      imixp = imix
      MIX => Mixture(imix)
      IF (.not.Fxr(ixsreg, iz)%lRes) cycle
      IF (lflag) then
        do ii=1,Fxr(ixsreg,iz)%niso
          id=Fxr(ixsreg,iz)%idiso(ii)
          if (ldiso(mapnucl(id))%lclad) exit
        enddo
        if (ii.le.Fxr(ixsreg,iz)%niso) then
            Fxr(ixsreg, iz)%lCLD=.true.
            lcladplane=.TRUE.
        else
            Fxr(ixsreg, iz)%lCLD=.false.
            if (.not.Fxr(ixsreg, iz)%lfuel) Fxr(ixsreg, iz)%lres=.false.
        endif
      ENDIF
    ENDDO
  ENDDO !End of Pin Sweep
  !Core%lFuelPlane(iz)=lFuelPlane
  IF( lFXR ) Core%lFuelPlane(iz)=lFuelPlane
  IF(lXsLib) Core%lcladplane(iz)=lcladplane
ENDDO !End of Plane Sweep
DO iz = myzb, myze
  DO ipin = 1, nxy
    FxrIdxSt = Pin(ipin)%FxrIdxSt
    icel = Pin(ipin)%Cell(iz)
    DO j = 1, CellInfo(icel)%nFxr
      ixsreg = FxrIdxSt + j - 1
      if (Fxr(ixsreg, iz)%lres) then
        if (.not.Core%lFuelPlane(iz)) Fxr(ixsreg, iz)%lres=.FALSE.
      endif
    ENDDO
  ENDDO
ENDDO

FMInfo%Fxr => Fxr
#ifdef MPI_ENV
!Update the lists of the resonance isotopes in the Active core
IF(lXsLib) CALL MakeResIsoList()
IF(nTRACERCntl%lPointRT) CALL MakePXSIsoList()
#endif
IF(lXeDyn) THEN
  CALL InitXeDynInfo(Fxr, nCoreFxr, myzb, myze)
ENDIF

IF(nTracerCntl%lInitBoron) THEN
  CALL SetBoronCoolant(Core, Fxr, nTracerCntl%BoronPPM , myzb, myze)
ENDIF
IF (.NOT. nTracerCntl%lnTIGRst) THEN
  IF(nTracerCntl%lXsLib) THEN
    Core%Hm_Mass0 = GetCoreHmMass(Core, FmInfo%Fxr, PE)
  ENDIF
END IF

IF (nTRACERCntl%lPointRT) CALL SetMGnPWInfo
IF (nTRACERCntl%lPSM) CALL CalcPSM_ISOPIN

!CALL  Calnum(Core, FmInfo, PE)
!CALL  ReadFxr(Core, FmInfo, PE)
!CALL SetStructure(Core, Fxr, 0.0859375, myzb, myze)
IF(nTracerCntl%lProblem .EQ. lTransient) THEN 
  DO iz = myzb, myze
    DO ifxr = 1, Core%nCoreFxr
      FmInfo%Fxr(ifxr, iz)%imix0 = FmInfo%Fxr(ifxr, iz)%imix
    ENDDO
  ENDDO
END IF
END SUBROUTINE

SUBROUTINE FindGdFxr(myFxr)
USE PARAM
USE TYPEDEF,    ONLY : FxrInfo_Type
INTEGER :: iso, idiso

TYPE(FxrInfo_Type) :: myFxr

niso = myFxr%niso
DO iso = 1, niso
  idiso = myFxr%idiso(iso)
  idiso = idiso/1000
  IF(idiso == 64) myFxr%lGd = .TRUE.
ENDDO
END SUBROUTINE

SUBROUTINE SetStructure(Core, Fxr, spaden, myzb, myze)
USE PARAM
USE TypeDef,       ONLY : CoreInfo_type, Fxrinfo_type
USE NuclidMap_mod, ONLY : AtomicWeight

IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
REAL :: spaden
INTEGER :: myzb, myze

TYPE(FxrInfo_Type), POINTER :: myFxr

INTEGER :: nFxr
INTEGER :: iz, ifxr, id
INTEGER :: i, j, k
REAL :: rho, aw, tp1,tp2
nFxr = Core%nCoreFxr

DO iz = myzb, myze
  IF(.NOT. Core%lFuelPlane(iz)) CYCLE
  DO ifxr = 1, nFxr
    myFxr => Fxr(ifxr, iz)
    IF(.NOT. myFxr%lH2o) CYCLE
    aw = AtomicWeight(40002)
    rho = spaden * avogadro / aw
    myFxr%niso=myFxr%niso + 1
    myFxr%idiso(myFxr%niso) = 40002
    myFxr%pnum(myFxr%niso) = rho
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE MakeResIsoList()
USE PARAM
USE TYPEDEF,         ONLY : coreinfo_type,  Pin_Type,    Cell_Type,             &
                            Mixture_Type
USE GEOM,            ONLY : Core,           Pin,         CellInfo,              &
                            CellInfo,       ng
USE Material_mod,    ONLY : Mixture
USE CNTL,            ONLY : nTracerCntl    
USE PE_MOD,          ONLY : PE
USE xslib_mod,       ONLY : CoreResIsoUpdate
IMPLICIT NONE
TYPE(Mixture_Type), POINTER :: Mix      
INTEGER :: j, iz, ipin, icel, imix, ifsrbeg
INTEGER :: nxy, niso, myzb, myze

myzb = PE%myzb; myze = PE%myze
nxy = Core%nxy

DO iz = myzb, myze
  DO ipin = 1, nxy
    icel = Pin(ipin)%Cell(iz)
    DO j = 1, CellInfo(icel)%nFxr
      ifsrbeg = CellInfo(icel)%MapFxr2FsrIdx(1, j)
      imix = CellInfo(icel)%ireg(ifsrbeg)
      MIX => Mixture(imix)
      niso = MIX%niso
      IF(MIX%lres) CALL CoreResIsoUpdate(mix%idiso(1:niso), niso, iz)
    ENDDO
  ENDDO
ENDDO
NULLIFY(Mix)
END SUBROUTINE
    
SUBROUTINE MakePXSIsoList()
USE PARAM
USE TYPEDEF,         ONLY : coreinfo_type,  Pin_Type,    Cell_Type,             &
                            Mixture_Type
USE GEOM,            ONLY : Core,           Pin,         CellInfo,              &
                            CellInfo,       ng
USE Material_mod,    ONLY : Mixture
USE CNTL,            ONLY : nTracerCntl    
USE PE_MOD,          ONLY : PE
USE PointXSRT_MOD,         ONLY : CorePXSIsoUpdate
IMPLICIT NONE
TYPE(Mixture_Type), POINTER :: Mix      
INTEGER :: j, iz, ipin, icel, imix, ifsrbeg
INTEGER :: nxy, niso, myzb, myze

myzb = PE%myzb; myze = PE%myze
nxy = Core%nxy

DO iz = myzb, myze
  DO ipin = 1, nxy
    icel = Pin(ipin)%Cell(iz)
    DO j = 1, CellInfo(icel)%nFxr
      ifsrbeg = CellInfo(icel)%MapFxr2FsrIdx(1, j)
      imix = CellInfo(icel)%ireg(ifsrbeg)
      MIX => Mixture(imix)
      niso = MIX%niso
      IF(MIX%lres) CALL CorePXSIsoUpdate(mix%idiso(1:niso), niso, iz)
    ENDDO
  ENDDO
ENDDO
NULLIFY(Mix)
END SUBROUTINE
    
SUBROUTINE AllocResIsoInfo()
! RESONANT NUCLIDE �� �����ϴ� index ���� �Ҵ� �� �ʱ�ȭ + category �ʱ�ȭ
USE PE_MOD,        ONLY : PE
USE CNTL,          ONLY : nTracerCntl    
USE XSLIB_MOD,     ONLY : nreshel, maxnid, nActiveCat, IDRES_COR, NRES_COR, mapnuclRes, mlgdata, npwxs
USE ALLOCS
USE nuclidmap_mod, ONLY : InitResoCat
USE PointXSRT_MOD,       ONLY : mapnuclPXS, IDPXS_COR, NPXS_COR
IMPLICIT NONE

INTEGER :: myzb, myze, iz, i

myzb = PE%myzb; myze = PE%myze
CALL Dmalloc0(NRES_COR, myzb, myze)
CALL Dmalloc0(IDRES_COR, 1, nreshel, myzb, myze)
CALL Dmalloc0(mapnuclRes, 1, maxnid, myzb, myze)
ALLOCATE(mlgdata(myzb:myze))
DO iz=myzb,myze
  NRES_COR(iz)=0
  DO i=1,nreshel
    IDRES_COR(i,iz)=0
  ENDDO
  DO i=1,maxnid
    mapnuclRes(i,iz)=0
  ENDDO
ENDDO

IF (.not.nTracerCntl%lMLG) THEN
  IF (nTracerCntl%lCAT) THEN
    CALL Dmalloc0(nActiveCat, myzb, myze)
    DO iz=myzb,myze
      nActiveCat(iz)=0
    ENDDO
    CALL InitResoCat()
  ENDIF
ENDIF
!
IF (nTRACERCntl%lPointRT) THEN
  CALL Dmalloc0(NPXS_COR, myzb, myze)
  CALL Dmalloc0(IDPXS_COR, 1, npwxs, myzb, myze)
  CALL Dmalloc0(mapnuclPXS, 1, maxnid, myzb, myze)
  DO iz=myzb,myze
    NPXS_COR(iz)=0
    DO i=1,npwxs
      IDPXS_COR(i,iz)=0
    ENDDO
    DO i=1,maxnid
      mapnuclPXS(i,iz)=0
    ENDDO
  ENDDO
END IF
!
END SUBROUTINE
    
SUBROUTINE AllocXsEQ()
!PrePare variables related to resonance treatment using equivalent XS
USE PARAM
USE TYPEDEF,       ONLY : Cell_TYpe, Pin_Type, ResVarPin_Type
USE GEOM,          ONLY : Core
USE CORE_mod,      ONLY : Fxr,            GroupInfo
USE CNTL,          ONLY : nTracerCntl    
USE PE_MOD,        ONLY : PE
USE BasicOperation,ONLY : CP_CA
USE XSLIB_MOD,     ONLY : nlvflxmax, nreshel, nelthel, mlgdata0, nCat
USE ALLOCS
USE CP_mod,        ONLY : ngauss
USE PointXSRT_MOD, ONLY : nPGrid
IMPLICIT NONE
TYPE(Cell_TYpe), POINTER :: CellInfo(:)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(ResVarPin_Type), POINTER :: ResVarPin(:,:)
INTEGER :: nFxr, iResGrpBeg, iResGrpEnd
INTEGER :: myzb, myze
INTEGER :: iz, ifxr, icel, ipin, nxy, FxrIdxSt, j, ig, igt, ngt, nlocalFxr, ir
REAL :: fuelrad(100), sumvol
  
IF (.not.nTracerCntl%lrestrmt) RETURN

nFxr = Core%nCoreFXR
iResGrpBeg = GroupInfo%nofg + 1;
iResGrpEnd = GroupInfo%nofg + GroupInfo%norg
myzb = PE%myzb; myze = PE%myze

nxy = Core%nxy
Pin => Core%Pin; CellInfo => Core%CellInfo
ALLOCATE(Core%ResVarPin(1:nxy,myzb:myze))

ResVarPin => Core%ResVarPin
DO iz = myzb, myze
  if (.not.Core%lFuelPlane(iz)) CYCLE  
  DO ipin = 1, nxy
    FxrIdxSt = Pin(ipin)%FxrIdxSt
    icel = Pin(ipin)%Cell(iz)
    DO j = 1, CellInfo(icel)%nFxr
        ifxr = FxrIdxSt + j - 1
        IF (Fxr(ifxr, iz)%lRes) then
            ResVarPin(ipin,iz)%lresA=.TRUE. ! unused..?
            IF (.not.Fxr(ifxr, iz)%lCLD) then
                ResVarPin(ipin,iz)%lres=.TRUE. ! existence of fuel or AIC?
            ELSE
                ResVarPin(ipin,iz)%lresC=.TRUE. ! unused..?
            ENDIF
        ELSE
          IF(Fxr(ifxr, iz)%lCrRes) THEN 
            ResVarPin(ipin, iz)%lCrRes = .TRUE.
          END IF
        ENDIF
        IF (FXR(ifxr,iz)%lfuel) ResVarPin(ipin, iz)%lfuel = .TRUE.
    ENDDO
  ENDDO
ENDDO

IF (nTracerCntl%lED) THEN
    IF (nTracerCntl%lMLG) CALL SetMLGrid_ED
ELSE
    IF (nTracerCntl%lMLG) CALL SetMLGrid
ENDIF

DO iz = myzb, myze    
  if (.not.Core%lFuelPlane(iz)) CYCLE
  DO ifxr = 1, nFxr
    IF (.NOT. (Fxr(ifxr,iz)%lres .OR. Fxr(ifxr, iz)%lCRRes)) CYCLE
    CALL Dmalloc0(Fxr(ifxr, iz)%fresoa, iResGrpBeg, iResGrpEnd)
    CALL Dmalloc0(Fxr(ifxr, iz)%fresof, iResGrpBeg, iResGrpEnd)    
    CALL Dmalloc0(Fxr(ifxr, iz)%fresonf, iResGrpBeg, iResGrpEnd)    
    CALL Dmalloc0(Fxr(ifxr, iz)%fresokf, iResGrpBeg, iResGrpEnd)    
    do ig = iResGrpBeg,iResGrpEnd
      Fxr(ifxr, iz)%fresoa(ig)=1.
      Fxr(ifxr, iz)%fresof(ig)=1.
      Fxr(ifxr, iz)%fresonf(ig)=1.
      Fxr(ifxr, iz)%fresokf(ig)=1.
    enddo
    
    ALLOCATE(Fxr(ifxr, iz)%fresoAIso(Fxr(ifxr, iz)%ndim,iResGrpBeg:iResGrpEnd)) 
    ALLOCATE(Fxr(ifxr, iz)%fresoFIso(Fxr(ifxr, iz)%ndim,iResGrpBeg:iResGrpEnd)) 
    do ig = iResGrpBeg,iResGrpEnd
      Fxr(ifxr, iz)%fresoAIso(:,ig)=1.
      Fxr(ifxr, iz)%fresoFIso(:,ig)=1.
    ENDDO
    IF(nTracerCntl%lGamma) THEN
      ALLOCATE(Fxr(ifxr, iz)%fresoCapIso(Fxr(ifxr, iz)%ndim,iResGrpBeg:iResGrpEnd)) 
      Fxr(ifxr, iz)%fresoCapIso = 1.
    END IF
    
    IF (nTracerCntl%lRST) THEN
      CALL Dmalloc0(Fxr(ifxr, iz)%fresos, iResGrpBeg, iResGrpEnd)
      CALL Dmalloc0(Fxr(ifxr, iz)%fresostr, iResGrpBeg, iResGrpEnd)
      do ig = iResGrpBeg,iResGrpEnd
        Fxr(ifxr, iz)%fresos(ig)=1.
        Fxr(ifxr, iz)%fresostr(ig)=1.
      enddo
    
      ALLOCATE(Fxr(ifxr, iz)%fresoSIso(Fxr(ifxr, iz)%ndim,iResGrpBeg:iResGrpEnd)) 
      ALLOCATE(Fxr(ifxr, iz)%fresoSSIso(Fxr(ifxr, iz)%ndim,iResGrpBeg:iResGrpEnd)) 
      ALLOCATE(Fxr(ifxr, iz)%fresoS1Iso(Fxr(ifxr, iz)%ndim,iResGrpBeg:iResGrpEnd)) 
      do ig = iResGrpBeg,iResGrpEnd
        Fxr(ifxr, iz)%fresoSIso(:,ig)=1.
        Fxr(ifxr, iz)%fresoSSIso(:,ig)=1.
        Fxr(ifxr, iz)%fresoS1Iso(:,ig)=1.
      ENDDO       
    ENDIF
    IF (nTracerCntl%lMLG) THEN
      IF (Fxr(ifxr,iz)%lCLD) THEN
        CALL Dmalloc(Fxr(ifxr, iz)%XsEq_c_1g, mlgdata0%c_nmaclv1G)
        Fxr(ifxr, iz)%XsEq_c_1g=0._8
      ELSEIF (Fxr(ifxr,iz)%lAIC .OR. Fxr(ifxr, iz)%lCRRes) THEN
        CALL Dmalloc(Fxr(ifxr, iz)%XsEq_f_1g, mlgdata0%f_nmaclv1G)
        Fxr(ifxr, iz)%XsEq_f_1g=0._8  
      ELSE
      ! fuel
        CALL Dmalloc0(Fxr(ifxr, iz)%XsEq_f_mg, 1, mlgdata0%f_nmaclv, iResGrpBeg, iResGrpEnd)   
        Fxr(ifxr, iz)%XsEq_f_mg=0._8
        CALL Dmalloc0(Fxr(ifxr, iz)%FnAdj, iResGrpBeg, iResGrpEnd)  
        Fxr(ifxr, iz)%FnAdj=1._8
        CALL Dmalloc0(Fxr(ifxr, iz)%FtAdj, 1, mlgdata0%f_nmaclv, iResGrpBeg, iResGrpEnd)  
        Fxr(ifxr, iz)%FtAdj=1._8  
      ENDIF
    ELSE
      IF (nTracerCntl%lCAT) THEN
        CALL Dmalloc0(Fxr(ifxr, iz)%XsEq, 1, nlvflxmax, 1, nCat, iResGrpBeg, iResGrpEnd)   
        CALL Dmalloc0(Fxr(ifxr, iz)%NDAF, 1, nlvflxmax, 1, nCat, iResGrpBeg, iResGrpEnd)    
      ELSE
        CALL Dmalloc0(Fxr(ifxr, iz)%XsEq, 1, nlvflxmax, 1, nreshel, iResGrpBeg, iResGrpEnd)   
        CALL Dmalloc0(Fxr(ifxr, iz)%NDAF, 1, nlvflxmax, 1, nreshel, iResGrpBeg, iResGrpEnd)    
      ENDIF
      Fxr(ifxr, iz)%XsEq = 0._8
      Fxr(ifxr, iz)%NDAF = 1._8
    ENDIF
  ENDDO
ENDDO

DO iz = myzb, myze
  if (.not.Core%lFuelPlane(iz)) CYCLE
  DO ipin = 1, nxy
    if (.not. (ResVarPin(ipin,iz)%lres .OR. ResVarPin(ipin, iz)%lCrRes)) CYCLE  
    icel = Pin(ipin)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
    ALLOCATE(ResVarPin(ipin,iz)%rad_cp(nlocalFxr))
    ALLOCATE(ResVarPin(ipin,iz)%delr(nlocalFxr))
    ALLOCATE(ResVarPin(ipin,iz)%vol(nlocalFxr))
    ALLOCATE(ResVarPin(ipin,iz)%Qsurfvol(nlocalFxr))
    ALLOCATE(ResVarPin(ipin,iz)%X(nlocalFxr,nlocalFxr,ngauss))
    
    ResVarPin(ipin,iz)%rad_cp(1:nlocalFxr)=CellInfo(icel)%rad_cp(1:nlocalFxr)    
    CALL Dmalloc0(ResVarPin(ipin,iz)%idiso, 1, nelthel)  
    CALL Dmalloc0(ResVarPin(ipin,iz)%pnum, 1, nelthel)  
    CALL Dmalloc0(ResVarPin(ipin,iz)%idiso_Res, 1, nreshel)  
    CALL Dmalloc0(ResVarPin(ipin,iz)%idx_Res, 1, nelthel)  
    CALL Dmalloc0(ResVarPin(ipin,iz)%pnum_Res, 1, nreshel)  
    !CALL Dmalloc0(ResVarPin(ipin,iz)%pnumrat, 1, nreshel, 1, nreshel)  
    
    ResVarPin(ipin,iz)%lbar = 2._8 * CellInfo(icel)%FuelRad0 ! mean chord length
    ResVarPin(ipin,iz)%XSEsc = 1._8 / ResVarPin(ipin,iz)%lbar ! escape XS in Equivalence Theory
    
    IF (.not.nTracerCntl%lMLG) THEN
      IF (nTracerCntl%lCAT) THEN
        CALL Dmalloc0(ResVarPin(ipin,iz)%avgxseq, 1, nlvflxmax, 1, nCat, iResGrpBeg, iResGrpEnd)  
      ELSE
        CALL Dmalloc0(ResVarPin(ipin,iz)%avgxseq, 1, nlvflxmax, 1, nreshel, iResGrpBeg, iResGrpEnd)
      ENDIF
      ResVarPin(ipin,iz)%avgxseq=0._8
    ELSE    
      CALL Dmalloc0(ResVarPin(ipin,iz)%avgxseq_mg, 1, mlgdata0%f_nmaclv, iResGrpBeg, iResGrpEnd)  
      ResVarPin(ipin,iz)%avgxseq_mg=0._8
      CALL Dmalloc0(ResVarPin(ipin,iz)%avgxseq_1g, 1, mlgdata0%f_nmaclv1G)  
      ResVarPin(ipin,iz)%avgxseq_1g=0._8
      CALL Dmalloc0(ResVarPin(ipin,iz)%FnAdj, iResGrpBeg, iResGrpEnd)  
      ResVarPin(ipin,iz)%FnAdj=1._8
    ENDIF
    IF (nTracerCntl%lRIF) THEN
      CALL Dmalloc0(ResVarPin(ipin,iz)%rifa, 1, nreshel, iResGrpBeg, iResGrpEnd)
      CALL Dmalloc0(ResVarPin(ipin,iz)%riff, 1, nreshel, iResGrpBeg, iResGrpEnd)   
      IF (nTracerCntl%lRST) CALL Dmalloc0(ResVarPin(ipin,iz)%rifs, 1, nreshel, iResGrpBeg, iResGrpEnd)
    ENDIF
  ENDDO
ENDDO
!
IF (nTRACERCntl%lPSM) THEN
DO iz = myzb, myze
  if (.not.Core%lFuelPlane(iz)) CYCLE
  DO ipin = 1, nxy
    if (.not.(ResVarPin(ipin,iz)%lres .OR. ResVarPin(ipin,iz)%lcrres)) CYCLE  
    ! for PSM EDIT JSU 2020/06
    ALLOCATE(ResVarPin(ipin,iz)%eta(nPGrid))
  END DO ! ipin
END DO ! iz
END IF
!
END SUBROUTINE
SUBROUTINE SetMLGrid
  USE PARAM
  USE TYPEDEF,       ONLY : Cell_TYpe, Pin_Type
  USE GEOM,          ONLY : Core
  USE CORE_mod,      ONLY : Fxr
  USE CNTL,          ONLY : nTracerCntl    
  USE PE_MOD,        ONLY : PE
  USE XSLIB_MOD
  IMPLICIT NONE
  TYPE(Cell_TYpe), POINTER :: CellInfo(:)
  TYPE(Pin_Type), POINTER :: Pin(:)
  INTEGER :: nFxr
  INTEGER :: myzb, myze
  INTEGER :: iz, ifxr, icel, ipin, nxy, nmlg, il
  REAL :: radmin,radmax,del,rad,sigmax,sigmin
  
  myzb = PE%myzb; myze = PE%myze
  nxy = Core%nxy
  Pin => Core%Pin; CellInfo => Core%CellInfo

  nmlg=nTracerCntl%nMLG ! Default : 8 / modifiable by input
  mlgdata0%f_nmaclv=nmlg
  mlgdata0%f_nmaclv1G=nmlg
  
!  mlgdata0%c_nmaclv1G=5
  mlgdata0%c_nmaclv1G=nTracerCntl%nMLGc ! Default : 5 / fixed
  allocate(mlgdata0%c_maclv1G(mlgdata0%c_nmaclv1G))
  allocate(mlgdata0%c_maclv1G_log(mlgdata0%c_nmaclv1G))
  mlgdata0%c_maclv1G=(/0.1_8,1._8,10._8,100._8,1000._8/)  
  mlgdata0%c_maclv1G_log=dlog(mlgdata0%c_maclv1G)
  mlgdata0%del_c1g=dlog(10._8);
  ! Set MLG of fuel for each plane
  do iz=myzb,myze
    if (.not.core%lfuelplane(iz)) cycle
    mlgdata(iz)%f_nmaclv=nmlg
    allocate(mlgdata(iz)%f_maclv(nmlg)) 
    allocate(mlgdata(iz)%f_maclv_log(nmlg))
    radmin=100._8; radmax=0._8
    DO ipin = 1, nxy  
        icel = Pin(ipin)%Cell(iz)      
        if (.not.CellInfo(icel)%lfuel) cycle
        if (CellInfo(icel)%FuelRad0.le.0.2_8) cycle
        radmin=min(radmin,CellInfo(icel)%FuelRad0)
        radmax=max(radmax,CellInfo(icel)%FuelRad0)
    ENDDO
    sigmax=82._8/radmin
    sigmin=0.164_8/radmax
    del=(sigmax/sigmin)**(1._8/nmlg)
    mlgdata(iz)%f_maclv(1)=sigmin; mlgdata(iz)%f_maclv(nmlg)=sigmax
    mlgdata(iz)%f_maclv_log(1)=dlog(sigmin); mlgdata(iz)%f_maclv_log(nmlg)=dlog(sigmax)
    DO il=2,nmlg-1
        mlgdata(iz)%f_maclv(il)=mlgdata(iz)%f_maclv(il-1)*del
        mlgdata(iz)%f_maclv_log(il)=dlog(mlgdata(iz)%f_maclv(il))
    ENDDO
    mlgdata(iz)%del_fMG = dlog(del);
    if (.not.core%lAICplane(iz) .AND. .NOT. nTracerCntl%lCrInfo ) cycle
    ! Set MLG of AIC for each plane    
    mlgdata(iz)%f_nmaclv1G=nmlg
    allocate(mlgdata(iz)%f_maclv1G(nmlg))
    allocate(mlgdata(iz)%f_maclv1G_log(nmlg))
    radmin=100._8; radmax=0._8
    DO ipin = 1, nxy  
        icel = Pin(ipin)%Cell(iz)      
        
        IF(CellInfo(icel)%lAIC) THEN
          rad=CellInfo(icel)%FuelRad0
        ELSEIF(CellInfo(icel)%lCrCell .AND. CellInfo(CellInfo(icel)%CrCell%CellCrIdx)%lAIC) THEN 
          rad = CellInfo(CellInfo(icel)%CrCell%CellCrIdx)%FuelRad0
        ELSE
          CYCLE
        END IF
        if (rad.le.0.2_8) cycle
        radmin=min(radmin,rad)
        radmax=max(radmax,rad)
    ENDDO
    sigmax=82._8/radmin
    sigmin=0.164_8/radmax
    del=(sigmax/sigmin)**(1._8/nmlg)
    mlgdata(iz)%f_maclv1G(1)=sigmin; mlgdata(iz)%f_maclv1G(nmlg)=sigmax
    mlgdata(iz)%f_maclv1G_log(1)=dlog(sigmin); mlgdata(iz)%f_maclv1G_log(nmlg)=dlog(sigmax)
    mlgdata(iz)%del_f1G = dlog(del);
    DO il=2,nmlg-1
        mlgdata(iz)%f_maclv1G(il)=mlgdata(iz)%f_maclv1G(il-1)*del
        mlgdata(iz)%f_maclv1G_log(il)=dlog(mlgdata(iz)%f_maclv1G(il))
    ENDDO
    
  enddo
  
END SUBROUTINE

SUBROUTINE SetMLGrid_ED
  USE PARAM
  USE TYPEDEF,       ONLY : Cell_TYpe, Pin_Type, ResVarPin_Type
  USE GEOM,          ONLY : Core
  USE CORE_mod,      ONLY : Fxr
  USE CNTL,          ONLY : nTracerCntl    
  USE PE_MOD,        ONLY : PE
  USE XSLIB_MOD
  USE UtilFunction,  ONLY : Array1DSORT
  IMPLICIT NONE
  TYPE(Cell_TYpe), POINTER :: CellInfo(:)
  TYPE(Pin_Type), POINTER :: Pin(:)
  TYPE(ResVarPin_Type), POINTER :: ResVarPin(:,:)
  INTEGER :: nFxr
  INTEGER :: myzb, myze
  INTEGER :: iz, ifxr, icel, ipin, nxy, nmlg, il, igt, ngt
  REAL :: del,rad,sigmax,sigmin
  REAL :: fuelrad(100)
  
  myzb = PE%myzb; myze = PE%myze
  nxy = Core%nxy
  Pin => Core%Pin; CellInfo => Core%CellInfo
  ResVarPin => Core%ResVarPin

  nmlg=nTracerCntl%nMLG
  mlgdata0%f_nmaclv=nmlg
  mlgdata0%f_nmaclv1G=nmlg
  
!  mlgdata0%c_nmaclv1G=5
  mlgdata0%c_nmaclv1G=nTracerCntl%nMLGc
  allocate(mlgdata0%c_maclv1G(mlgdata0%c_nmaclv1G))
  allocate(mlgdata0%c_maclv1G_log(mlgdata0%c_nmaclv1G))
  mlgdata0%c_maclv1G=(/0.1_8,1._8,10._8,100._8,1000._8/)  
  mlgdata0%c_maclv1G_log=dlog(mlgdata0%c_maclv1G)
  mlgdata0%del_c1g = dlog(10._8);
  do iz=myzb,myze
    if (.not.core%lfuelplane(iz)) cycle
    
    mlgdata(iz)%f_nmaclv=nmlg
    mlgdata(iz)%f_nmaclv1G=nmlg
    ! count the number of fuel pellet sizes
    ngt = 0; fuelrad=0._8
    DO ipin = 1, nxy
      if (.NOT.ResVarPin(ipin,iz)%lres) cycle
      icel = Pin(ipin)%Cell(iz)
      do igt = 1, ngt
          if (CellInfo(icel)%fuelrad0.eq.fuelrad(igt)) exit
      enddo
      if (igt.gt.ngt) then
          ngt = ngt + 1
          fuelrad(ngt) = CellInfo(icel)%fuelrad0
      endif
    ENDDO
    CALL Array1DSORT(fuelrad(1:ngt), ngt, ngt, .true., .false.)
    DO ipin = 1, nxy
      IF (.not.ResVarPin(ipin,iz)%lres) CYCLE
      icel = Pin(ipin)%Cell(iz)
      DO igt = 1, ngt
          if (fuelrad(igt).eq.CellInfo(icel)%fuelrad0) exit
      ENDDO        
      ResVarPin(ipin,iz)%igt = igt
    ENDDO    
    mlgdata(iz)%ngeomtype = ngt
  
    allocate(mlgdata(iz)%f_maclv_pin(nmlg,ngt)) 
    allocate(mlgdata(iz)%f_maclv_pin_log(nmlg,ngt))
    DO igt = 1, ngt  
      rad = fuelrad(igt)
      sigmax=82._8/rad
      sigmin=0.164_8/rad  
      del=(sigmax/sigmin)**(1._8/nmlg)
      mlgdata(iz)%f_maclv_pin(1,igt)=sigmin; mlgdata(iz)%f_maclv_pin(nmlg,igt)=sigmax
      mlgdata(iz)%f_maclv_pin_log(1,igt)=dlog(sigmin)
      mlgdata(iz)%f_maclv_pin_log(nmlg,igt)=dlog(sigmax)
      DO il=2,nmlg-1
        mlgdata(iz)%f_maclv_pin(il,igt)=mlgdata(iz)%f_maclv_pin(il-1,igt)*del
        mlgdata(iz)%f_maclv_pin_log(il,igt)=dlog(mlgdata(iz)%f_maclv_pin(il,igt))
      ENDDO  
      mlgdata(iz)%del_fmg = dlog(del);
      mlgdata(iz)%del_f1g = dlog(del);
    ENDDO
  enddo
  
END SUBROUTINE

!Set Group Sweeping and Scattering Range Setting
SUBROUTINE SetGroupInfo(GroupInfo,lXsLib)
USE PARAM
USE TYPEDEF,           ONLY : GroupInfo_Type
!USE CORE_MOD,          ONLY : GroupInfo
USE GEOM,              ONLY : ng,              nPrec,         ngg
USE BenchXs,           ONLY : nXslType,        XsSmBen,       DynMacXsBen,  XssmDynBen
USE ALLOCS             
USE XsLib_Mod,         ONLY : nelthel,         noghel,        nelrhel,      nchihel,              &
                              ldiso ,          noggphl,       enbhel,       enbgam,  phatom, nelmGAM
USE TranMacXsLib_Mod,  ONLY : SetTranXsLibInfo
USE cntl,              ONLY : nTracerCntl
IMPLICIT NONE

TYPE(GroupInfo_Type) :: GroupInfo
LOGICAL :: lXsLib
REAL, POINTER :: XsMacSm(:, :)
INTEGER :: i, j, k, ig, ig2, iso
INTEGER :: nid, ntemp, itemp
INTEGER :: igg, igg2

ALLOCATE(XsMacSm(ng, ng))

IF(lXslib) ng = noghel
IF(.NOT. lXsLib) GroupInfo%nprec = nprec
IF(nTracerCntl%lProblem .EQ. lTransient .AND. lXsLib) THEN
  CALL SetTranXsLibInfo(GroupInfo)
  nprec = GroupInfo%nprec
ENDIF
GroupInfo%ng = ng
GroupInfo%iresoGrp1 = GroupInfo%nofg + 1; GroupInfo%iresoGrp2 = GroupInfo%nofg + GroupInfo%norg

CALL Dmalloc(GroupInfo%InScatRange, 2, ng)
CALL Dmalloc(GroupInfo%OutScatRange, 2, ng)
GroupInfo%InScatRange(1, :) = ng;  GroupInfo%InScatRange(2, :) = 1;
GroupInfo%OutScatRange(1, :) = ng; GroupInfo%OutScatRange(2, :) = 1;
!InScattering Info

IF(lXsLib) THEN
  GroupInfo%UpScatRange(:) = noghel
  GroupInfo%nChi = nChihel
  DO iso = 1, nelthel
    nid = iso; ntemp = ldiso(nid)%ntemp
    DO ig = 1, noghel
      DO itemp = 1, ntemp
        !Inscattering Range
        GroupInfo%InScatRange(1, ig) = min(GroupInfo%InScatRange(1, ig), ldiso(nid)%sm(ig, itemp)%ib)
        GroupInfo%InScatRange(2, ig) = max(GroupInfo%InScatRange(2, ig), ldiso(nid)%sm(ig, itemp)%ie)          
        !Out Scattering Range
        DO ig2 = ldiso(nid)%sm(ig, itemp)%ib, ldiso(nid)%sm(ig, itemp)%ie
          GroupInfo%OutScatRange(1, ig2) = min(ldiso(nid)%sm(ig, itemp)%ib, GroupInfo%OutScatRange(1, ig2))
          GroupInfo%OutScatRange(2, ig2) = max(ldiso(nid)%sm(ig, itemp)%ie, GroupInfo%OutScatRange(2, ig2))
        ENDDO
        !Upscatering Range
        IF(ldiso(nid)%sm(ig,itemp)%ie .GT. ig) THEN
          GroupInfo%UpScatRange(1) = min(ig, GroupInfo%UpScatRange(1))
          CONTINUE
        ENDIF
      ENDDO !Temperature Sweep
    ENDDO  !Group Sweep
  ENDDO  !Isotope Sweep
  
  IF (nTracerCntl%lGamma) THEN
    ngg = noggphl
    GroupInfo%ngg = ngg
    
    ALLOCATE(GroupInfo%InScatRange_ph(2, ngg))
    ALLOCATE(GroupInfo%OutScatRange_ph(2, ngg))
    GroupInfo%InScatRange_ph(1, :) = ngg;  GroupInfo%InScatRange_ph(2, :) = 1;
    GroupInfo%OutScatRange_ph(1, :) = ngg; GroupInfo%OutScatRange_ph(2, :) = 1;
    
    !InScattering Info
    GroupInfo%UpScatRange_ph = ngg
    DO nid = 1, nelmGAM
      DO igg = 1, ngg
        !Inscattering Range
        GroupInfo%InScatRange_ph(1, igg) = min(GroupInfo%InScatRange_ph(1, igg), phatom(nid)%sm(igg)%ib)
        GroupInfo%InScatRange_ph(2, igg) = max(GroupInfo%InScatRange_ph(2, igg), phatom(nid)%sm(igg)%ie)          
        !Out Scattering Range
        DO igg2 = phatom(nid)%sm(igg)%ib, phatom(nid)%sm(igg)%ie
          GroupInfo%OutScatRange_ph(1, igg2) = min(phatom(nid)%sm(igg)%ib, GroupInfo%OutScatRange_ph(1, igg2))
          GroupInfo%OutScatRange_ph(2, igg2) = max(phatom(nid)%sm(igg)%ie, GroupInfo%OutScatRange_ph(2, igg2))
        ENDDO
        !Upscatering Range
        IF(phatom(nid)%sm(igg)%ie .GT. igg) THEN
          GroupInfo%UpScatRange_ph(1) = min(igg, GroupInfo%UpScatRange_ph(1))
          CONTINUE
        ENDIF
      ENDDO  !Group Sweep
    ENDDO  !Isotope Sweep

    GroupInfo%UpScatRange_ph(2) = ngg
    GroupInfo%lUpScat_ph = FALSE
    DO igg = 1, ngg
      IF(GroupInfo%InScatRange_ph(2, igg) .GT. igg) THEN
        GroupInfo%lUpScat_ph = TRUE
        GroupInfo%UpScatRange_ph(1) = igg
        EXIT
      ENDIF
    ENDDO
    
  END IF
ELSE
  DO i = 1, nXslType
    IF(nTracerCntl%lDynamicBen) THEN
      ntemp = DynMacXsBen(i)%ntemp
      DO itemp = 1, ntemp
        CALL XsSmDynBen(i, DynMacXsBen(i)%temp(itemp), 1, ng, 1, ng, XsMacSm)
        !InScattering
        DO ig = 1, ng
          !Upper Bound
          DO ig2 = 1, ig
            IF(XsMacSm(ig2, ig) .NE. ZERO) EXIT
          ENDDO
          IF(ig2 .GT. ig) ig2 = ig
          GroupInfo%InScatRange(1, ig) = min(GroupInfo%InScatRange(1, ig), ig2)
          !Lower Bound
          DO ig2 = ng, ig, -1
            IF(XsMacSm(ig2, ig) .NE. ZERO) EXIT
          ENDDO
          IF(ig2 .LT. ig) ig2 = ig
          GroupInfo%InScatRange(2, ig) = MAX(GroupInfo%InScatRange(2, ig), ig2)
        ENDDO

        DO ig = 1, ng
          !Upper Bound
          DO ig2 = 1, ig
            IF(XsMacSm(ig, ig2) .NE. ZERO) EXIT
          ENDDO
          IF(ig2 .GT. ig) ig2 = ig
          GroupInfo%OutScatRange(1, ig) = min(GroupInfo%OutScatRange(1, ig), ig2)
          !Lower Bound
          DO ig2 = ng, ig, -1
            IF(XsMacSm(ig, ig2) .NE. ZERO) EXIT
          ENDDO
          IF(ig2 .LT. ig) ig2 = ig
          GroupInfo%OutScatRange(2, ig) = MAX(GroupInfo%OutScatRange(2, ig), ig2)
        ENDDO
      END DO
    ELSE
      CALL XsSmBen(i, 1, ng, 1, ng, XsMacSm)
      !InScattering
      DO ig = 1, ng
        !Upper Bound
        DO ig2 = 1, ig
          IF(XsMacSm(ig2, ig) .NE. ZERO) EXIT
        ENDDO
        IF(ig2 .GT. ig) ig2 = ig
        GroupInfo%InScatRange(1, ig) = min(GroupInfo%InScatRange(1, ig), ig2)
        !Lower Bound
        DO ig2 = ng, ig, -1
          IF(XsMacSm(ig2, ig) .NE. ZERO) EXIT
        ENDDO
        IF(ig2 .LT. ig) ig2 = ig
        GroupInfo%InScatRange(2, ig) = MAX(GroupInfo%InScatRange(2, ig), ig2)
      ENDDO

      DO ig = 1, ng
        !Upper Bound
        DO ig2 = 1, ig
          IF(XsMacSm(ig, ig2) .NE. ZERO) EXIT
        ENDDO
        IF(ig2 .GT. ig) ig2 = ig
        GroupInfo%OutScatRange(1, ig) = min(GroupInfo%OutScatRange(1, ig), ig2)
        !Lower Bound
        DO ig2 = ng, ig, -1
          IF(XsMacSm(ig, ig2) .NE. ZERO) EXIT
        ENDDO
        IF(ig2 .LT. ig) ig2 = ig
        GroupInfo%OutScatRange(2, ig) = MAX(GroupInfo%OutScatRange(2, ig), ig2)
      ENDDO
    END IF
  ENDDO  !End of Base MacroScopic XS sweep
ENDIF
GroupInfo%UpScatRange(2) = ng
GroupInfo%lUpScat = FALSE
DO ig = 1, ng
  IF(GroupInfo%InScatRange(2, ig) .GT. ig) THEN
    GroupInfo%lUpScat = TRUE
    GroupInfo%UpScatRange(1) = ig
    EXIT
  ENDIF
ENDDO
IF(ng .EQ. 47) THEN
  GroupInfo%UpScatRange = (/24, 47/)
ENDIF
END SUBROUTINE

SUBROUTINE DcplPrepFxr()
USE PARAM
USE TYPEDEF,         ONLY : coreinfo_type      ,Pin_Type       ,Cell_Type             &
                           ,Mixture_Type       ,FxrInfo_Type   ,DcplFxrInfo_Type      &
                           ,PE_Type
USE GEOM,            ONLY : Core               ,ng
USE Core_Mod,        ONLY : GroupInfo
USE DcplCore_Mod,    ONLY : DcplInfo           ,DcplFxr            
USE Material_mod,    ONLY : Mixture
USE CNTL,            ONLY : nTracerCntl
USE xslib_mod,       ONLY : CoreResIsoUpdate, ldiso, nelthel, mapnucl
USE BasicOperation,  ONLY : CP_VA
USE PE_MOD,          ONLY : PE
USE ALLOCS
IMPLICIT NONE

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(Mixture_Type), POINTER :: Mix       

INTEGER :: irefpln
INTEGER :: nRefPln, nRefTemp
INTEGER :: nCoreFxr, nCoreFsr, nxy, FsrIdxSt, FxrIdxSt, nFsrInFxr
INTEGER :: i, j, k, ii
INTEGER :: iz, ipin, icel, ireg, ifsr, ifxr, ixsreg, ifsrbeg, ig, imix, imix1, imixp
INTEGER :: niso, ntiso, nchi, id

REAL :: AREA
LOGICAL :: lXsLib, lres, lflag, lcladplane

lxsLib = nTracerCntl%lXsLib
IF(lXsLib) lres = nTracerCntl%lrestrmt
nRefPln = DcplInfo%nRefPln; nRefTemp = DcplInfo%nRefTemp
ALLOCATE(DcplFxr(nRefPln))
ng = GroupInfo%ng
!PRINT *, PE%myrank, PE%lmyRefPln(1:4)
DO iRefPln = 1, nRefPln
  IF(.NOT. PE%lmyRefPln(iRefPln)) CYCLE
  iz = DcplInfo%RefPln(iRefPln)
  Pin => Core%Pin
  CellInfo => Core%CellInfo
  nCoreFxr = Core%nCoreFxr; nCoreFsr = Core%nCoreFsr
  nxy = Core%nxy 
  ALLOCATE(DcplFxr(iRefPln)%Fxr(nCoreFxr, iz:iz))
  Fxr => DcplFxr(iRefPln)%Fxr
  
  DO ipin = 1, nxy
    FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
    icel = Pin(ipin)%Cell(iz)
    lflag=.true.
    DO j = 1, CellInfo(icel)%nFxr
      ixsreg = FxrIdxSt + j - 1
      nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
      ifsrbeg = CellInfo(icel)%MapFxr2FsrIdx(1, j)
      area = 0
      DO i = 1, nFsrInFxr
        area = area + CellInfo(icel)%vol(i)
      ENDDO
      
      !ifsrend = CellInfo(icel)%MapFxr2FsrIdx(nFsrInFxr, j)
      !area = SUM(CellInfo(icel)%vol(ifsrbeg:ifsrend))
      
      imix = CellInfo(icel)%ireg(ifsrbeg)
      if (j.eq.1) imix1=imix
      Fxr(ixsreg, iz)%imix = imix !celtyp(it)%ireg(ir1)
      Fxr(ixsreg, iz)%area = area
      Fxr(ixsreg, iz)%nFsrInFxr = nFsrInFxr
      Fxr(ixsreg, iz)%FsrIdxSt = FsrIdxSt + ifsrbeg - 1
      Fxr(ixsreg, iz)%lFuel = CellInfo(icel)%lFuel
      IF(.NOT. lXsLib) CYCLE
      !Only for XSlibrary 
      MIX => Mixture(imix)
      Fxr(ixsreg, iz)%lfuel = Mix%lfuel; Fxr(ixsreg, iz)%ldepl = Mix%ldepl
      Fxr(ixsreg, iz)%lCLD = Mix%lCLD
      Fxr(ixsreg, iz)%lAIC = Mix%lAIC
      IF(lres .and. Mixture(imix)%lres) THEN
          IF (CellInfo(icel)%geom%lCircle) Fxr(ixsreg, iz)%lres = TRUE
      ENDIF
      Fxr(ixsreg, iz)%lh2o = Mix%lh2o;  Fxr(ixsreg, iz)%temp = Mix%temp
        !Allocation Isotope list and thems compostion list
      niso = Mix%niso; ntiso = GroupInfo%ntiso; nchi = GroupInfo%nchi
      IF(Fxr(ixsreg, iz)%ldepl) THEN
        CALL Dmalloc(Fxr(ixsreg, iz)%idiso, ntiso)
        CALL Dmalloc(Fxr(ixsreg, iz)%pnum, ntiso)
        CALL Dmalloc(Fxr(ixsreg, iz)%chi, nchi)
      ELSEIF(Fxr(ixsreg, iz)%lh2O .AND. nTracerCntl%lInitBoron) THEN
        CALL Dmalloc(Fxr(ixsreg, iz)%idiso, niso + 1)
        CALL Dmalloc(Fxr(ixsreg, iz)%pnum, niso + 1)        
      ELSE
        CALL Dmalloc(Fxr(ixsreg, iz)%idiso, niso + 1)
        CALL Dmalloc(Fxr(ixsreg, iz)%pnum, niso + 1)
      ENDIF
      Fxr(ixsreg, iz)%niso = niso
      CALL CP_VA(Fxr(ixsreg,iz)%idiso(1:niso),mix%idiso(1:niso),niso)
      CALL CP_VA(Fxr(ixsreg,iz)%pnum(1:niso),mix%pnum(1:niso),niso)
      IF (imix.eq.imix1.and.lflag) then
          Fxr(ixsreg, iz)%lRes = .FALSE.
      else
          lflag = .false.
      endif
      !Updae the lists of the resonance isotopes in the Active core
      IF(Fxr(ixsreg, iz)%lRes) CALL CoreResIsoUpdate(Fxr(ixsreg,iz)%idiso(1:niso), niso, iz)

      CALL Dmalloc0(Fxr(ixsreg, iz)%Dcpl_Temp, 1, nRefTemp)
      CALL Dmalloc0(Fxr(ixsreg, iz)%Dcpl_pnum, 1,niso+1, 1, nRefTemp)
      NULLIFY(MIX) 
    ENDDO    
    ! Define Cladding region
    IF ((.not.lres).OR.(.NOT. lXsLib)) CYCLE
    lflag=.false.
    DO j = CellInfo(icel)%nFxr, 1, -1
      ixsreg = FxrIdxSt + j - 1
      ifsrbeg = CellInfo(icel)%MapFxr2FsrIdx(1, j)

      imix = CellInfo(icel)%ireg(ifsrbeg)
      if (j.eq.CellInfo(icel)%nFxr) then
          imixp = imix
          cycle
      endif
      if (imix.ne.imixp) then
        lflag=.true.
      endif
      imixp = imix
      MIX => Mixture(imix)
      IF (.not.Fxr(ixsreg, iz)%lRes) cycle      
      IF (lflag) then
        do ii=1,Fxr(ixsreg,iz)%niso
          id=Fxr(ixsreg,iz)%idiso(ii)
          if (ldiso(mapnucl(id))%lclad) exit
        enddo
        if (ii.le.Fxr(ixsreg,iz)%niso) then
            Fxr(ixsreg, iz)%lCLD=.true.
            lcladplane=.TRUE.
        endif
      ENDIF
    ENDDO
  ENDDO
  
  
ENDDO
CONTINUE
NULLIFY(CellInfo, Pin, Fxr)
NULLIFY(Mix)
ENDSUBROUTINE

SUBROUTINE AllocDcplXsEQ()
!PrePare
USE PARAM
USE TYPEDEF,      ONLY : coreinfo_type, FxrInfo_Type
USE GEOM,         ONLY : Core
USE CORE_mod,     ONLY : GroupInfo
USE DcplCore_mod, ONLY : DcplInfo, DcplFxr
USE CNTL,         ONLY : nTracerCntl    
USE PE_MOD,       ONLY : PE
USE XsLib_Mod,    ONLY : nlvflxmax, nreshel
USE ALLOCS
USE BasicOperation, ONLY : CP_CA
IMPLICIT NONE
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)

INTEGER :: nFxr, nSubFlx, ig, iResGrpBeg, iResGrpEnd, nCatMg
INTEGER :: nRefPln, nRefTemp, myRefPlnBeg, myRefPlnEnd
INTEGER :: myzb, myze
INTEGER :: iz, ifxr, iRefPln

nFxr = Core%nCoreFXR;
nRefPln = DcplInfo%nRefPln; nRefTemp = DcplInfo%nRefTemp
iResGrpBeg = GroupInfo%nofg + 1;
iResGrpEnd = GroupInfo%nofg + GroupInfo%norg

myzb = PE%myzb; myze = PE%myze

#ifndef MPI_ENV
myRefPlnBeg = 1; myRefPlnEnd = DcplInfo%nRefPln
#else
myRefPlnBeg = PE%myRefPlnBeg; myRefPlnEnd = PE%myRefPlnEnd
#endif
DO iRefPln = myRefPlnBeg, myRefPlnEnd
  myzb = DcplInfo%RefPln(iRefPln); myze = myzb
  Fxr => DcplFxr(iRefPln)%Fxr
  DO iz = myzb, myze
    DO ifxr = 1, nFxr
      
      IF (.NOT. Fxr(ifxr,iz)%lres) CYCLE
      !IF (Fxr(ifxr,iz)%lCLD) THEN
      !  CALL Dmalloc(Fxr(ifxr, iz)%XsEq_c, mlgdata%c_nmaclv, mlgdata%c_nmaclp, mlgdata%c_nmaclpf)
      !ELSE
      !  CALL Dmalloc(Fxr(ifxr, iz)%XsEq_f, mlgdata%f_nmaclv, mlgdata%f_nmaclp)
      !ENDIF 
      CALL Dmalloc0(Fxr(ifxr, iz)%XsEq, 1, nlvflxmax, 1, nreshel, iResGrpBeg, iResGrpEnd)   
      Fxr(ifxr, iz)%XsEq = 0._8
      CALL Dmalloc0(Fxr(ifxr, iz)%fresoa, iResGrpBeg, iResGrpEnd)
      CALL Dmalloc0(Fxr(ifxr, iz)%fresoS, iResGrpBeg, iResGrpEnd)
      CALL Dmalloc0(Fxr(ifxr, iz)%fresoStr, iResGrpBeg, iResGrpEnd)
      CALL Dmalloc0(Fxr(ifxr, iz)%fresoF, iResGrpBeg, iResGrpEnd)
      CALL Dmalloc0(Fxr(ifxr, iz)%fresoNF, iResGrpBeg, iResGrpEnd)  ! nu-fission resonance treatment  EDIT JSU 20190810
      CALL Dmalloc0(Fxr(ifxr, iz)%fresokF, iResGrpBeg, iResGrpEnd)  ! kappa-fission resonance treatment  EDIT JSU 20190816
      !CALL CP_CA(Fxr(ifxr, iz)%fresoA(iResGrpBeg:iResGrpEnd), 1._8, iResGrpEnd-iResGrpBeg+1)
      !CALL CP_CA(Fxr(ifxr, iz)%fresoS(iResGrpBeg:iResGrpEnd), 1._8, iResGrpEnd-iResGrpBeg+1)
      !CALL CP_CA(Fxr(ifxr, iz)%fresoStr(iResGrpBeg:iResGrpEnd), 1._8, iResGrpEnd-iResGrpBeg+1)
      !CALL CP_CA(Fxr(ifxr, iz)%fresoF(iResGrpBeg:iResGrpEnd), 1._8, iResGrpEnd-iResGrpBeg+1)
      CALL Dmalloc0(Fxr(ifxr, iz)%Dcpl_fresoa, iResGrpBeg, iResGrpEnd, 0, nRefTemp)
      CALL Dmalloc0(Fxr(ifxr, iz)%Dcpl_fresoS, iResGrpBeg, iResGrpEnd, 0, nRefTemp)
      CALL Dmalloc0(Fxr(ifxr, iz)%Dcpl_fresoStr, iResGrpBeg, iResGrpEnd, 0, nRefTemp)
      CALL Dmalloc0(Fxr(ifxr, iz)%Dcpl_fresoF, iResGrpBeg, iResGrpEnd, 0, nRefTemp)
      Fxr(ifxr, iz)%Dcpl_fresoA = 1._8
      Fxr(ifxr, iz)%Dcpl_fresoS = 1._8
      Fxr(ifxr, iz)%Dcpl_fresoStr = 1._8
      Fxr(ifxr, iz)%Dcpl_fresoF = 1._8
      
      ALLOCATE(Fxr(ifxr, iz)%fresoAIso(Fxr(ifxr, iz)%ndim,iResGrpBeg:iResGrpEnd)) 
      ALLOCATE(Fxr(ifxr, iz)%fresoFIso(Fxr(ifxr, iz)%ndim,iResGrpBeg:iResGrpEnd)) 
      Fxr(ifxr, iz)%fresoAIso=1._8
      Fxr(ifxr, iz)%fresoFIso=1._8
      IF (nTRACERCNTL%lGamma) THEN
        ALLOCATE(Fxr(ifxr, iz)%fresocapIso(Fxr(ifxr, iz)%ndim,iResGrpBeg:iResGrpEnd)) 
        Fxr(ifxr, iz)%fresocapIso=1._8
      END IF
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE GetFxrHmMass(Fxr, hz)
USE PARAM
USE TYPEDEF,   ONLY : FxrInfo_Type
USE nuclidmap_mod, ONLY : AtomicWeight
TYPE(FxrInfo_Type) :: Fxr
REAL :: hz
REAL :: HmReg
INTEGER :: i, id

HmReg = 0
DO i = 1, Fxr%niso
  id = Fxr%idiso(i)
  IF(id .LT. 90000)  CYCLE
  HmReg = HmReg + Fxr%pnum(i) * AtomicWeight(id) / AVOGADRO
ENDDO
HmReg = HmReg * Fxr%Area * Hz / 1000._8
Fxr%HmKg0 = HmReg
END SUBROUTINE