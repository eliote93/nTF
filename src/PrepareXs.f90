#include <defines.h>
SUBROUTINE PrepFxr(lXsLib, lfxr)
!PrePare
USE PARAM
USE TYPEDEF,         ONLY : coreinfo_type,  Pin_Type,    Cell_Type,             &
                            Mixture_Type,   Fxrinfo_type
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
USE XSLIB_MOD,       ONLY : ldiso, mapnucl, CoreResIsoUpdate
IMPLICIT NONE

TYPE(Mixture_Type), POINTER :: Mix       
TYPE(Fxrinfo_type), POINTER :: FXR_Loc
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
      
      FXR_Loc => Fxr(ixsreg, iz)
      
      FXR_Loc%area = area
      FXR_Loc%ipin = ipin
      
      imix = CellInfo(icel)%ireg(ifsrbeg)
      if(j.eq.1) imix1=imix
      IF( lfxr )THEN
          FXR_Loc%imix = ixsreg !celtyp(it)%ireg(ir1)
          FXR_Loc%lFuel = MacXsBen(ixsreg)%lfuel
          IF(MacXsBen(ixsreg)%lfuel)THEN
              lFuelPin=.TRUE.      
              lFuelPlane=.TRUE.
              IF(.NOT. Core%AsyInfo(iasytype)%lfuel) Core%AsyInfo(iasytype)%lfuel=.TRUE.
          ENDIF          
      ELSE
          FXR_Loc%imix = imix !celtyp(it)%ireg(ir1)
          FXR_Loc%lFuel = CellInfo(icel)%lFuel
      ENDIF
      FXR_Loc%nFsrInFxr = nFsrInFxr
      FXR_Loc%FsrIdxSt = FsrIdxSt + ifsrbeg - 1
      IF(.NOT. lXsLib) CYCLE
      !Only for XSlibrary 
      MIX => Mixture(imix)
      FXR_Loc%lfuel = Mix%lfuel; FXR_Loc%ldepl = Mix%ldepl
      FXR_Loc%lCLD = Mix%lCLD
      FXR_Loc%lAIC = Mix%lAIC
      IF(lres .and. Mixture(imix)%lres) then
          IF (CellInfo(icel)%Geom%lcircle) FXR_Loc%lres = TRUE  !! Only for circular geometry
      ENDIF
      FXR_Loc%lh2o = Mix%lh2o; FXR_Loc%temp = Mix%temp
        !Allocation Isotope list and thems compostion list
      niso = Mix%niso; nchi = GroupInfo%nchi; ntiso = GroupInfo%ntiso
      IF(FXR_Loc%ldepl) THEN
        CALL Dmalloc(FXR_Loc%idiso, ntiso)
        CALL Dmalloc(FXR_Loc%pnum, ntiso)
        CALL Dmalloc(FXR_Loc%chi, nchi)
        FXR_Loc%ndim = ntiso
      ELSEIF(FXR_Loc%lh2o .AND. nTracerCntl%lInitBoron) THEN
        CALL Dmalloc(FXR_Loc%idiso, niso + 4)
        CALL Dmalloc(FXR_Loc%pnum, niso + 4)
        FXR_Loc%ndim = niso + 4
      ELSEIF(FXR_Loc%lh2o) THEN
        CALL Dmalloc(FXR_Loc%idiso, niso + 4)
        CALL Dmalloc(FXR_Loc%pnum, niso + 4)
        FXR_Loc%ndim = niso + 4
      ELSE
        CALL Dmalloc(FXR_Loc%idiso, niso)
        CALL Dmalloc(FXR_Loc%pnum, niso)
        FXR_Loc%ndim = niso
      ENDIF
      FXR_Loc%niso = niso; FXR_Loc%niso_depl = niso
      CALL CP_VA(Fxr(ixsreg,iz)%idiso(1:niso),mix%idiso(1:niso),niso)
      CALL CP_VA(Fxr(ixsreg,iz)%pnum(1:niso),mix%pnum(1:niso),niso)

      CALL FIndGdFxr(FXR_Loc)
      IF(FXR_Loc%lDepl) CALL GetFxrHmMass(Fxr(ixsreg,iz), Core%hz(iz))
      if (imix.eq.imix1.and.lflag) then
          FXR_Loc%lRes = .FALSE.
      else
          lflag=.false.
      endif
      IF (FXR_Loc%lRes .AND. FXR_Loc%lCLD) lcladplane=.TRUE.
#ifndef MPI_ENV
      !Updae the lists of the resonance isotopes in the Active core
      IF(FXR_Loc%lRes) CALL CoreResIsoUpdate(Fxr(ixsreg,iz)%idiso(1:niso), niso, iz)
#endif
      NULLIFY(MIX)
    ENDDO    
    !Pin(ipin)%lFuel=lFuelPin
    IF( lFXR )    Pin(ipin)%lFuel=lFuelPin
    ! Define Cladding region
    IF ((.not.lres).OR.(.NOT. lXsLib)) CYCLE
    lflag=.false.
    DO j = CellInfo(icel)%nFxr, 1, -1
      ixsreg = FxrIdxSt + j - 1
      ifsrbeg = CellInfo(icel)%MapFxr2FsrIdx(1, j)
      
      FXR_Loc => Fxr(ixsreg, iz)
      
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
      IF (.not.FXR_Loc%lRes) cycle
      IF (lflag) then
        do ii=1,Fxr(ixsreg,iz)%niso
          id=Fxr(ixsreg,iz)%idiso(ii)
          if (ldiso(mapnucl(id))%lclad) exit
        enddo
        if (ii.le.Fxr(ixsreg,iz)%niso) then
            FXR_Loc%lCLD=.true.
            lcladplane=.TRUE.
        else
            FXR_Loc%lCLD=.false.
            if (.not.FXR_Loc%lfuel) FXR_Loc%lres=.false.
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
!Updae the lists of the resonance isotopes in the Active core
IF(lXsLib) CALL MakeResIsoList()
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
!CALL  Calnum(Core, FmInfo, PE)
!CALL  ReadFxr(Core, FmInfo, PE)
!CALL SetStructure(Core, Fxr, 0.0859375, myzb, myze)
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
    
SUBROUTINE AllocResIsoInfo()
USE PE_MOD,        ONLY : PE
USE CNTL,          ONLY : nTracerCntl    
USE XSLIB_MOD,     ONLY : nreshel, maxnid, nActiveCat, IDRES_COR, NRES_COR, mapnuclRes, mlgdata
USE ALLOCS
USE nuclidmap_mod,  ONLY : InitResoCat
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

END SUBROUTINE
    
SUBROUTINE AllocXsEQ()
!PrePare
USE PARAM
USE TYPEDEF,       ONLY : Cell_TYpe, Pin_Type, ResVarPin_Type, Fxrinfo_type
USE GEOM,          ONLY : Core
USE CORE_mod,      ONLY : Fxr,            GroupInfo
USE CNTL,          ONLY : nTracerCntl    
USE PE_MOD,        ONLY : PE
USE BasicOperation,ONLY : CP_CA
USE XSLIB_MOD,     ONLY : nlvflxmax, nreshel, nelthel, mlgdata0, nCat
USE ALLOCS
USE CP_mod,        ONLY : ngauss
IMPLICIT NONE
TYPE(Cell_TYpe), POINTER :: CellInfo(:)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(ResVarPin_Type), POINTER :: ResVarPin(:,:)
TYPE(ResVarPin_Type), POINTER :: rvPin_Loc
TYPE(Fxrinfo_type), POINTER :: FXR_Loc
INTEGER :: nFxr, iResGrpBeg, iResGrpEnd
INTEGER :: myzb, myze
INTEGER :: iz, ifxr, icel, ipin, nxy, FxrIdxSt, j, ig, igt, ngt, nlocalFxr
REAL :: fuelrad(100)
  
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
            ResVarPin(ipin,iz)%lresA=.TRUE.
            IF (.not.Fxr(ifxr, iz)%lCLD) then
                ResVarPin(ipin,iz)%lres=.TRUE.
            ELSE
                ResVarPin(ipin,iz)%lresC=.TRUE.
            ENDIF
        ENDIF
    ENDDO
  ENDDO
ENDDO

IF (nTracerCntl%lED) THEN
    IF (nTracerCntl%lMLG) CALL SetMLGrid_ED
ELSE
    IF (nTracerCntl%lMLG) CALL SetMLGrid
ENDIF

!$OMP PARALLEL PRIVATE(iz, ifxr, FXR_Loc, ig)
!$OMP DO SCHEDULE(GUIDED)
DO iz = myzb, myze    
  if (.not.Core%lFuelPlane(iz)) CYCLE
  
  DO ifxr = 1, nFxr
    FXR_Loc => Fxr(ifxr, iz)
    
    IF (.NOT. FXR_Loc%lres) CYCLE
    
    CALL Dmalloc0(FXR_Loc%fresoa, iResGrpBeg, iResGrpEnd)
    CALL Dmalloc0(FXR_Loc%fresof, iResGrpBeg, iResGrpEnd)    
    do ig = iResGrpBeg,iResGrpEnd
      FXR_Loc%fresoa(ig)=1.
      FXR_Loc%fresof(ig)=1.
    enddo
    
    ALLOCATE(FXR_Loc%fresoAIso(FXR_Loc%ndim,iResGrpBeg:iResGrpEnd)) 
    ALLOCATE(FXR_Loc%fresoFIso(FXR_Loc%ndim,iResGrpBeg:iResGrpEnd)) 
    do ig = iResGrpBeg,iResGrpEnd
      FXR_Loc%fresoAIso(:,ig)=1.
      FXR_Loc%fresoFIso(:,ig)=1.
    ENDDO

    IF (nTracerCntl%lRST) THEN
      CALL Dmalloc0(FXR_Loc%fresos, iResGrpBeg, iResGrpEnd)
      CALL Dmalloc0(FXR_Loc%fresostr, iResGrpBeg, iResGrpEnd)
      do ig = iResGrpBeg,iResGrpEnd
        FXR_Loc%fresos(ig)=1.
        FXR_Loc%fresostr(ig)=1.
      enddo
    
      ALLOCATE(FXR_Loc%fresoSIso(FXR_Loc%ndim,iResGrpBeg:iResGrpEnd)) 
      ALLOCATE(FXR_Loc%fresoSSIso(FXR_Loc%ndim,iResGrpBeg:iResGrpEnd)) 
      ALLOCATE(FXR_Loc%fresoS1Iso(FXR_Loc%ndim,iResGrpBeg:iResGrpEnd)) 
      do ig = iResGrpBeg,iResGrpEnd
        FXR_Loc%fresoSIso(:,ig)=1.
        FXR_Loc%fresoSSIso(:,ig)=1.
        FXR_Loc%fresoS1Iso(:,ig)=1.
      ENDDO       
    ENDIF
    IF (nTracerCntl%lMLG) THEN
      IF (FXR_Loc%lCLD) THEN
        CALL Dmalloc(FXR_Loc%XsEq_c_1g, mlgdata0%c_nmaclv1G)
        FXR_Loc%XsEq_c_1g=0._8
      ELSEIF (FXR_Loc%lAIC) THEN
        CALL Dmalloc(FXR_Loc%XsEq_f_1g, mlgdata0%f_nmaclv1G)
        FXR_Loc%XsEq_f_1g=0._8  
      ELSE
        CALL Dmalloc0(FXR_Loc%XsEq_f_mg, 1, mlgdata0%f_nmaclv, iResGrpBeg, iResGrpEnd)   
        FXR_Loc%XsEq_f_mg=0._8
        CALL Dmalloc0(FXR_Loc%FnAdj, iResGrpBeg, iResGrpEnd)  
        FXR_Loc%FnAdj=1._8
        CALL Dmalloc0(FXR_Loc%FtAdj, 1, mlgdata0%f_nmaclv, iResGrpBeg, iResGrpEnd)  
        FXR_Loc%FtAdj=1._8  
      ENDIF
    ELSE
      IF (nTracerCntl%lCAT) THEN
        CALL Dmalloc0(FXR_Loc%XsEq, 1, nlvflxmax, 1, nCat, iResGrpBeg, iResGrpEnd)   
        CALL Dmalloc0(FXR_Loc%NDAF, 1, nlvflxmax, 1, nCat, iResGrpBeg, iResGrpEnd)    
      ELSE
        CALL Dmalloc0(FXR_Loc%XsEq, 1, nlvflxmax, 1, nreshel, iResGrpBeg, iResGrpEnd)   
        CALL Dmalloc0(FXR_Loc%NDAF, 1, nlvflxmax, 1, nreshel, iResGrpBeg, iResGrpEnd)    
      ENDIF
      FXR_Loc%XsEq = 0._8
      FXR_Loc%NDAF = 1._8
    ENDIF
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL PRIVATE(iz, ipin, rvPin_Loc, icel, nlocalfxr)
!$OMP DO SCHEDULE(GUIDED)
DO iz = myzb, myze
  if (.not.Core%lFuelPlane(iz)) CYCLE
  DO ipin = 1, nxy
    rvPin_Loc => ResVarPin(ipin,iz)
    
    if (.not.rvPin_Loc%lres) CYCLE  
    icel = Pin(ipin)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
    ALLOCATE(rvPin_Loc%rad_cp(nlocalFxr))
    ALLOCATE(rvPin_Loc%delr(nlocalFxr))
    ALLOCATE(rvPin_Loc%vol(nlocalFxr))
    ALLOCATE(rvPin_Loc%Qsurfvol(nlocalFxr))
    ALLOCATE(rvPin_Loc%X(nlocalFxr,nlocalFxr,ngauss))
    
    rvPin_Loc%rad_cp(1:nlocalFxr)=CellInfo(icel)%rad_cp(1:nlocalFxr)    
    CALL Dmalloc0(rvPin_Loc%idiso, 1, nelthel)  
    CALL Dmalloc0(rvPin_Loc%pnum, 1, nelthel)  
    
    rvPin_Loc%lbar = 2._8 * CellInfo(icel)%FuelRad0
    
    IF (.not.nTracerCntl%lMLG) THEN
      IF (nTracerCntl%lCAT) THEN
        CALL Dmalloc0(rvPin_Loc%avgxseq, 1, nlvflxmax, 1, nCat, iResGrpBeg, iResGrpEnd)  
      ELSE
        CALL Dmalloc0(rvPin_Loc%avgxseq, 1, nlvflxmax, 1, nreshel, iResGrpBeg, iResGrpEnd)
      ENDIF
      rvPin_Loc%avgxseq=0._8
    ELSE    
      CALL Dmalloc0(rvPin_Loc%avgxseq_mg, 1, mlgdata0%f_nmaclv, iResGrpBeg, iResGrpEnd)  
      rvPin_Loc%avgxseq_mg=0._8
      CALL Dmalloc0(rvPin_Loc%avgxseq_1g, 1, mlgdata0%f_nmaclv1G)  
      rvPin_Loc%avgxseq_1g=0._8
      CALL Dmalloc0(rvPin_Loc%FnAdj, iResGrpBeg, iResGrpEnd)  
      rvPin_Loc%FnAdj=1._8
    ENDIF
    IF (nTracerCntl%lRIF) THEN
      CALL Dmalloc0(rvPin_Loc%rifa, 1, nreshel, iResGrpBeg, iResGrpEnd)
      CALL Dmalloc0(rvPin_Loc%riff, 1, nreshel, iResGrpBeg, iResGrpEnd)   
      IF (nTracerCntl%lRST) CALL Dmalloc0(rvPin_Loc%rifs, 1, nreshel, iResGrpBeg, iResGrpEnd)
    ENDIF
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

NULLIFY (FXR_Loc)
NULLIFY (rvPin_Loc)

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

  nmlg=nTracerCntl%nMLG
  mlgdata0%f_nmaclv=nmlg
  mlgdata0%f_nmaclv1G=nmlg
  
  mlgdata0%c_nmaclv1G=5
  allocate(mlgdata0%c_maclv1G(mlgdata0%c_nmaclv1G))
  allocate(mlgdata0%c_maclv1G_log(mlgdata0%c_nmaclv1G))
  mlgdata0%c_maclv1G=(/0.1_8,1._8,10._8,100._8,1000._8/)  
  mlgdata0%c_maclv1G_log=dlog(mlgdata0%c_maclv1G)
  
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
    
    if (.not.core%lAICplane(iz)) cycle
    
    mlgdata(iz)%f_nmaclv1G=nmlg
    allocate(mlgdata(iz)%f_maclv1G(nmlg))
    allocate(mlgdata(iz)%f_maclv1G_log(nmlg))
    radmin=100._8; radmax=0._8
    DO ipin = 1, nxy  
        icel = Pin(ipin)%Cell(iz)      
        if (.not.CellInfo(icel)%lAIC) cycle
        rad=CellInfo(icel)%FuelRad0
        if (rad.le.0.2_8) cycle
        radmin=min(radmin,rad)
        radmax=max(radmax,rad)
    ENDDO
    sigmax=82._8/radmin
    sigmin=0.164_8/radmax
    del=(sigmax/sigmin)**(1._8/nmlg)
    mlgdata(iz)%f_maclv1G(1)=sigmin; mlgdata(iz)%f_maclv1G(nmlg)=sigmax
    mlgdata(iz)%f_maclv1G_log(1)=dlog(sigmin); mlgdata(iz)%f_maclv1G_log(nmlg)=dlog(sigmax)
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
  
  mlgdata0%c_nmaclv1G=5
  allocate(mlgdata0%c_maclv1G(mlgdata0%c_nmaclv1G))
  allocate(mlgdata0%c_maclv1G_log(mlgdata0%c_nmaclv1G))
  mlgdata0%c_maclv1G=(/0.1_8,1._8,10._8,100._8,1000._8/)  
  mlgdata0%c_maclv1G_log=dlog(mlgdata0%c_maclv1G)
  
  do iz=myzb,myze
    if (.not.core%lfuelplane(iz)) cycle
    
    mlgdata(iz)%f_nmaclv=nmlg
    mlgdata(iz)%f_nmaclv1G=nmlg
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
    ENDDO
  enddo
  
END SUBROUTINE

!Set Group Sweeping and Scattering Range Setting
SUBROUTINE SetGroupInfo(GroupInfo,lXsLib)
USE PARAM
USE TYPEDEF,           ONLY : GroupInfo_Type
!USE CORE_MOD,          ONLY : GroupInfo
USE GEOM,              ONLY : ng,              nPrec
USE BenchXs,           ONLY : nXslType,        XsSmBen
USE ALLOCS             
USE XsLib_Mod,         ONLY : nelthel,         noghel,        nelrhel,      nchihel,              &
                              ldiso       
USE TranMacXsLib_Mod,  ONLY : SetTranXsLibInfo
USE cntl,              ONLY : nTracerCntl
IMPLICIT NONE

TYPE(GroupInfo_Type) :: GroupInfo
LOGICAL :: lXsLib
REAL, POINTER :: XsMacSm(:, :)
INTEGER :: i, j, k, ig, ig2, iso
INTEGER :: nid, ntemp, itemp
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
ELSE
  DO i = 1, nXslType
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