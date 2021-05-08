#include <defines.h>
#include <DefDBG.h>
#include <Depletion.h>
#ifdef __PGI

MODULE CUDA_DEPL
  USE HPDeplMod
  USE HPDeplType
  USE DEPL_MOD, ONLY : DeplCntl, decayXe135, decayI135, yieldXe135, yieldI135
  USE CUDA_MASTER, ONLY : cuDepl
  IMPLICIT NONE
  TYPE(DeplFxrBundle_Type) :: DeplFxrBundle
  TYPE(DeplSysBundle_Type) :: DeplSysBundle
  TYPE(DeplLib_Type) :: DeplLib
  INTEGER :: DeplSolTyp=1, nDeplTh=1, PC_OPT=1
  INTEGER :: nSubStep, SysByte, Scale

  INTEGER, POINTER :: MapXs2Dep(:), MapDep2XS(:)
  INTEGER, POINTER :: DeplFxrMap(:,:)
  INTEGER, PARAMETER :: nIsoGd = 7
  INTEGER :: IdIsoGd(7) = (/641520, 641540, 641550, 641560, 641570, 641580, 641600/)

  ! Conditioners
  LOGICAL :: lCallInitDepl = .FALSE., lInitMapIso = .FALSE.
  CONTAINS

  SUBROUTINE Init_Depl(PE)
  USE files, ONLY : Filename, DeplFileidx, io_quick
  USE TYPEDEF, ONLY : PE_Type
  IMPLICIT NONE
  TYPE(PE_Type) :: PE

  IF (lCallInitDepl) RETURN
  DeplSolTyp = DeplCntl%SOLVER; nDeplTh = PE%nDeplThread; PC_OPT = DeplCntl%PC_OPT
  nSubStep = cuDepl%nSubStep; SysByte = cuDepl%SysByte; Scale = cuDepl%Scale
#ifdef GDQUAD_SIMPLE
  nSubStep = nSubStep*2
#endif
  CALL DeplLibInit(io_quick, FILENAME(deplFileIdx), DeplLib)
  lCallInitDepl = .TRUE.
  END SUBROUTINE

  SUBROUTINE AllocDeplFxrMem(Core, Fxr, GroupInfo, PE)
  USE PARAM
  USE allocs
  USE TYPEDEF,    ONLY : CoreInfo_Type,     FxrInfo_Type,         PE_TYPE,          &
                         GroupInfo_Type,    Pin_Type,             Cell_Type
  USE nuclidmap_mod, ONLY : iposiso
  IMPLICIT NONE
  TYPE(CoreInfo_Type) :: Core
  TYPE(FxrInfo_Type), POINTER :: FXR(:, :)
  TYPE(GroupInfo_Type) :: GroupInfo
  TYPE(PE_Type) :: PE

  TYPE(Pin_Type), POINTER :: Pin(:)
  TYPE(Cell_Type), POINTER :: Cell(:)
  TYPE(FxrInfo_Type), POINTER :: myFxr

  INTEGER :: ntiso, nFxr, nxy, myzb, myze
  INTEGER :: FxrIdxSt, nLocalFxr
  INTEGER :: ixy, iz, ifxr, icell, i, j, thisoid, ntisoid
  INTEGER :: nIsoDepl   ! 16/02/11 Depletion timestep bug fixed
  INTEGER :: ndfxr, ndnet, ndgd      ! 19/01/22 Edit by LHG

  TYPE(DeplFxr_Type), POINTER :: mydFxr
  TYPE(GdFxrInfo_Type), POINTER :: mygdFxr

  myzb = PE%myzb; myze = PE%myze
  nFxr = Core%nCoreFxr; nxy = Core%nxy
  ntiso = GroupInfo%ntiso
  nIsoDepl = DeplLib%nofiso
  Pin => Core%Pin; Cell => Core%CellInfo

  DeplFxrBundle%nIsoGd = nIsoGd;
  ALLOCATE(DeplFxrBundle%IdIsoGd(nIsoGd))
  DO i = 1, nIsoGd
    DeplFxrBundle%IdIsoGd(i) = IsoidSrch(DeplLib%IdIso, nIsoDepl, IdIsoGd(i))
  END DO
  !DeplFxrBundle%IdIsoGd = IdIsoGd
  DeplFxrBundle%i155 = IsoidSrch(DeplLib%IdIso, nIsoDepl, 641550)
  ndgd = 0; ndnet = 0;

  !DeplFxrBundle%Hmkg0 = DeplCntl%Hm_Mass0_kg

  ! Original nTRACER Data -- Edit by LHG on 19/01/22
  DO iz = myzb, myze
    DO ixy = 1, nxy
      icell = Pin(ixy)%Cell(iz)
      FxrIdxSt = Pin(ixy)%FxrIdxSt; nLocalFxr = Cell(icell)%nFxr
      DO i = 1, nLocalFxr
        iFxr = FxrIdxSt + i - 1
        IF(.NOT. Fxr(iFxr, iz)%lDepl) CYCLE
        myFxr => Fxr(iFxr, iz)
        ndnet = ndnet+1; ! 10/01/22 edit by LHG

        CALL Dmalloc(myFxr%idiso_past, ntiso)
        CALL Dmalloc(myFxr%pnum_past, ntiso)
        !N.D. Vectors for saving all isotopes specified in Depl. Lib.
        CALL Dmalloc(myFxr%pnum_all, nisodepl)
        CALL Dmalloc(myFxr%pnum_past_all, nisodepl)
        !--- end  ! 16/02/11 Depletion timestep bug fixed
        myFxr%l_pnum_all = .FALSE.
        !myFxr%burnup = 0; myFxr%burnup_past = 0
#ifdef GdHO
        IF(.NOT. myFxr%lGD) CYCLE
#endif
        ndgd = ndgd+1; ! 19/01/22 edit by LHG
        ALLOCATE(myFxr%DeplPCP)
        ALLOCATE(myFxr%DeplXS1g(-1:0))
        DO j = -1, 0
          CALL Dmalloc(myFxr%DeplXs1g(j)%idiso, ntiso)
          CALL Dmalloc(myFxr%DeplXs1g(j)%xsa,   ntiso)
          CALL Dmalloc(myFxr%DeplXs1g(j)%xsf,   ntiso)
          CALL Dmalloc(myFxr%DeplXs1g(j)%xsn2n, ntiso)
          CALL Dmalloc(myFxr%DeplXs1g(j)%xsn3n, ntiso)
        ENDDO
      ENDDO
    ENDDO
  ENDDO


  ! New DM Data -- Edit by LHG on 19/01/22
  DeplFxrBundle%nfxr = nFxr*(myze-myzb+1)
  DeplFxrBundle%nTrueDepl = ndnet
  DeplFxrBundle%nTrueGd = ndgd
  ALLOCATE(DeplFxrBundle%FxrBundle(nFxr*(myze-myzb+1)))
  IF (ndnet.GT.0) ALLOCATE(DeplFxrBundle%IdTrueDepl(ndnet))
  IF (ndgd.GT.0) THEN
    ALLOCATE(DeplFxrBundle%GdFxrBundle(ndgd))
    ALLOCATE(DeplFxrBundle%IdTrueGd(ndgd))
  END IF
  DeplFxrBundle%delT = 0.; DeplFxrBundle%Hmkg0 = 0.; DeplFxrBundle%Power = 0.;

  ALLOCATE(DeplFxrMap(nFxr, myzb:myze))
  DeplFxrMap = 0

  ndfxr = 0; ndnet = 0; ndgd = 0;
  DO iz = myzb, myze
    DO ixy = 1, nxy
      icell = Pin(ixy)%Cell(iz)
      FxrIdxSt = Pin(ixy)%FxrIdxSt; nLocalFxr = Cell(icell)%nFxr
      DO i = 1, nLocalFxr
        ndfxr = ndfxr+1; mydFxr => DeplFxrBUndle%FxrBundle(ndfxr)
        iFxr = FxrIdxSt + i - 1
        DeplFxrMap(iFxr, iz) = ndfxr
        mydFxr%lDepl = .FALSE.; mydFxr%lFuel = .FALSE.; mydFxr%lGd = .FALSE.;
        IF(.NOT. Fxr(iFxr, iz)%lDepl) CYCLE
        ndnet = ndnet+1;
        myFxr => Fxr(iFxr, iz)
        DeplFxrBundle%IdTrueDepl(ndnet) = ndfxr
        mydFxr%iDeplFxr = ndnet
        mydFxr%NisoEig = ntiso
        mydFxr%lDepl = .TRUE.
        IF(Fxr(iFxr,iz)%lfuel) mydFxr%lFuel = .TRUE.
        ALLOCATE(mydFxr%pnum_pre(nisodepl), mydFxr%pnum_cor(nisodepl))
        mydFxr%pnum_pre = 0.; mydFxr%pnum_cor = 0.;
        ALLOCATE(mydFxr%pnum_sseig(ntiso), mydFxr%pnum_depl(nisodepl))
        mydFxr%pnum_depl = 0._8; mydFxr%pnum_sseig = 0._8;
        !mydFxr%pnum_sseig = myFxr%pnum;
        ALLOCATE(mydFxr%xs1g(DeplLib%NofRct, nisodepl), mydFxr%Kappa(nisodepl))
        ALLOCATE(mydFxr%idisoEig(ntiso))
        mydFxr%burnup = 0.; mydFxr%Tmpt = myFxr%temp; mydFxr%Vol = myFxr%area*Core%hz(iz)
        mydFxr%Hmkg0 = myFxr%Hmkg0; mydFxr%NormFlux1g = 0.; mydFxr%idisoEig = 0
        mydFxr%idisoEig = MapXs2Dep;

        mydFxr%xs1g(:,:) = DeplLib%RctXS(:,:)
        mydFxr%Kappa(:) = DeplLib%Kappa(:)

        DO j = 1, ntiso
          thisoid = myFxr%idiso(j)
          IF (thisoid .EQ. 0) CYCLE
          ntisoid = iposiso(thisoid)
          IF (ntisoid.GT.0 .AND. ntisoid.LE.ntiso) THEN
            mydFxr%pnum_sseig(ntisoid) = myFxr%pnum(j)
          END IF
          thisoid = MapXS2Dep(ntisoid)
          IF (thisoid.GT.0 .AND. thisoid.LE.nisodepl) THEN
            mydfxr%pnum_depl(thisoid) = myFxr%pnum(j)
          END IF
        END DO
#ifdef GdHO
        IF(.NOT. myFxr%lGD) CYCLE
#endif
        ndgd = ndgd+1;
        DeplFxrBundle%IdTrueGd(ndgd) = ndfxr
        mygdFxr => DeplFxrBundle%GdFxrBundle(ndgd)
        mygdFxr%aFxr => DeplFxrBundle%FxrBundle(ndfxr)
        ALLOCATE(mygdFxr%Gdpnum(2,-1:1,nIsoGd))
        ALLOCATE(mygdFxr%f_pc(nIsoGd))
        ALLOCATE(mygdFxr%GDRR(-1:1,nIsoGd))
        ALLOCATE(mygdFxr%Gd155(-1:1))
        ALLOCATE(mygdFxr%c_qd(3,nIsoGd))
        mygdFxr%f_pc = 0.; mygdFxr%Gdpnum = 0.;
        mygdFxr%c_qd = 0.; mygdFxr%GDRR = 0.; myGdFxr%Gd155 = 0.
      ENDDO
    ENDDO
  ENDDO

  IF(Associated(myFXR)) NULLIFY(myFXR)

  END SUBROUTINE

  SUBROUTINE MakeDeplXs1g(Fxr, ifxr, ng, phis, NormFactor, Core, ipin, iz, GroupInfo, PE)
  USE PARAM
  USE TypeDef,          ONLY : CoreInfo_Type,     FxrInfo_Type,      Cell_Type,       PE_TYPE,    &
                               GroupInfo_Type,    XsMac_Type
  USE MacXsLib_Mod,     ONLY : MacXsBase,  MacXsBase_Gen, MacXsAF, EffMacXs,  MacXS1gBase
  USE BasicOperation,   ONLY : CP_VA
  USE XsUtil_mod,       ONLY : FreeXsMac, FreeXsIsoMac
  USE nuclidmap_mod,    ONLY : iposiso
  USE Timer,            ONLY : nTracer_dclock, TimeChk
  IMPLICIT NONE
  TYPE(CoreInfo_Type) :: Core
  TYPE(FxrInfo_Type),POINTER :: Fxr(:,:)
  TYPE(FxrInfo_Type),POINTER :: myFxr
  TYPE(PE_TYPE) :: PE
  TYPE(GroupInfo_Type) :: GroupInfo
  REAL :: NormFactor
  INTEGER :: ng, ipin, iz, ifxr

  REAL :: phis(:)
  REAL, POINTER :: xs1g(:,:)
  REAL(8), ALLOCATABLE :: xsa(:), xsf(:), xsn2n(:), xsn3n(:)
  REAL, POINTER :: IsoXsMacf(:, :), IsoXsMacA(:, :)
  REAL :: Temp
  REAL :: Phi1g
  INTEGER, POINTER :: idiso(:)
  REAL, POINTER :: pnum(:)
  INTEGER :: nFsrInFxr, niso, iResoGrp1, iResoGrp2
  INTEGER :: ifsr, ig
  INTEGER:: i, j
  LOGICAL :: FlagF

  TYPE(XsMac_Type) :: XSMac
  TYPE(DeplFxr_Type), POINTER :: mydFxr
  INTEGER, POINTER :: IdIsoEig(:)
  REAL(8) :: frac, xsthis, sigcapm, sigcap, signa, signp, sign2n, sign2nm

  REAL :: Tb, Te
  REAL(4), POINTER :: fresoAIso(:,:), fresoFIso(:,:)

  myFxr => Fxr(ifxr,iz)
  mydFxr => DeplFxrBundle%FxrBundle(DeplFxrMap(ifxr,iz))

  IF (.NOT.mydFxr%ldepl) RETURN

  !Resonance Range
  iresoGrp1 = GroupInfo%nofg + 1; iresoGrp2 = GroupInfo%nofg + GroupInfo%norg
#define Depl_mod1
#ifdef Depl_mod0
  Temp = myFxr%Temp; niso = myFxr%niso_past
  idiso => myFxr%idiso_past; pnum => myFxr%pnum_past
#endif
#ifdef Depl_mod1
  Temp = myFxr%Temp; niso = myFxr%niso
  idiso => myFxr%idiso; pnum => myFxr%pnum
#endif

  IdIsoEig => mydFxr%IdIsoEig

  !TempSubgrp = Fxr%SubGrpTemp
  !1G XS
  phi1g = SUM(Phis(1:ng))

  xs1g => mydFxr%xs1g
  ALLOCATE(xsa(niso), xsf(niso), xsn2n(niso), xsn3n(niso))
  xsa = 0; xsf = 0; xsn2n = 0; xsn3n = 0;
  Tb = nTracer_dclock(FALSE, FALSE)
#ifdef OPTDBG
  CALL MacXS1gBase(xsa, xsf, xsn2n, xsn3n, Temp, niso, idiso(1:niso), pnum(1:niso), phis(1:ng), ng)
#else
  DO i = 1, niso
    CALL MacXS1gBase(xsa(i:i), xsf(i:i), xsn2n(i:i), xsn3n(i:i), Temp, 1, idiso(i:i), pnum(i:i), phis(1:ng), ng)
  ENDDO
#endif
  Te = nTracer_dclock(FALSE, FALSE)
  TimeChk%DeplBase1GTime = TimeChk%DeplBase1GTime + (Te-Tb)

  Tb = Te
  ALLOCATE(IsoXsMacA(niso, ng), IsoXsMacf(niso,ng))
  CALL MacXsAF(IsoXsMacA,IsoXsMacf,Temp,niso,idiso,pnum,1, ng, ng)
  Te = nTracer_dclock(FALSE, FALSE)
  TimeChk%DeplBaseMGTime = TimeChk%DeplBaseMGTime + (Te-Tb)
  !fresoAIso = myFxr%fresoAIso; fresoFIso => myFxr%fresoFIso;
  IF(myFxr%lRes) THEN
    DO ig = iResoGrp1, iResoGrp2
          DO i = 1, niso
              IsoXsMacF(i,ig) = IsoXsMacF(i,ig) * myFxr%fresoFIso(i,ig)
              IsoXsMacA(i,ig) = IsoXsMacA(i,ig) * myFxr%fresoAIso(i,ig)
          ENDDO
    ENDDO
  ENDIF
  DO i = 1, niso
    IF(pnum(i) .GT. epsm20) THEN
      xsa(i) = 0; xsf(i) = 0
      DO ig = 1, ng
        xsa(i) = xsa(i) + Phis(ig) * IsoXsMacA(i, ig)
        xsf(i) = xsf(i) + Phis(ig) * IsoXsMacF(i, ig)
      ENDDO
      xsa(i) = xsa(i) / phi1g; xsf(i) = xsf(i) / phi1g
    ENDIF
    xsa(i) = xsa(i) / pnum(i); xsf(i) = xsf(i) / pnum(i)
    xsn2n(i) = xsn2n(i) / pnum(i); xsn3n(i) = xsn3n(i) / pnum(i)
    xsa(i)  = xsa(i)  + xsn2n(i) + 2.*xsn3n(i)
  ENDDO
  mydFxr%NormFlux1g = phi1g * NormFactor * 1.E-24_8

  ! From TranDepXS2DeplLib - 19/01/24 eidt by LHG
  DO i = 1, niso
    j = idiso(i);
    j = iposiso(j);
    j = MapXs2Dep(j)
    !j = IdIsoEig(i)

!    if (ANY(xsa.LT.0)) print*, 'XSA Negative', idiso(i)
!    if (ANY(xsf.LT.0)) print*, 'XSF Negative', idiso(i)
!    if (ANY(xsn2n.LT.0)) print*, 'XSN2N Negative', idiso(i)
!    if (ANY(xsn3n.LT.0)) print*, 'XSN3N Negative', idiso(i)

    IF (j.LT.1 .OR. j.GT.DeplLib%NofIso) CYCLE
    sigcap = DeplLib%RctXs(RctIdCap,j); sigcapm = DeplLib%RctXs(RctIdCapm,j)
    sign2n = DeplLib%RctXs(RctIdN2N,j); sign2nm = DeplLib%RctXs(RctIdN2Nm,j)
    IF (.NOT. DeplLib%lActinide(j)) THEN
      sigNA = DeplLib%RctXs(RctIdNA,j); sigNP = DeplLib%RctXs(RctIdNP,j)
    END IF
    ! Capture, Fission, (n, 2n), and Minor Reactions
    IF ((sigcap+sigcapm).GT.1.D-30 .AND. (xsa(i)-xsf(i)).GT.1.D-30) THEN
      IF (DeplLib%lActinide(j)) THEN
        frac = sigcap/(sigcap+sigcapm); xsthis = xsa(i)-xsf(i)

        xs1g(RctIdCap,j) = frac*xsthis
        xs1g(RctIdCapm,j) = (1.-frac)*xsthis

        xs1g(RctIdN3N,j) = xsn3n(i); xs1g(RctIdFis,j) = xsf(i)
      ELSE
        frac = 1./(sigcap+sigcapm+signa+signp); xsthis = xsa(i)
        frac = frac*xsthis

        xs1g(RctIdCap,j) = sigcap*frac; xs1g(RctIdCapm,j) = sigcapm*frac
        xs1g(RctIdNA,j) = signa*frac; xs1g(RctIdNP,j) = signp*frac
      END IF
    END IF
    ! (n,2n)
    IF ((sign2n+sign2nm).GT.1.D-30 .AND. xsn2n(i).GT.1.D-30) THEN
      frac = sign2n/(sign2n+sign2nm); xsthis = xsn2n(i)
      xs1g(RctIdN2N,j) = frac*xsthis
      xs1g(RctIdN2Nm,j) = (1.-frac)*xsthis
    END IF
  END DO


  NULLIFY(idiso, pnum)
  DEALLOCATE(xsa, xsf, xsn2n, xsn3n)
  DEALLOCATE(IsoXsMacA, IsoXsMacf)
  END SUBROUTINE MakeDeplXs1g

  SUBROUTINE DepletionStep(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, PE, lSavePre)
  USE PARAM
  USE TYPEDEF,          ONLY : CoreInfo_Type,       FmInfo_Type,       GroupInfo_Type,    &
                               PE_Type,             ThInfo_Type,                          &
                               FxrInfo_Type,        Pin_Type,          Cell_Type
  USE CNTL,             ONLY : nTracerCntl_Type
  USE Depl_Mod,         ONLY : FluxNormalizeFactor
  USE Timer,            ONLY : nTracer_dclock, TimeChk
  USE XSLIB_MOD,        ONLY : ldiso
  IMPLICIT NONE
  TYPE(CoreInfo_Type) :: Core
  TYPE(FmInfo_Type) ::  FmInfo
  TYPE(ThInfo_Type) :: ThInfo
  TYPE(GroupInfo_Type) :: GroupInfo
  TYPE(nTracerCntl_Type) :: nTracerCntl
  TYPE(PE_TYPE) :: PE
  LOGICAL :: lSavePre

  TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
  !TYPE(FxrInfo_Type), POINTER :: myFXR
  TYPE(Pin_Type), POINTER :: Pin(:)
  TYPE(Cell_Type), POINTER :: CellInfo(:)

  REAL, POINTER :: phis(:, :, :), SpecConv(:)

  REAL, POINTER :: AvgPhi(:)
  REAL :: NormFactor, BurnUpTime

  INTEGER :: nxy, nFxr, nFsr, ntiso, myzb, myze, ng
  INTEGER :: nIsoLib, nIsoDepl
  INTEGER :: FsrIdxSt, FxrIdxSt, nLocalFxr, nLocalFsr, nFsrInFxr

  LOGICAL :: lPredictStep
  LOGICAL :: lCorrectStep
  LOGICAL :: lCritSpec
  LOGICAL :: lHighOd, lGdQuadDepl
  LOGICAL :: lfirst

  REAL :: Tbeg, Tend, Telapsed

  lPredictStep = DeplCntl%lPredict
  lCorrectStep = .NOT. lPredictStep
  lCritSpec = nTracerCntl%lCritSpec
  lGdQuadDepl = .FALSE.
  lfirst = (DeplCntl%NowStep .EQ. 1)
  !lHighOd = .TRUE.

  IF (DeplCntl%NowStep .GT. 2) lGdQuadDepl = .TRUE.

  Fxr=>FmInfo%Fxr; myzb = PE%myzb; myze = PE%myze

  !IF (DeplFxrBundle%nTrueDepl .EQ. 0) RETURN
  Tbeg = nTracer_dclock(FALSE, FALSE)
  CALL PreSetDeplSys(Core, FmInfo, GroupInfo, nTracerCntl, PE)
  CALL SetGdVars(DeplFxrBundle, lCorrectStep, nDeplTh, ldiso(:)%crit_nd)
  Tend = nTracer_dclock(FALSE, FALSE)
  TimeChk%DeplSetTime = TimeChk%DeplSetTime+Tend-Tbeg
  Tbeg = Tend
  CALL FxrBurnUp(lSavePre, lGdQuadDepl)
  CALL PostCorrection(DeplFxrBundle, lCorrectStep, lfirst, nDeplTh)
  Tend = nTracer_dclock(FALSE, FALSE)
  TimeChk%DeplSolTime = TimeChk%DeplSolTime+Tend-Tbeg
  Tbeg = Tend
  CALL UpdateDeplFxrInfo(Fxr, myzb, myze, GroupInfo, lCorrectStep, DeplCntl%lXeDyn)
  Tend = nTracer_dclock(FALSE, FALSE)
  TimeChk%DeplPostTime = TimeChk%DeplPostTime+Tend-Tbeg

  END SUBROUTINE DepletionStep

  SUBROUTINE PreSetDeplSys(Core, FmInfo, GroupInfo, nTracerCntl, PE)
    USE PARAM
    USE TYPEDEF,          ONLY : CoreInfo_Type,       FmInfo_Type,       GroupInfo_Type,    &
                                 PE_Type,                                                   &
                                 FxrInfo_Type,        Pin_Type,          Cell_Type
    USE CNTL,             ONLY : nTracerCntl_Type
    USE Depl_mod,         ONLY : FluxNormalizeFactor
    USE BasicOperation,   ONLY : MULTI_VA
    USE MOC_Mod,          ONLY : FxrAvgPhi
    USE OMP_LIB
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    TYPE(FmInfo_Type) ::  FmInfo
    TYPE(GroupInfo_Type) :: GroupInfo
    TYPE(nTracerCntl_Type) :: nTracerCntl
    TYPE(PE_TYPE) :: PE

    TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
    TYPE(Pin_Type), POINTER :: Pin(:)
    TYPE(Cell_Type), POINTER :: CellInfo(:)

    REAL, POINTER :: phis(:, :, :), SpecConv(:)

    REAL, POINTER :: AvgPhi(:)
    REAL :: NormFactor, BurnUpTime

    INTEGER :: nxy, nFxr, nFsr, ntiso, myzb, myze, ng
    INTEGER :: nIsoLib, nIsoDepl
    INTEGER :: FsrIdxSt, FxrIdxSt, nLocalFxr, nLocalFsr, nFsrInFxr

    LOGICAL :: lPredictStep
    LOGICAL :: lCorrectStep
    LOGICAL :: lCritSpec
    LOGICAL :: lHighOd, lGdQuadDepl

    INTEGER :: iz, ipin, icel, ifxr, tid
    INTEGER :: i, j

    lPredictStep = DeplCntl%lPredict
    lCorrectStep = .NOT. lPredictStep
    lCritSpec = nTracerCntl%lCritSpec
    !lHighOd = .TRUE.
    Pin => Core%Pin; CellInfo => Core%CellInfo;
    Fxr => FmInfo%Fxr; Phis => FmInfo%Phis

    ng = GroupInfo%ng
    nxy = Core%nxy;  ntiso = GroupInfo%ntiso
    nFsr = Core%nCoreFsr; nFxr = Core%nCoreFxr

    myzb = PE%myzb; myze = PE%myze


    NormFactor = FluxNormalizeFactor(Core, FmInfo, GroupInfo, DeplCntl%PowerCore, lCritSpec, TRUE, PE)
    NormFactor = NormFactor * nTracerCntl%PowerLevel
    CALL SetLocalBurnup(Core, FmInfo%FXR, FmInfo%Power, NormFactor, DeplCntl%Tsec, lCorrectStep, PE)

    IF(lCritSpec) SpecConv => FmInfo%SpecConv

    call OMP_SET_NUM_THREADS(nDeplTh);
    !$OMP PARALLEL DEFAULT(SHARED)&
    !$OMP PRIVATE(ipin,icel,j,ifxr,nFsrInFxr,AvgPhi,FsrIdxSt,FxrIdxSt,nLocalFxr,nLocalFsr)
    ALLOCATE(avgphi(ng))
    !$OMP DO SCHEDULE(GUIDED)
    DO iz = myzb, myze
      DO ipin = 1, nxy
        FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
        icel = Pin(ipin)%Cell(iz);
        nlocalFxr = CellInfo(icel)%nFxr; nlocalFsr = CellInfo(icel)%nFsr
        !Save Isotope list and number density prior to run depletion calculation
        IF(lPredictStep) CALL SaveFxrIsoInfo(Fxr(FxrIdxSt : FxrIdxst + nLocalFxr - 1, iz), nLocalFxr)

        DO j = 1, nLocalFxr
          ifxr = FxrIdxSt + j -1
          IF(.NOT. Fxr(ifxr, iz)%lDepl) CYCLE
          nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)

          AvgPhi = FxrAvgPhi(Core, Fxr, Phis, ipin, j, iz, ng, PE)
          !Conversion Critical Spectrum
          IF(lCritSpec) CALL MULTI_VA(SpecConv(1:ng),AvgPhi, ng)
          !Calculate 1G XS for Depletion Calculation
          CALL ConstDeplVasFirst(Fxr(ifxr,iz), ifxr, iz, ntiso, DeplLib%NofIso)
          CALL MakeDeplXs1g(Fxr, ifxr, ng, AvgPhi, NormFactor, Core, ipin, iz, GroupInfo, PE)
        ENDDO
      ENDDO
    ENDDO
    !$OMP END DO
    DEALLOCATE(avgphi)
    !$OMP END PARALLEL
    IF(lCritSpec) NULLIFY(SpecConv)
  END SUBROUTINE PreSetDeplSys

  SUBROUTINE FxrBurnUp(lSavePre, lQuad)
  USE PARAM
  USE CUDA_MASTER, ONLY : cuDevice
  USE CUDA_Workspace, ONLY : SetMatInfo, DeleteBatchMem, DeleteBlockMem
  USE cuMatExponential, ONLY : IsNaNval
  USE Timer,            ONLY : nTracer_dclock, TimeChk
  IMPLICIT NONE
  LOGICAL :: lSavePre
  LOGICAL :: lCorrector, lQuad
  LOGICAL :: lGd, IsCRAM

  INTEGER :: nfxr, NBun, NsysB, Nsys
  INTEGER :: i, ifxrbeg

  INTEGER :: NofIso, nIsoGd, i155
  INTEGER :: j, k, l
  REAL(8) :: f, xsabs, delT
  REAL(8), ALLOCATABLE :: pnums(:), pnums_RK(:,:), pnums_Gd155(:,:)
  INTEGER, POINTER :: IdIsoGd(:), IdTrueGd(:)

  TYPE(GdFxrInfo_Type), POINTER :: GdFxrs(:)
  TYPE(DeplFxr_Type), POINTER :: aFxr

  REAL(8), POINTER :: Nvec(:), SolVec(:)
  INTEGER :: Nsys0, ifxrbeg0

  REAL :: Tb, Te

  LOGICAL :: IsNaNHost(16384)

  lCorrector = (.NOT. DeplCntl%lPredict)
  IsCRAM  = (DeplSolTyp .EQ. CRAMSolTyp)

  DeplFxrBundle%delT = DeplCntl%Tsec

  lGd = .FALSE.
  nfxr = DeplFxrBundle%nTrueDepl
  IF (nfxr.LT.1) RETURN
  Nsys = CalSizeSysBun(DeplLib, SysByte, Scale, nfxr, IsCRAM)
  Nsys = MAX(Nsys,1)
  NBun = nfxr/Nsys
!  print*, Nsys, NBun
  IF (NBun .LT. 1) NBun = 1
  NsysB = FLOOR(dble(nfxr)/dble(NBun))

  NofIso = DeplLib%NofIso
  CALL SetMatInfo(SIZE(DeplLib%YldMapColIdx), NofIso, DeplLib%YldMapRowPtr, DeplLib%YldMapColIdx)

  ifxrbeg = 1
  DO i = 1, NBun
    IF (i .EQ. NBun) NsysB = nfxr-NsysB*(NBun-1)
    IF (i .NE. 1) CALL DestroyVecs_wCopy(DeplFxrBundle, ifxrbeg0, Nsys0, NofIso, Nvec, Solvec, lCorrector, .FALSE., .FALSE., cuDevice%mystream)
    Tb = nTracer_dclock(FALSE, FALSE)
    CALL SetDeplSys_woCSRT(DeplLib, DeplFxrBundle, NsysB, ifxrbeg, DeplSysBundle, lGd, (lQuad.AND.lCorrector), nDeplTh)
    Te = nTracer_dclock(FALSE, FALSE)
    TimeChk%DeplSysTime = TimeChk%DeplSysTime + (Te-Tb)
    TimeChk%DeplSolTime = TimeChk%DeplSolTime - (Te-Tb)
    CALL cuSolveDeplSys(DeplSysBundle, DeplFxrBundle, lCorrector, lGd, DeplSolTyp, nDeplTh)
    Tb = nTracer_dclock(FALSE, FALSE)
    CALL DestroySys_wPoint(DeplSysBundle, Nvec, Solvec, .FALSE.)
    Te = nTracer_dclock(FALSE, FALSE)
    TimeChk%DeplSysTime = TimeChk%DeplSysTime + (Te-Tb)
    TimeChk%DeplSolTime = TimeChk%DeplSolTime - (Te-Tb)
    !CALL DestroySys(DeplSysBundle)
    IsNaNHost = IsNaNval
    IF (ANY(IsNaNHost)) THEN
      DO j = 1, 16384
        IF (IsNaNHost(j)) THEN
          WRITE(*, '(A,I2,A,I6,A)') "NaN during CRAM at ", i, "th batch, and ", j,"th fxr"
          STOP
        END IF
      END DO
    END IF
    Nsys0 = NsysB; ifxrbeg0 = ifxrbeg
    ifxrbeg = ifxrbeg+NsysB
  END DO
  CALL DestroyVecs_wCopy(DeplFxrBundle, ifxrbeg0, Nsys0, NofIso, Nvec, Solvec, lCorrector, .FALSE., .FALSE., cuDevice%mystream)

#ifdef GDQUAD_SIMPLE
  IF (lQuad .AND. lCorrector .AND. (DeplFxrBundle%nTrueGd.GT.0)) THEN
!    print*, "WOWOWOWOWOWOW"
    lGd = .TRUE.
    NofIso = DeplLib%NofIso; nIsoGd = DeplFxrBundle%nIsoGd; i155 = DeplFxrBundle%i155
    IdIsoGd => DeplFxrBundle%IdIsoGd; IdTrueGd => DeplFxrBundle%IdTrueGd

    nfxr = DeplFxrBundle%nTrueGd
    ALLOCATE(pnums_RK(nfxr*7, 4)); pnums_RK = 0.;
    ALLOCATE(pnums_Gd155(nfxr, 4)); pnums_Gd155 = 0.;
    ALLOCATE(pnums(nfxr*7)); pnums = 0.;
    GdFxrs => DeplFxrBundle%GdFxrBundle(:)
    DO k = 1, nfxr
      pnums(1+(k-1)*7:7*k) = GdFxrs(k)%aFxr%pnum_depl(IdIsoGd(:))
    END DO
    DO k = 1, 4
      pnums_RK(:,k) = pnums
    END DO
    DO j = 1, nSubStep
      DO k = 1, nfxr
        pnums_Gd155(k,1) = pnums(7*(k-1)+3)
      END DO
      DO k = 1, nfxr
        aFxr => GdFxrs(k)%aFxr
        DO l = 1, nIsoGd
          xsabs = QuadFunc(GdFxrs(k)%c_qd(:,l),pnums_Gd155(k,1))
          IF (.NOT.(xsabs .GE. 1.e-40)) CYCLE
          f = xsabs
          xsabs = aFxr%xs1g(RctIdCAP,IdIsoGd(l))+aFxr%xs1g(RctIdCAPm,IdIsoGd(l))+&
            aFxr%xs1g(RctIdNP,IdIsoGd(l))+aFxr%xs1g(RctIdNA,IdIsoGd(l))
          f = f/xsabs
          aFxr%xs1g(RctIdCAP,IdIsoGd(l)) = aFxr%xs1g(RctIdCAP,IdIsoGd(l))*f
          aFxr%xs1g(RctIdCAPm,IdIsoGd(l)) = aFxr%xs1g(RctIdCAPm,IdIsoGd(l))*f
          aFxr%xs1g(RctIdNP,IdIsoGd(l)) = aFxr%xs1g(RctIdNP,IdIsoGd(l))*f
          aFxr%xs1g(RctIdNA,IdIsoGd(l)) = aFxr%xs1g(RctIdNA,IdIsoGd(l))*f
        END DO
      END DO
      ! ------------------------------------------ 1st --------------------------------------------
      CALL SolveAnalyticGd(DeplLib,DeplFxrBundle,pnums_RK(:,1),nDeplTh,nSubStep)
      DO k = 1, nfxr
        pnums_Gd155(k,2) = pnums_RK(7*(k-1)+3,1)
        pnums_Gd155(k,2) = 0.5*(pnums_Gd155(k,2)+pnums_Gd155(k,1))
      END DO

       DO k = 1, nfxr
        aFxr => GdFxrs(k)%aFxr
        DO l = 1, nIsoGd
          xsabs = QuadFunc(GdFxrs(k)%c_qd(:,l),pnums_Gd155(k,2))
          IF (.NOT.(xsabs .GE. 1.e-40)) CYCLE
          f = xsabs
          xsabs = aFxr%xs1g(RctIdCAP,IdIsoGd(l))+aFxr%xs1g(RctIdCAPm,IdIsoGd(l))+&
            aFxr%xs1g(RctIdNP,IdIsoGd(l))+aFxr%xs1g(RctIdNA,IdIsoGd(l))
          f = f/xsabs
          !IF (.NOT.(f.LT.1.5 .AND. f.GT.0.5)) print*, xsabs, f*xsabs, GdFxrs(k)%GdRR(:,l), GdFxrs(k)%Gd155(:)
          aFxr%xs1g(RctIdCAP,IdIsoGd(l))  = aFxr%xs1g(RctIdCAP,IdIsoGd(l))*f
          aFxr%xs1g(RctIdCAPm,IdIsoGd(l)) = aFxr%xs1g(RctIdCAPm,IdIsoGd(l))*f
          aFxr%xs1g(RctIdNP,IdIsoGd(l))   = aFxr%xs1g(RctIdNP,IdIsoGd(l))*f
          aFxr%xs1g(RctIdNA,IdIsoGd(l))   = aFxr%xs1g(RctIdNA,IdIsoGd(l))*f
        END DO
      END DO
      ! ------------------------------------------ 2nd --------------------------------------------
      CALL SolveAnalyticGd(DeplLib,DeplFxrBundle,pnums_RK(:,2),nDeplTh,nSubStep)
      DO k = 1, nfxr
        pnums_Gd155(k,3) = pnums_RK(7*(k-1)+3,2)
        pnums_Gd155(k,3) = 0.5*(pnums_Gd155(k,3)+pnums_Gd155(k,1))
      END DO

      DO k = 1, nfxr
        aFxr => GdFxrs(k)%aFxr
        DO l = 1, nIsoGd
          xsabs = QuadFunc(GdFxrs(k)%c_qd(:,l),pnums_Gd155(k,3))
          IF (.NOT.(xsabs .GE. 1.e-40)) CYCLE
          f = xsabs
          xsabs = aFxr%xs1g(RctIdCAP,IdIsoGd(l))+aFxr%xs1g(RctIdCAPm,IdIsoGd(l))+&
            aFxr%xs1g(RctIdNP,IdIsoGd(l))+aFxr%xs1g(RctIdNA,IdIsoGd(l))
          f = f/xsabs
          !IF (.NOT.(f.LT.1.5 .AND. f.GT.0.5)) print*, xsabs, f*xsabs, GdFxrs(k)%GdRR(:,l), GdFxrs(k)%Gd155(:)
          aFxr%xs1g(RctIdCAP,IdIsoGd(l)) = aFxr%xs1g(RctIdCAP,IdIsoGd(l))*f
          aFxr%xs1g(RctIdCAPm,IdIsoGd(l)) = aFxr%xs1g(RctIdCAPm,IdIsoGd(l))*f
          aFxr%xs1g(RctIdNP,IdIsoGd(l)) = aFxr%xs1g(RctIdNP,IdIsoGd(l))*f
          aFxr%xs1g(RctIdNA,IdIsoGd(l)) = aFxr%xs1g(RctIdNA,IdIsoGd(l))*f
        END DO
      END DO
      ! ------------------------------------------ 3rd --------------------------------------------
      CALL SolveAnalyticGd(DeplLib,DeplFxrBundle,pnums_RK(:,3),nDeplTh,nSubStep)
      DO k = 1, nfxr
        pnums_Gd155(k,4) = pnums_RK(7*(k-1)+3,3)
      END DO

      DO k = 1, nfxr
        aFxr => GdFxrs(k)%aFxr
        DO l = 1, nIsoGd
          xsabs = QuadFunc(GdFxrs(k)%c_qd(:,l),pnums_Gd155(k,4))
          IF (.NOT.(xsabs .GE. 1.e-40)) CYCLE
          f = xsabs
          xsabs = aFxr%xs1g(RctIdCAP,IdIsoGd(l))+aFxr%xs1g(RctIdCAPm,IdIsoGd(l))+&
            aFxr%xs1g(RctIdNP,IdIsoGd(l))+aFxr%xs1g(RctIdNA,IdIsoGd(l))
          f = f/xsabs
          !IF (.NOT.(f.LT.1.5 .AND. f.GT.0.5)) print*, xsabs, f*xsabs, GdFxrs(k)%GdRR(:,l), GdFxrs(k)%Gd155(:)
          aFxr%xs1g(RctIdCAP,IdIsoGd(l)) = aFxr%xs1g(RctIdCAP,IdIsoGd(l))*f
          aFxr%xs1g(RctIdCAPm,IdIsoGd(l)) = aFxr%xs1g(RctIdCAPm,IdIsoGd(l))*f
          aFxr%xs1g(RctIdNP,IdIsoGd(l)) = aFxr%xs1g(RctIdNP,IdIsoGd(l))*f
          aFxr%xs1g(RctIdNA,IdIsoGd(l)) = aFxr%xs1g(RctIdNA,IdIsoGd(l))*f
        END DO
      END DO
      ! ------------------------------------------ 4th --------------------------------------------
      CALL SolveAnalyticGd(DeplLib,DeplFxrBundle,pnums_RK(:,4),nDeplTh,nSubStep)
      ! -------------------------------------------------------------------------------------------

      pnums(:) = pnums_RK(:,1)+2.*pnums_RK(:,2)+2.*pnums_RK(:,3)+pnums_RK(:,4)
      pnums(:) = pnums(:)/6.
      DO k = 1, nfxr*7
        IF (pnums(k) .LT. 1.e-30) pnums(k) = 0.
      END DO
    END DO
    DO l = 1, nfxr
      aFxr => GdFxrs(l)%aFxr
      aFxr%pnum_cor(IdIsoGd(:)) = pnums((l-1)*7+1:l*7)
    END DO
    DEALLOCATE(pnums, pnums_Gd155, pnums_RK)
  END IF
#else
  IF (lQuad .AND. lCorrector .AND. (DeplFxrBundle%nTrueGd.GT.0)) THEN
!    print*, "WOWOWOWOWOWOW"
    lGd = .TRUE.
    NofIso = DeplLib%NofIso; nIsoGd = DeplFxrBundle%nIsoGd; i155 = DeplFxrBundle%i155
    IdIsoGd => DeplFxrBundle%IdIsoGd; IdTrueGd => DeplFxrBundle%IdTrueGd

    nfxr = DeplFxrBundle%nTrueGd
    Nsys = CalSizeSysBun(DeplLib, SysByte, Scale, nfxr, IsCRAM)
    NBun = nfxr/Nsys
    IF (NBun .LT. 1) NBun = 1
    NsysB = FLOOR(dble(nfxr)/dble(NBun))
    ifxrbeg = 1
    DO i = 1, NBun
      IF (i .EQ. NBun) NsysB = nfxr-NsysB*(NBun-1)
      ALLOCATE(pnums_RK(NsysB*NofIso, 4)); pnums_RK = 0.;
      ALLOCATE(pnums_Gd155(NsysB, 4)); pnums_Gd155 = 0.;
      ALLOCATE(pnums(NsysB*NofIso)); pnums = 0.;
      GdFxrs => DeplFxrBundle%GdFxrBundle(ifxrbeg:ifxrbeg-1+NsysB)
      DO  k = 1, NsysB
        pnums(NofIso*(k-1)+1:NofIso*k) = GdFxrs(k)%aFxr%pnum_depl(:)
      END DO
      DO j = 1, nSubStep
        Tb = nTracer_dclock(FALSE, FALSE)
        DO k = 1, NsysB
          pnums_Gd155(k,1) = pnums(NofIso*(k-1)+i155)
        END DO
          DO k = 1, NsysB
            aFxr => GdFxrs(k)%aFxr
            DO l = 1, nIsoGd
              xsabs = QuadFunc(GdFxrs(k)%c_qd(:,l),pnums_Gd155(k,1))
              IF (.NOT.(xsabs .GE. 1.e-40)) CYCLE
              f = xsabs
              xsabs = aFxr%xs1g(RctIdCAP,IdIsoGd(l))+aFxr%xs1g(RctIdCAPm,IdIsoGd(l))+&
                aFxr%xs1g(RctIdNP,IdIsoGd(l))+aFxr%xs1g(RctIdNA,IdIsoGd(l))
              f = f/xsabs
              aFxr%xs1g(RctIdCAP,IdIsoGd(l)) = aFxr%xs1g(RctIdCAP,IdIsoGd(l))*f
              aFxr%xs1g(RctIdCAPm,IdIsoGd(l)) = aFxr%xs1g(RctIdCAPm,IdIsoGd(l))*f
              aFxr%xs1g(RctIdNP,IdIsoGd(l)) = aFxr%xs1g(RctIdNP,IdIsoGd(l))*f
              aFxr%xs1g(RctIdNA,IdIsoGd(l)) = aFxr%xs1g(RctIdNA,IdIsoGd(l))*f
            END DO
          END DO
        !END IF
        ! ------------------------------------------ 1st --------------------------------------------
        CALL SetDeplSys_woCSRT(DeplLib, DeplFxrBundle, NsysB, ifxrbeg, DeplSysBundle, lGd, nDeplTh, nSubStep)
        Te = nTracer_dclock(FALSE, FALSE);
        TimeChk%DeplSysTime = TimeChk%DeplSysTime + (Te-Tb);
        TimeChk%DeplSolTime = TimeChk%DeplSolTime - (Te-Tb);
        DeplSysBundle%pnums = pnums
        CALL cuSolveDeplSys(DeplSysBundle, DeplFxrBundle, lCorrector, lGd, deplsoltyp, nDeplTh)
        Tb = nTracer_dclock(FALSE, FALSE)
        CALL CopySolVec(DeplSysBundle, NofIso, pnums_RK(:,1), .FALSE., cuDevice%mystream)
        DO k = 1, NsysB
          pnums_Gd155(k,2) = pnums_RK(NofIso*(k-1)+i155,1)
          pnums_Gd155(k,2) = 0.5*(pnums_Gd155(k,2)+pnums_Gd155(k,1))
        END DO
        CALL DestroySysnVec(DeplSysBundle,.FALSE.)

         DO k = 1, NsysB
          aFxr => GdFxrs(k)%aFxr
          DO l = 1, nIsoGd
            xsabs = QuadFunc(GdFxrs(k)%c_qd(:,l),pnums_Gd155(k,2))
            IF (.NOT.(xsabs .GE. 1.e-40)) CYCLE
            f = xsabs
            xsabs = aFxr%xs1g(RctIdCAP,IdIsoGd(l))+aFxr%xs1g(RctIdCAPm,IdIsoGd(l))+&
              aFxr%xs1g(RctIdNP,IdIsoGd(l))+aFxr%xs1g(RctIdNA,IdIsoGd(l))
            f = f/xsabs
            !IF (.NOT.(f.LT.1.5 .AND. f.GT.0.5)) print*, xsabs, f*xsabs, GdFxrs(k)%GdRR(:,l), GdFxrs(k)%Gd155(:)
            aFxr%xs1g(RctIdCAP,IdIsoGd(l))  = aFxr%xs1g(RctIdCAP,IdIsoGd(l))*f
            aFxr%xs1g(RctIdCAPm,IdIsoGd(l)) = aFxr%xs1g(RctIdCAPm,IdIsoGd(l))*f
            aFxr%xs1g(RctIdNP,IdIsoGd(l))   = aFxr%xs1g(RctIdNP,IdIsoGd(l))*f
            aFxr%xs1g(RctIdNA,IdIsoGd(l))   = aFxr%xs1g(RctIdNA,IdIsoGd(l))*f
          END DO
        END DO
        ! ------------------------------------------ 2nd --------------------------------------------
        CALL SetDeplSys_woCSRT(DeplLib, DeplFxrBundle, NsysB, ifxrbeg, DeplSysBundle, lGd, nDeplTh, nSubStep)
        Te = nTracer_dclock(FALSE, FALSE);
        TimeChk%DeplSysTime = TimeChk%DeplSysTime + (Te-Tb);
        TimeChk%DeplSolTime = TimeChk%DeplSolTime - (Te-Tb);
        DeplSysBundle%pnums = pnums
!        print*, "BOWOWOWOWOWOW2", NsysB
        CALL cuSolveDeplSys(DeplSysBundle, DeplFxrBundle, lCorrector, lGd, deplsoltyp, nDeplTh)
!        print*, "COWOWOWOWOWOW2"
        Tb = nTracer_dclock(FALSE, FALSE)
        CALL CopySolVec(DeplSysBundle, NofIso, pnums_RK(:,2), .FALSE., cuDevice%mystream)
        !pnums_RK(:,1) = DeplSysBundle%pnums1(:)
        DO k = 1, NsysB
          pnums_Gd155(k,3) = pnums_RK(NofIso*(k-1)+i155,2)
          pnums_Gd155(k,3) = 0.5*(pnums_Gd155(k,3)+pnums_Gd155(k,1))
        END DO
        CALL DestroySysnVec(DeplSysBundle,.FALSE.)

        DO k = 1, NsysB
          aFxr => GdFxrs(k)%aFxr
          DO l = 1, nIsoGd
            xsabs = QuadFunc(GdFxrs(k)%c_qd(:,l),pnums_Gd155(k,3))
            IF (.NOT.(xsabs .GE. 1.e-40)) CYCLE
            f = xsabs
            xsabs = aFxr%xs1g(RctIdCAP,IdIsoGd(l))+aFxr%xs1g(RctIdCAPm,IdIsoGd(l))+&
              aFxr%xs1g(RctIdNP,IdIsoGd(l))+aFxr%xs1g(RctIdNA,IdIsoGd(l))
            f = f/xsabs
            !IF (.NOT.(f.LT.1.5 .AND. f.GT.0.5)) print*, xsabs, f*xsabs, GdFxrs(k)%GdRR(:,l), GdFxrs(k)%Gd155(:)
            aFxr%xs1g(RctIdCAP,IdIsoGd(l)) = aFxr%xs1g(RctIdCAP,IdIsoGd(l))*f
            aFxr%xs1g(RctIdCAPm,IdIsoGd(l)) = aFxr%xs1g(RctIdCAPm,IdIsoGd(l))*f
            aFxr%xs1g(RctIdNP,IdIsoGd(l)) = aFxr%xs1g(RctIdNP,IdIsoGd(l))*f
            aFxr%xs1g(RctIdNA,IdIsoGd(l)) = aFxr%xs1g(RctIdNA,IdIsoGd(l))*f
          END DO
        END DO
        ! ------------------------------------------ 3rd --------------------------------------------
        CALL SetDeplSys_woCSRT(DeplLib, DeplFxrBundle, NsysB, ifxrbeg, DeplSysBundle, lGd, nDeplTh, nSubStep)
        Te = nTracer_dclock(FALSE, FALSE);
        TimeChk%DeplSysTime = TimeChk%DeplSysTime + (Te-Tb);
        TimeChk%DeplSolTime = TimeChk%DeplSolTime - (Te-Tb);
        DeplSysBundle%pnums = pnums
        CALL cuSolveDeplSys(DeplSysBundle, DeplFxrBundle, lCorrector, lGd, deplsoltyp, nDeplTh)
!        print*, "BOWOWOWOWOWOW3"
        Tb = nTracer_dclock(FALSE, FALSE)
        CALL CopySolVec(DeplSysBundle, NofIso, pnums_RK(:,3), .FALSE., cuDevice%mystream)
        !pnums_RK(:,1) = DeplSysBundle%pnums1(:)
        DO k = 1, NsysB
          pnums_Gd155(k,4) = pnums_RK(NofIso*(k-1)+i155,3)
        END DO
        CALL DestroySysnVec(DeplSysBundle,.FALSE.)

        DO k = 1, NsysB
          aFxr => GdFxrs(k)%aFxr
          DO l = 1, nIsoGd
            xsabs = QuadFunc(GdFxrs(k)%c_qd(:,l),pnums_Gd155(k,4))
            IF (.NOT.(xsabs .GE. 1.e-40)) CYCLE
            f = xsabs
            xsabs = aFxr%xs1g(RctIdCAP,IdIsoGd(l))+aFxr%xs1g(RctIdCAPm,IdIsoGd(l))+&
              aFxr%xs1g(RctIdNP,IdIsoGd(l))+aFxr%xs1g(RctIdNA,IdIsoGd(l))
            f = f/xsabs
            !IF (.NOT.(f.LT.1.5 .AND. f.GT.0.5)) print*, xsabs, f*xsabs, GdFxrs(k)%GdRR(:,l), GdFxrs(k)%Gd155(:)
            aFxr%xs1g(RctIdCAP,IdIsoGd(l)) = aFxr%xs1g(RctIdCAP,IdIsoGd(l))*f
            aFxr%xs1g(RctIdCAPm,IdIsoGd(l)) = aFxr%xs1g(RctIdCAPm,IdIsoGd(l))*f
            aFxr%xs1g(RctIdNP,IdIsoGd(l)) = aFxr%xs1g(RctIdNP,IdIsoGd(l))*f
            aFxr%xs1g(RctIdNA,IdIsoGd(l)) = aFxr%xs1g(RctIdNA,IdIsoGd(l))*f
          END DO
        END DO
        ! ------------------------------------------ 4th --------------------------------------------
        CALL SetDeplSys_woCSRT(DeplLib, DeplFxrBundle, NsysB, ifxrbeg, DeplSysBundle, lGd, nDeplTh, nSubStep)
        DeplSysBundle%pnums = pnums
        Te = nTracer_dclock(FALSE, FALSE);
        TimeChk%DeplSysTime = TimeChk%DeplSysTime + (Te-Tb);
        TimeChk%DeplSolTime = TimeChk%DeplSolTime - (Te-Tb);
        CALL cuSolveDeplSys(DeplSysBundle, DeplFxrBundle, lCorrector, lGd, deplsoltyp, nDeplTh)
!        print*, "BOWOWOWOWOWOW4"
        Tb = nTracer_dclock(FALSE, FALSE);
        CALL CopySolVec(DeplSysBundle, NofIso, pnums_RK(:,4), .FALSE., cuDevice%mystream)
        !pnums_RK(:,1) = DeplSysBundle%pnums1(:)
        !pnums_Gd155(k,1) = pnums_RK(NofIso*(k-1)+i155,4)
        CALL DestroySysnVec(DeplSysBundle,.FALSE.)
        ! -------------------------------------------------------------------------------------------

        pnums(:) = pnums_RK(:,1)+2.*pnums_RK(:,2)+2.*pnums_RK(:,3)+pnums_RK(:,4)
        pnums(:) = pnums(:)/6.
        DO k = 1, NsysB*NofIso
          IF (pnums(k) .LT. 1.e-30) pnums(k) = 0.
        END DO
        Te = nTracer_dclock(FALSE, FALSE);
        TimeChk%DeplSysTime = TimeChk%DeplSysTime + (Te-Tb);
        TimeChk%DeplSolTime = TimeChk%DeplSolTime - (Te-Tb);
        !IF (j.eq.1) print*, pnums
      END DO
      DO l = 1, NsysB
        aFxr => GdFxrs(l)%aFxr
        aFxr%pnum_cor(:) = pnums((l-1)*NofIso+1:l*NofIso)
      END DO
      ifxrbeg = ifxrbeg+NsysB
      DEALLOCATE(pnums, pnums_Gd155, pnums_RK)
    END DO
  END IF
#endif
  CALL UpdatePnumDepl(DeplFxrBundle, lCorrector, lSavePre)
  IF (lQuad .AND. lCorrector) CALL UpdatePnumDeplGd(DeplFxrBundle, lCorrector)
  CALL UpdatePnumSS(DeplFxrBundle, lCorrector)
  CALL DeleteBatchMem
  CALL DeleteBlockMem
  END SUBROUTINE FxrBurnUp

  SUBROUTINE UpdateDeplFxrInfo(Fxr, myzb, myze, GroupInfo, lCorrectStep, lXeDyn)
  USE PARAM
  USE TYPEDEF,            ONLY : FxrInfo_Type,       GroupInfo_Type
  USE BasicOperation,     ONLY : MULTI_CA,  CP_CA
  USE XsUtil_mod,         ONLY : SetXeDynEnv
  USE XSLIB_MOD,          ONLY : mapnucl
  USE nuclidmap_mod,      ONLY : iposiso,   PndCrit, nuclide
  IMPLICIT NONE
  TYPE(FxrInfo_Type), POINTER :: Fxr(:,:)
  INTEGER :: myzb, myze
  TYPE(GroupInfo_Type) :: GroupInfo
  LOGICAL :: lCorrectStep, lXeDyn

  REAL, POINTER :: IsoNum(:), pnum(:), pnum_all(:), pnum_depl(:), pnum_sseig(:)

  INTEGER :: iz, ifxr, i, j, IdXs, IdDep
  INTEGER :: nfxr, ntiso, niso, nisoDepl, tniso, niso_wolib
  INTEGER, POINTER :: IdIsoEig(:), IdIso(:)
  INTEGER :: IdIso_wolib(100)
  REAL :: pnum_wolib(100)
  TYPE(FxrInfo_Type), POINTER :: myFxr
  TYPE(DEplFxr_Type), POINTER :: mydFxr

  !INTEGER, POINTER :: CpIdIsoEig(:)
  !INTEGER, POINTER :: CpIdIso(:)

  ntiso = GroupInfo%ntiso
  nisoDepl = DeplLib%NofIso
  !nfxr = DeplFxrBundle%nfxr
  nfxr = SIZE(DeplFxrMap,1)

  DO iz = myzb, myze
    DO ifxr = 1, nfxr
      myFxr => Fxr(ifxr, iz)
      IF (DeplFxrMap(ifxr,iz) .EQ. 0) CYCLE
      IF (.NOT. myFxr%lDepl) CYCLE
      mydFxr => DeplFxrBundle%FxrBundle(DeplFxrMap(ifxr,iz))
      IdIso => myFxr%IdIso
      IF (.NOT. mydFxr%ldepl) CYCLE

      pnum => myFxr%pnum
      pnum_all => myFxr%pnum_all
      IF(DeplCntl%lPredict) THEN
        pnum_depl => mydFxr%pnum_pre
      ELSE
        pnum_depl => mydFxr%pnum_depl
      END IF
      pnum_sseig => mydFxr%pnum_sseig

      !pnum = mydFxr%pnum_sseig
      pnum_all = mydFxr%pnum_depl
      IdIsoEig => mydFxr%IdIsoEig(:)

      niso_wolib = 0
      DO i = 1, myFxr%niso
        j = IdIso(i)
        j = iposiso(j)
        j = MapXs2Dep(j)
        IF (j.EQ.0) THEN
          niso_wolib = niso_wolib+1
          IdIso_wolib(niso_wolib) = IdIso(i)
          pnum_wolib(niso_wolib) = pnum(i)
        END IF
      END DO

      !ALLOCATE(CpIdIsoEig(ntiso), CpIdIso(ntiso))
      !CpIdIsoEig = 0; CpIdIso = 0;
      !CpIdIso = IdIso
      !CALL CP_CA(IdIso, 0, ntiso)
      !CALL CP_CA(pnum, 0._8, ntiso)
      IdIso(1:ntiso) = 0.;
      pnum(1:ntiso) = 0.;


      niso = 0
      DO i = 1, nIsoDepl
        pnum_all(i) = pnum_depl(i)
        j = MapDep2Xs(i)
        IF(j .EQ. 0) CYCLE
        IF(abs(pnum_depl(i)) .LT. epsm20) CYCLE
        IF(pnum_depl(i) .LT. pndcrit(j)) CYCLE
        niso = niso + 1
        IdXs = nuclide(j)
        pnum(niso) = pnum_sseig(j)
        idiso(niso) = IdXs
        !CpIdIsoEig(niso) = i
      ENDDO

      DO i = 1, niso_wolib
        niso = niso+1
        pnum(niso) = pnum_wolib(i)
        Idiso(niso) = IdIso_wolib(i)
      END DO

      tniso = niso
      myFxr%niso = niso
      myFxr%niso_depl = niso
      !mydFxr%NisoEig = niso

      !DEALLOCATE(mydFxr%IdIsoEig)
      !ALLOCATE(mydFxr%IdIsoEig(1:niso))
      !mydFxr%IdIsoEig(1:niso) = CpIdIsoEig(1:niso)
      !DEALLOCATE(CpIdIsoEig, CpIdIso)

      IF(lXeDyn) CALL SetXeDynEnv(IdIso, pnum, myFxr%niso, myFxr%niso_depl, ntiso)
    END DO
  END DO

  END SUBROUTINE

  SUBROUTINE ConstDeplVasFirst(Fxr, ifxr, iz, nIsoLib, nIsoDepl)
  USE PARAM
  USE TYPEDEF,        ONLY : FxrInfo_type
  USE nuclidmap_mod,  ONLY : iposiso
  IMPLICIT NONE
  TYPE(FxrInfo_Type) :: Fxr
  INTEGER :: ifxr, iz
  INTEGER :: nIsoLib, nIsoDepl
  INTEGER :: i, j, k, niso, l
  INTEGER, POINTER :: IdIso(:)!, idisoEig(:), cpidiso(:)
  REAL, POINTER :: pnum(:)

  !DO i = 1, nFxr
    IdIso => Fxr%IdIso_past; pnum => Fxr%pnum_past
    nIso = Fxr%nIso_past;
    IF (.NOT. Fxr%l_pnum_all) THEN
      !DEALLOCATE(DeplFxrBundle%FxrBundle(DeplFxrMap(ifxr,iz))%idisoEig)
      !ALLOCATE(cpidiso(nisolib)); cpidiso = 0;
      l =0;
      DO j = 1, nIso
        k = IdIso(j); k = iposIso(k)
        k = MapXS2Dep(k)
        IF (k .EQ. 0) CYCLE
        l = l+1
        Fxr%pnum_past_all(k) = pnum(j)
        !cpidiso(l) = k
      END DO
      !ALLOCATE(DeplFxrBundle%FxrBundle(DeplFxrMap(ifxr,iz))%idisoEig(l))
      !DeplFxrBUndle%FxrBundle(DeplFxrMap(ifxr,iz))%idisoEig = cpidiso(1:l)
      !DeplFxrBUndle%FxrBundle(DeplFxrMap(ifxr,iz))%NisoEig = l
      !DEALLOCATE(cpidiso)
!      DeplFxrBundle%FxrBundle(DeplFxrMap(ifxr,iz))%pnum_depl = Fxr%pnum_past_all;
      Fxr%l_pnum_all = .TRUE.
    END IF
  !END DO

  END SUBROUTINE

  SUBROUTINE SaveFxrIsoInfo(Fxr, nFxr)
  USE PARAM
  USE TYPEDEF,     ONLY : FxrInfo_Type
  IMPLICIT NONE
  TYPE(FxrInfo_Type) :: Fxr(1:nFxr)
  INTEGER :: nFxr
  INTEGER :: ifxr, iz, niso

    DO ifxr = 1, nFxr
      IF(.NOT. Fxr(ifxr)%ldepl) CYCLE
      niso = Fxr(ifxr)%niso_depl
      Fxr(ifxr)%niso_past = niso
      Fxr(ifxr)%idiso_past(1:niso) = Fxr(ifxr)%idiso(1:niso)
      Fxr(ifxr)%pnum_past(1:niso) = Fxr(ifxr)%pnum(1:niso)
      IF(Fxr(ifxr)%l_pnum_all) Fxr(ifxr)%pnum_past_all = Fxr(ifxr)%pnum_all  ! 16/02/11 Depletion timestep bug fixed
    ENDDO
  END SUBROUTINE

  SUBROUTINE SetLocalBurnup(Core, FXR, Power, normalizer, Tsec, lCorrectStep, PE)
  USE PARAM
  USE TYPEDEF,       ONLY : CoreInfo_Type,   FxrInfo_Type, PE_Type,      &
                           Pin_Type,        Cell_Type
  IMPLICIT NONE

  TYPE(CoreInfo_Type) :: Core
  TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
  REAL, POINTER :: Power(:, :)
  REAL :: Normalizer, Tsec
  LOGICAL :: lCorrectStep
  TYPE(PE_TYPE) :: PE

  TYPE(Pin_Type), POINTER :: Pin(:)
  TYPE(Cell_Type), POINTER :: CellInfo(:)
  TYPE(FxrInfo_Type), POINTER :: myFxr

  REAL :: Area,norm_Mw, powersum, burnup, Tday
  INTEGER :: FsrIdxSt, FxrIdxSt
  INTEGER :: ifxr, ifsr, icel, iz, ipin, iLocalFsr, nFsrInFxr, nLocalFxr
  INTEGER :: i, j

  CellInfo => Core%CellInfo
  Pin => Core%Pin
  powersum = 0; Tday = Tsec / 86400._8
  norm_Mw = normalizer / 1.e+6_8
  DO iz = PE%myzb, PE%myze
    DO ipin = 1, Core%nxy
      FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
      icel = Pin(ipin)%Cell(iz)
      nLocalFxr = CellInfo(icel)%nFxr
      DO j =1, nLocalFxr
        ifxr = FxrIdxSt + j - 1; nFsrInFxr = Fxr(ifxr, iz)%nFsrInFxr
        IF(.NOT. Fxr(ifxr, iz)%lfuel) CYCLE
        burnup =0
        DO i = 1, nFsrInFxr
          iLocalFsr = CellInfo(icel)%MapFxr2FsrIdx(i, j)
          iFsr = FsrIdxSt + iLocalFsr - 1
          Area = CellInfo(icel)%vol(iLocalFsr);
          !FxrAvgPhi = FxrAvgPhi + Area * Phis(ifsr, iz, 1:ng)
          burnup = burnup + area * core%hz(iz) * power(ifsr, iz)
          powersum  = powersum + area * core%hz(iz) * power(ifsr, iz) * normalizer
        ENDDO
        burnup = burnup * norm_Mw * Tday / Fxr(ifxr, iz)%Hmkg0
        IF(.NOT. lCorrectStep) THEN
          Fxr(ifxr, iz)%burnup = Fxr(ifxr, iz)%burnup_past + burnup
        ELSE
          burnup = 0.5_8*(burnup + Fxr(ifxr, iz)%burnup - Fxr(ifxr, iz)%burnup_past)
          Fxr(ifxr, iz)%burnup = Fxr(ifxr, iz)%burnup_past + burnup
          Fxr(ifxr, iz)%burnup_past = Fxr(ifxr, iz)%burnup
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  NULLIFY(CellInfo, Pin)
  END SUBROUTINE

  FUNCTION GetCoreHmMass(Core, FXR, PE)
  USE PARAM
  USE TYPEDEF,       ONLY : CoreInfo_Type,   FxrInfo_TYPE,   PE_TYPE,    &
                            Pin_Type,        Cell_Type
  USE nuclidmap_mod, ONLY : iposnucl,        AtomicWeight
#ifdef MPI_ENV
  USE MPIComm_mod,   ONLY : REDUCE
#endif
  IMPLICIT NONE

  TYPE(CoreInfo_Type) :: Core
  TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
  TYPE(PE_TYPE) :: PE

  REAL :: GetCoreHmMass

  TYPE(Pin_Type), POINTER :: Pin(:)
  TYPE(Cell_Type), POINTER :: Cell(:)
  TYPE(FxrInfo_Type), POINTER :: myFxr

  REAL :: area, HmReg, HmCell, HmPln, HmCore, buf

  INTEGER :: iz, ixy, ifxr, icell, iso, id
  INTEGER :: i, j, k
  INTEGER :: nFSR, nFXR, nxy
  INTEGER :: FxrIdxSt, nLocalFxr, niso
  INTEGER :: myzb, myze

  Pin => Core%Pin; Cell => Core%CellInfo
  nFSR = Core%nCoreFsr; nFxr = Core%nCoreFxr
  nxy = Core%nxy
  myzb = PE%myzb; myze = PE%myze
  GetCoreHmMass = 0
  HmCore = 0
  DO iz = myzb, myze
    HmPln = 0
    DO ixy = 1, nxy
      icell = Pin(ixy)%cell(iz)
      FxrIdxSt = Pin(ixy)%FxrIdxSt; nLocalFxr = Cell(icell)%nFxr
      HmCell = 0
      DO i = 1, nLocalFxr
        ifxr = FxrIdxSt + i - 1
        IF(.NOT. Fxr(ifxr, iz)%lfuel) CYCLE
        myFxr => Fxr(ifxr, iz); niso = myFxr%niso
        area = myFxr%area; HmReg = 0
        DO iso = 1, niso
          id = myFxr%idiso(iso)
          IF(id .LT. 90000) CYCLE
          HmReg = HmReg + myFxr%pnum(iso) * AtomicWeight(id) / AVOGADRO
        ENDDO  !Isotope Sweep
        HmCell = HmCell + HmReg * Area
        !myFxr%Hmkg0 = HmReg * Area * Core%Hz(iz) / 1000._8
      ENDDO  !Local Fxr Sweep
      HmPln = HmPln + HmCell
    ENDDO
    HmCore = HmCore + HmPln * Core%Hz(iz)
  ENDDO

#ifdef MPI_ENV
  CALL REDUCE(HmCore, buf, PE%MPI_RTMASTER_COMM, .TRUE.)
  HmCore = buf
#endif

  GetCoreHmMass = HmCore
  GetCoreHmMass = GetCoreHmMass * epsm6    !Ton Unit Conversion

  END FUNCTION

  SUBROUTINE setXeDynVar()
  USE XsLib_mod, only : nelthel,ldiso
  USE nuclidmap_mod, ONLY : iposiso
  implicit none

  integer :: idXe, idI, iyXe, iyI, nid, i, ifis
  INTEGER :: NIsoDep, NofYld, NofFis

  NIsoDep = DeplLib%NofIso
  NofYld = DeplLib%NofYld
  NofFIs = DeplLib%NofFis

  idXe = 546350; idI = 536350

  idXe = IsoIdSrch(DeplLib%Idiso, NIsoDep, idXe)
  iyXe = IsoIdSrch(DeplLib%IdYld, NofYld, idXe)
  idI = IsoIdSrch(DeplLib%Idiso, NIsoDep, idI)
  iyI = IsoIdSrch(DeplLib%IdYld, NofYld, idI)

  allocate(yieldXe135(nelthel), yieldI135(nelthel))
  do i = 1, nelthel
      yieldXe135(i) = 0._8; yieldI135(i) = 0._8
      if (ldiso(i)%ifis.eq.0) cycle
      nid = ldiso(i)%nid
      nid = iposiso(nid)
      ifis = MapXS2Dep(nid)
      IF (ifis .LT. 1 .OR. ifis .GT. NIsoDep) CYCLE
      yieldXe135(i) = DeplLib%FisFrac(ifis,iyXe)
      yieldI135(i) = DeplLib%FisFrac(ifis,iyI)
  enddo
  decayxe135 = SUM(DeplLib%DecFrac(:,idXe))
  decayI135 = SUM(DeplLib%DecFrac(:,idI))
  END SUBROUTINE

  SUBROUTINE MakeXsDeplMap(GroupInfo)
  USE TYPEDEF, ONLY : GroupInfo_TYPE
  USE nuclidmap_mod, ONLY : iposiso
  IMPLICIT NONE
  TYPE(GroupInfo_Type) :: GroupInfo

  INTEGER :: ntiso, nisoDepl
  INTEGER, POINTER :: idiso(:)
  INTEGER :: ixs, i

  IF (lInitMapIso) RETURN

  ntiso = GroupInfo%ntiso
  nisoDepl = DeplLib%NofIso

  ALLOCATE(MapXs2Dep(ntiso), MapDep2XS(nisoDepl))
  MapXs2Dep = 0; MapDep2Xs = 0
  idiso => DeplLib%idiso

  DO i = 1, nisoDepl
    ixs = idiso(i)
    ixs = (ixs/10)+mod(ixs,10)*100
    ixs = iposiso(ixs)
    IF (ixs .NE. 0) THEN
      MapDep2XS(i) = ixs
      MapXS2Dep(ixs) = i
    END IF
  END DO
  lInitMapIso = .TRUE.

  END SUBROUTINE

  SUBROUTINE Depletion_Driver_CUDA()
  USE PARAM
  USE GEOM,          ONLY : Core, ng
  USE RAYS,          ONLY : RayInfo
  USE Core_mod,      ONLY : CMInfo,            FmInfo,               THInfo,        &
                            GroupInfo,         eigv,                 peigv,         &
                            xkconv,            xkcrit
  USE SubGrp_mod,    ONLY : SubGrpEffXsGen
  USE Th_Mod,        ONLY : SteadyStateTH,     THFeedBack,          ThVar
  USE Boron_mod,     ONLY : CalB10XS,          SaveBoronPPM,        PostB10Depl,   &
                            PostB10Depl0, ppmd
  USE Depl_Mod,      ONLY : SetBurnUpStepInfo,                                      &
                            SetDeplCoreState,  SaveCoreState,        EffBurnUpTime, &
                            SetDeplTimeStep,                                        &
                            UpdtFluxSolution,  DeplBanner
  USE CritSpec_mod,  ONLY : XsKinf
  USE XeDyn_Mod,     ONLY : SetDeplXeDynMode,  UpdtXeDyn,                           &
                            InitXe,            FinalXe,              UpdtXe
  USE PE_MOD,        ONLY : PE
  USE CNTL,          ONLY : nTracerCntl
  USE ItrCNTL_mod,   ONLY : ItrCntl
  USE FILES,         ONLY : io8
  USE IOUTIL,        ONLY : message,           ShowHbar1,            ShowHbar2
  USE Timer,         ONLY : nTracer_dclock, TimeChk
  USE Restart_mod,   ONLY : WriteRestartFile
#ifdef MPI_ENV
  USE MPIComm_mod,   ONLY : MPI_SYNC, MPI_MAX_REAL
#endif
  USE CUDA_Workspace, ONLY : WriteConstCRAM
  USE XS_COMMON
  USE AuxilDriver,    ONLY : CUDASubGrpEffXSGen
  USE SUBGRP_MOD,     ONLY : UpdtCoreIsoInfo
  IMPLICIT NONE

  REAL :: PrevBurnUp(2)
  REAL :: ppm
  REAL :: efftime, delpow
  REAL :: Tbeg, Tend, Telapsed
  INTEGER :: nBurnUpStep
  INTEGER :: iBurnUp, itype

  LOGICAL :: Master, Slave
  LOGICAL, SAVE :: lFirst


  Data lFirst /.TRUE./


  Master = PE%Master; Slave = PE%Slave

  WRITE(mesg, '(A)') 'Performing Depletion ...'
  IF(Master) CALL message(io8, TRUE, TRUE, mesg)

  IF(lFirst) THEN
    CALL WriteConstCRAM
    IF(.NOT. DeplCntl%lInitDepl) THEN
      DeplCntl%lInitDepl = .TRUE.
      CALL Init_Depl(PE)
      CALL MakeXsDeplMap(GroupInfo)
      CALL AllocDeplFxrMem(Core, FmInfo%Fxr, GroupInfo, PE)
      IF (nTracerCntl%lXsAlign) CALL AlignDeplFxrMem(FmInfo%Fxr, GroupInfo)
    ENDIF
    CALL DeplBanner(1)
    lFirst = .FALSE.
  ENDIF


  itype = DeplCntl%BurnUpType
  nBurnUpStep = DeplCntl%nBurnUpStep
  ppm = 0

  !Set Xenon Dynamics Mode
  CALL SetDeplXeDynMode(Core, FmInfo, DeplCntl, nTracerCntl, PE)

  CALL SetBurnUpStepInfo(nTracerCntl%PowerCore, DeplCntl)
  IF(DeplCntl%lCoreFollow) CALL InitCoreFollow(DeplCntl, nTracerCntl)

  !
  IF(DeplCntl%NowBurnUp(1) .LE. epsm8) THEN
    !Core%DelT = DeplCntl%TSec
    CALL SetDeplTimeStep(0, Core, DeplCntl, nTracerCntl ,PE)
    !Initialize Xenon Dynamics
    CALL UpdtXeDyn(Core, FmInfo, THInfo, GroupInfo, nTracerCntl, DeplCntl, PE, Mode = InitXe)
    !Solve Depletion Chain
    CALL DepletionStep(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, PE, .TRUE.)
    !Effective Cross Section Updates
    IF (nTracerCntl%lXsAlign) THEN
      CALL UpdtCoreIsoInfo(Core, FmInfo%Fxr, PE)
      CALL UpdateBlockFxr(FmInfo%Fxr, GroupInfo, DEPLFeed)
      CALL XsLinIntPolTemp(nTracerCntl%ScatOd)
!      CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, eigv, GroupInfo, nTracerCntl, PE)
      CALL CUDASubGrpEffXSGen(core,FmInfo%fxr,nTracerCntl,PE)
    ELSE
      CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, eigv, GroupInfo, nTracerCntl, PE)
    END IF
    !Update Xe Dynamics
    CALL UpdtXeDyn(Core, FmInfo, THInfo, GroupInfo, nTracerCntl, DeplCntl, PE, Mode = FinalXe)
  ENDIF

  IF(DeplCntl%B10DeplMod .NE. 0) CALL CalB10XS(Core, FmInfo, 0, GroupInfo%ng, GroupInfo, nTracerCntl, PE)

  DO iBurnUp = 1, nBurnUpStep
    nTracerCntl%CalNoId = nTracerCntl%CalNoId + 1
    IF(DeplCntl%lCoreFollow) THEN
      CALL SetDeplCoreState(iburnup, DeplCntl%CoreState, ThInfo, ThVar, nTracerCntl, PE)
      IF(DeplCntl%CoreState%lStateChg(iburnup)) CALL UpdtFluxSolution(iburnup, DeplCntl, nTracerCntl, PE)
    ENDIF
    !Burnup Time
    CALL SetDeplTimeStep(iburnup, Core,  DeplCntl, nTracerCntl, PE)

    WRITE(mesg, '(a)') 'Predictor Step...'
    IF(Master) CALL message(io8, TRUE, TRUE, mesg)

    !Time Check
    Tbeg = nTracer_dclock(FALSE, FALSE)


    !--- BYS edit 15/12/30 --- depletion test

    !Prediction Step
    DeplCntl%lPredict = .TRUE.
  !!#define H2H_DBG
#ifdef H2H_DBG
    WRITE(361,*) "Predictor Step...", iburnup
#endif
    !Effective Cross Section
    !CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, Eigv, GroupInfo, nTracerCntl, PE)
    !Depletion Matrix Calculation
!#define HGDBG
#ifdef HGDBG
    WRITE(120,*) 'Step',iBurnUp,'Predictor',DeplCntl%lPredict
#endif
    CALL DepletionStep(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, PE, .FALSE.)
    !Time Check
    Tend = nTracer_dclock(FALSE, FALSE)
    Telapsed = Tend-Tbeg
    CALL MPI_MAX_REAL(Telapsed, PE%MPI_RTMASTER_COMM, TRUE)
    TimeChk%DeplTime = TimeChk%DeplTime + Telapsed
    CALL MPI_MAX_REAL(TimeChk%DeplSetTime, PE%MPI_RTMASTER_COMM, TRUE)
    CALL MPI_MAX_REAL(TimeChk%DeplSysTime, PE%MPI_RTMASTER_COMM, TRUE)
    CALL MPI_MAX_REAL(TimeChk%DeplSolTime, PE%MPI_RTMASTER_COMM, TRUE)
    CALL MPI_MAX_REAL(TimeChk%DeplPostTime, PE%MPI_RTMASTER_COMM, TRUE)
    IF (nTracerCntl%lXsAlign) THEN
      CALL UpdtCoreIsoInfo(Core, FmInfo%Fxr, PE)
      CALL UpdateBlockFxr(FmInfo%Fxr, GroupInfo, DEPLFeed)
      CALL XsLinIntPolTemp(nTracerCntl%ScatOd)
    END IF

    IF(DeplCntl%B10DeplMod .NE. 0) CALL CalB10XS(Core, FmInfo, iBurnUp, GroupInfo%ng, GroupInfo, nTracerCntl, PE)
    IF(DeplCntl%B10DeplMod .EQ. 1 .AND.ppmd .GT. 0) CALL OnlineB10Depl(iBurnUp, .TRUE., DeplCntl, PE)

    !Flux Calculation
    CALL SSEIG
    IF (nTracerCntl%OutpCntl%nDeplRegXsOut.gt.0) CALL ProcessEffXsDepl(Core, FmInfo, CmInfo, ThInfo, GroupInfo, nTracerCntl, DeplCntl,PE,iburnup)

    IF(Master) CALL ShowHbar1(io8)

    !Correction Step
    WRITE(mesg, '(a)') 'Corrector Step'
    IF(Master) CALL message(io8, TRUE, TRUE, mesg)
    DeplCntl%lPredict = .FALSE.

    Tbeg = nTracer_dclock(FALSE, FALSE)
#ifdef HGDBG
    WRITE(120,*) 'Step',iBurnUp,'Predictor',DeplCntl%lPredict
#endif
    CALL DepletionStep(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, PE, .FALSE.)
    Tend = nTracer_dclock(FALSE, FALSE)
    Telapsed = Tend-Tbeg
    CALL MPI_MAX_REAL(Telapsed, PE%MPI_RTMASTER_COMM, TRUE)
    TimeChk%DeplTime = TimeChk%DeplTime + Telapsed

    CALL MPI_MAX_REAL(TimeChk%DeplSetTime, PE%MPI_RTMASTER_COMM, TRUE)
    CALL MPI_MAX_REAL(TimeChk%DeplSysTime, PE%MPI_RTMASTER_COMM, TRUE)
    CALL MPI_MAX_REAL(TimeChk%DeplSolTime, PE%MPI_RTMASTER_COMM, TRUE)
    CALL MPI_MAX_REAL(TimeChk%DeplPostTime, PE%MPI_RTMASTER_COMM, TRUE)
    !Changhyun Test - Number Density Out for Each Fxr
    !CALL Write_ND_Debug(Core, FmInfo, DeplCntl, PE)
    !Effective Cross Section
    IF (nTracerCntl%lXsAlign) THEN
      CALL UpdtCoreIsoInfo(Core, FmInfo%Fxr, PE)
      CALL UpdateBlockFxr(FmInfo%Fxr, GroupInfo, DEPLFeed)
      CALL XsLinIntPolTemp(nTracerCntl%ScatOd)
!      CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, eigv, GroupInfo, nTracerCntl, PE)
      CALL CUDASubGrpEffXSGen(core,FmInfo%fxr,nTracerCntl,PE)
    ELSE
      CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, eigv, GroupInfo, nTracerCntl, PE)
    END IF
    CALL UpdtXeDyn(Core, FmInfo, THInfo, GroupInfo, nTracerCntl, DeplCntl, PE, Mode = FinalXe)

    !--- BYS edit END 15/12/30 --- depletion test

    !2nd Flux Calculation
    IF(DeplCntl%PC_OPT .EQ. 2) CALL SSEIG() !Changhyun Test
    IF(DeplCntl%lCoreFollow .OR. DeplCntl%PC_OPT .EQ. 2) CALL SaveCoreState(iburnup, DeplCntl%CoreState, nTracerCntl, PE)

    xkcrit = XsKinf(Core, FmInfo, GroupInfo, .TRUE., nTracerCntl%lCritSpec, PE)
    CALL DeplBanner(2)
    IF(PE%lCmfdGrp) CALL OutputEdit()
    IF(nTRACERCntl%lWriteDeplRst) CALL  WriteRestartFile(Core, FmInfo, .TRUE., DeplCntl%NowStep, nTracerCntl, PE)
    IF(MASTER) CALL ShowHbar2(io8)
  ENDDO

  IF(DeplCntl%B10DeplMod .EQ. 2) CALL PostB10Depl(DeplCntl, PE)
  END SUBROUTINE
END MODULE

#endif
