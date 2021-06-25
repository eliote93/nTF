#include <defines.h>
SUBROUTINE SubGrpFsp(Core, Fxr, THInfo, RayInfo,  GroupInfo, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,      ONLY : CoreInfo_Type,    FxrInfo_Type,   THInfo_Type,  RayInfo_Type,  GroupInfo_Type,   PE_Type
USE CNTL,         ONLY : nTracerCntl_Type
USE SUBGRP_MOD,   ONLY : UpdtCoreIsoInfo, SubGrpFsp_CAT, SubGrpFsp_ISO, SubGrpFsp_MLG, CalcDancoff, CalcEscXSCP
USE TH_Mod,       ONLY : Cal_RefFuelTemp
USE FILES,        ONLY : IO8
USE IOUTIL,       ONLY : message, terminate
USE SubGrpFspNM,  ONLY : SubGrpFSP_MLG_NM   !--- CNJ Edit : Node Majors
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type),POINTER :: Fxr(:, :)
TYPE(THInfo_Type) :: THInfo
TYPE(RayInfo_Type) :: RayInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE

INTEGER :: iz,myzb,myze
LOGICAL :: lfuel, RTmaster, master

IF(.NOT. nTracerCntl%lXsLib) RETURN

myzb = PE%myzb; myze = PE%myze
lFuel=.false.
DO iz = myzb, myze
  lFuel = lFuel .OR. Core%lFuelPlane(iz)
ENDDO

#ifndef MPI_ENV
IF(.NOT. lFuel) RETURN
#endif

master = PE%master; RTmaster = PE%RTmaster

!Calculate the Reference Temperature
IF(RTmaster) THEN
  CALL Cal_RefFuelTemp(THInfo%RefFuelTemp, Core, Fxr, nTracerCntl, PE)
  CALL UpdtCoreIsoInfo(Core, Fxr, PE)
ENDIF

IF (nTracerCntl%lED) THEN
  CALL CalcDancoff(Core, Fxr, RayInfo, THInfo, nTracerCntl, PE)
  CALL CalcEscXSCP(Core, Fxr, GroupInfo, THInfo, nTracerCntl, PE)
ELSE
  IF (nTracerCntl%lMLG) THEN
    !--- CNJ Edit : Node Majors
    IF (nTracerCntl%lNodeMajor) THEN
      CALL SubGrpFsp_MLG_NM(Core, Fxr, THInfo, RayInfo, GroupInfo)
    ELSE
      CALL SubGrpFsp_MLG(Core, Fxr, THInfo, RayInfo,  GroupInfo, nTracerCntl, PE)
    ENDIF
  ELSE
    IF (nTracerCntl%lNodeMajor) CALL terminate("ONLY NODE MAJOR")
    
    IF (nTracerCntl%lCAT) THEN
      IF(RTmaster) CALL UpdtResoCat(PE)
      CALL SubGrpFsp_CAT(Core, Fxr, THInfo, RayInfo,  GroupInfo, nTracerCntl, PE)
    ELSE
      CALL SubGrpFsp_ISO(Core, Fxr, THInfo, RayInfo,  GroupInfo, nTracerCntl, PE)
    ENDIF
  ENDIF
ENDIF

END SUBROUTINE

SUBROUTINE SubGrpFsp_CAT(Core, Fxr, THInfo, RayInfo,  GroupInfo, nTracerCntl, PE)  !!!!!!!!!!****************************************!!!!!!!!!!
USE PARAM
USE TYPEDEF,        ONLY : coreinfo_type,    Fxrinfo_type,   GroupInfo_Type,    RayInfo_Type,   PE_TYPE, THInfo_Type
USE CNTL,           ONLY : nTracerCntl_Type
#ifdef MPI_ENV
USE MPICOMM_MOD,    ONLY : MPI_SYNC,    MPI_MAX_REAL,   MPI_MAX_INT,     BCAST
#endif
USE FILES,          ONLY : IO8
USE SUBGRP_MOD,     ONLY : SetPlnLsigP_cat,    SubGrpFspErr,    EquipXSGen,    SetSubGrpSrc1g, UpdtNDAF_CAT
USE MOC_MOD,        ONLY : RayTraceGM_OMP
USE BasicOperation, ONLY : CP_CA,                CP_VA
USE IOUTIL,         ONLY : message
USE TIMER,          ONLY : nTracer_dclock,       TimeChk
USE XSLib_mod,      ONLY : LIBDATA, nCat, ResoCat, mapnucl, ldiso, nActiveCat, ResoCatUse

IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type),POINTER :: Fxr(:, :)
TYPE(THInfo_Type) :: THInfo
TYPE(RayInfo_Type) :: RayInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE
TYPE(LIBDATA), POINTER :: isodata

INTEGER :: nfsr, nfxr, myzb, myze, iz, nlvflx, ilv, nid, id, icat
INTEGER :: ig, ig1, ig2, nofg, norg, niter, iter, nPhiAngSv, nPolar

REAL :: lv, errmax, errmaxlv, niteravg, Tbeg, Tend, rt1, rt2, rtnet
REAL, POINTER :: phis1g(:), phis1gd(:), PhiAngIn1g(:, :), jout(:,:,:)
REAL, POINTER :: siglamPot(:), xstr1g(:), src1g(:)

LOGICAL :: master, RTmaster, lDcpl, lSilent

nofg = GroupInfo%nofg; norg = GroupInfo%norg
ig1 = nofg + 1; ig2 = nofg + norg
nFxr = Core%nCoreFxr; nFsr = Core%nCoreFsr
myzb = PE%myzb; myze = PE%myze
nPolar = RayInfo%nPolarAngle; nPhiAngSv = RayInfo%nPhiAngSv
lSilent = .FALSE.
master = PE%master; RTmaster = PE%RTmaster

#ifdef MPI_ENV
lDcpl = nTracerCntl%lDcpl; lSilent = lDcpl
#endif

WRITE(mesg,'(a24)') 'Solving Subgroup FSP... '
IF(master) CALL message(io8, TRUE, TRUE, mesg)
Tbeg = nTracer_Dclock(FALSE, FALSE)

ALLOCATE(phis1g(nFsr)); ALLOCATE(phis1gd(nFsr))
ALLOCATE(xstr1g(nFsr)); ALLOCATE(src1g(nFsr))
ALLOCATE(PhiAngIn1g(nPolar, nPhiAngSv))
ALLOCATE(SigLamPot(nFxr))

write(mesg,'(a,f10.2,a)') "Reference Fuel Temperature",THInfo%RefFuelTemp(0)," C"
IF(master) call message(io8,TRUE,TRUE,mesg)
WRITE(mesg,'(a,i3)') 'Solving Subgroup FSP'
IF(master) CALL message(io8, TRUE, TRUE, mesg)

DO iz = myzb, myze
  IF(.NOT. Core%lFuelPlane(iz)) CYCLE
  WRITE(mesg,'(a12,i3,a30,i3)') 'SGFSP : Pln ',iz,'  # Of Resonance Categories : ',nActiveCat(iz)
  IF(RtMaster) CALL message(io8, TRUE, TRUE, mesg)
  rtnet = 0
  DO icat = 1, nCat
    if (.not.ResoCatUse(icat,iz)) CYCLE
    nid = ResoCat(icat)%repid;  id = mapnucl(nid)
    isodata => ldiso(id)
    nlvflx = isodata%nlv
    IF (nTracerCntl%l4Lv) nlvflx = isodata%nlvflx

    niteravg = 0
    DO ig = ig1, ig2

        errmax = 0._8; niter = 0
        DO ilv = 1, nlvflx
          IF (nTracerCntl%l4Lv) THEN
              lv = isodata%lvflx(ilv,ig)
          ELSE
              lv = isodata%lvabs(ilv,ig)
          ENDIF
          !Initialize the flux variable
          CALL CP_CA(phis1g, 1._8, nFsr); CALL CP_CA(PhiAngIn1g, 1._8, nPolar, nPhiAngSv)
          PhiAngIn1g(:, 1) = 0;
          !transport xs and source setup for Subgroup Fixed Source Problem
          IF(RTMASTER) THEN
            CALL SetPlnLsigP_cat(SigLamPot, xstr1g, lv, icat, Core, Fxr, ilv, iz, ig, PE)
            CALL SetSubGrpSrc1g(src1g, SigLamPot, xstr1g, Core, Fxr, iz, PE)
          ENDIF
          DO iter = 1, 100
            !Solving FSP
            IF(RtMaster) CALL CP_VA(phis1gd, phis1g, nfsr) !Save Previous data
            rt1 = nTracer_dclock(FALSE, FALSE)
            CALL RayTraceGM_OMP(RayInfo, Core, phis1g, PhiAngIn1g, xstr1g, src1g, jout, iz, .FALSE., nTracerCntl%FastMocLv)
            rt2 = nTracer_dclock(FALSE, FALSE)
            rtnet = rtnet + (rt2 - rt1)
            !Update Equivalence XS
            errmaxlv = SubGrpFspErr(phis1g, phis1gd, nfsr, PE)
            niter = niter + 1
            IF(RTMASTER) THEN
              CALL EquipXSGen(Phis1g, SigLamPot, xstr1g, ilv, icat, Core, Fxr, iz, ig, PE)
              CALL UpdtNDAF_CAT(ilv, icat, Core, Fxr, iz, ig, PE)
            ENDIF
            IF(errmaxlv .lt. epsm3) EXIT
            IF(RTMASTER) THEN
              CALL SetPlnLsigP_cat(SigLamPot, xstr1g, lv, icat, Core, Fxr, ilv, iz, ig, PE)
              CALL SetSubGrpSrc1g(src1g, SigLamPot, xstr1g, Core, Fxr, iz, PE)
            ENDIF
          ENDDO ! DO iter = 1, 100
          errmax = max(errmax,errmaxlv)
          TimeChk%SubGrpFSPNum = TimeChk%SubGrpFSPNum + 1
        ENDDO ! DO ilv = 1, nlvflx
        niteravg = niteravg + niter
        TimeChk%SubGrpRTniter = TimeChk%SubGrpRTniter + niter
    ENDDO ! DO ig = ig1, ig2
    niteravg = niteravg / norg
    WRITE(mesg,'(6x,a,i3,a,i6,a,f7.2,1p,e13.3)') 'Pln', iz,'  Istp.', nid, '  In_Itr_Avg', niteravg, errmax
    IF(master.and..NOT.lSilent) CALL MESSAGE(io8, TRUE, TRUE, mesg)

  ENDDO ! DO icat = 1, nCat
ENDDO ! DO iz = myzb, myze

DEALLOCATE(xstr1g, SigLamPot, phis1g, phis1gd, PhiAngIn1g)

nTracerCntl%lSubGrpSweep = TRUE
#ifdef MPI_ENV
CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
#endif
Tend = nTracer_Dclock(FALSE, FALSE)
TimeChk%SubGrpTime = TimeChk%SubGrpTime + (Tend - Tbeg)
TimeChk%NetRTSubGrpTime = TimeChk%NetRTSubGrpTime + rtnet

CONTINUE
END SUBROUTINE
SUBROUTINE SubGrpFsp_ISO(Core, Fxr, THInfo, RayInfo,  GroupInfo, nTracerCntl, PE)  !!!!!!!!!!****************************************!!!!!!!!!!
USE PARAM
USE TYPEDEF,        ONLY : coreinfo_type,    Fxrinfo_type,   GroupInfo_Type,    RayInfo_Type,   PE_TYPE, THInfo_Type
USE CNTL,           ONLY : nTracerCntl_Type
#ifdef MPI_ENV
USE MPICOMM_MOD,    ONLY : MPI_SYNC,    MPI_MAX_REAL,   MPI_MAX_INT,     BCAST
#endif
USE FILES,          ONLY : IO8
USE SUBGRP_MOD,     ONLY : SetPlnLsigP,    SubGrpFspErr,    EquipXsGen,    SetSubGrpSrc1g, UpdtNDAF
USE MOC_MOD,        ONLY : RayTraceGM_OMP
USE BasicOperation, ONLY : CP_CA,                CP_VA
USE IOUTIL,         ONLY : message
USE TIMER,          ONLY : nTracer_dclock,       TimeChk
USE XSLib_mod,      ONLY : LIBDATA, nRes_Cor, idRes_Cor, mapnucl, ldiso

IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type),POINTER :: Fxr(:, :)
TYPE(THInfo_Type) :: THInfo
TYPE(RayInfo_Type) :: RayInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE
TYPE(LIBDATA), POINTER :: isodata

INTEGER :: nfsr, nfxr, myzb, myze, iz, nlvflx, ilv, nid, id, irc
INTEGER :: ig, ig1, ig2, nofg, norg, niter, iter, nPhiAngSv, nPolar

REAL :: lv, niteravg, errmax, errmaxlv, Tbeg, Tend, rt1, rt2, rtnet
REAL, POINTER :: phis1g(:), phis1gd(:), PhiAngIn1g(:, :), jout(:,:,:)
REAL, POINTER :: siglamPot(:), xstr1g(:), src1g(:)

LOGICAL :: master, RTmaster, lDcpl, lSilent

nofg = GroupInfo%nofg; norg = GroupInfo%norg
ig1 = nofg + 1; ig2 = nofg + norg
nFxr = Core%nCoreFxr; nFsr = Core%nCoreFsr
myzb = PE%myzb; myze = PE%myze
nPolar = RayInfo%nPolarAngle; nPhiAngSv = RayInfo%nPhiAngSv
lSilent = .FALSE.
master = PE%master; RTmaster = PE%RTmaster

#ifdef MPI_ENV
lDcpl = nTracerCntl%lDcpl; lSilent = lDcpl
#endif

WRITE(mesg,'(a24)') 'Solving Subgroup FSP... '
IF(master) CALL message(io8, TRUE, TRUE, mesg)
Tbeg = nTracer_Dclock(FALSE, FALSE)

ALLOCATE(phis1g(nFsr)); ALLOCATE(phis1gd(nFsr))
ALLOCATE(xstr1g(nFsr)); ALLOCATE(src1g(nFsr))
ALLOCATE(PhiAngIn1g(nPolar, nPhiAngSv))
ALLOCATE(SigLamPot(nFxr))

write(mesg,'(a,f10.2,a)') "Reference Fuel Temperature",THInfo%RefFuelTemp(0)," C"
IF(master) call message(io8,TRUE,TRUE,mesg)
WRITE(mesg,'(a,i3)') 'Solving Subgroup FSP'
IF(master) CALL message(io8, TRUE, TRUE, mesg)

DO iz = myzb, myze
  IF(.NOT. Core%lFuelPlane(iz)) CYCLE

  WRITE(mesg,'(a12,i3,a28,i3)') 'SGFSP : Pln ',iz,'  # Of Resonance Isotopes : ',nRes_Cor(iz)
  IF(RTmaster) CALL message(io8,TRUE,TRUE,mesg)

  rtnet = 0
  DO irc = 1, nRes_Cor(iz)
    nid = idRes_Cor(irc,iz);  id = mapnucl(nid)
    if (id.eq.0) id = mapnucl(nid+500)
    isodata => ldiso(id)
    nlvflx = isodata%nlv
    IF (nTracerCntl%l4Lv) nlvflx = isodata%nlvflx

    niteravg = 0
    DO ig = ig1, ig2

        errmax = 0._8; niter = 0
        DO ilv = 1, nlvflx
          IF (nTracerCntl%l4Lv) THEN
              lv = isodata%lvflx(ilv,ig)
          ELSE
              lv = isodata%lvabs(ilv,ig)
          ENDIF
          !Initialize the flux variable
          CALL CP_CA(phis1g, 1._8, nFsr); CALL CP_CA(PhiAngIn1g, 1._8, nPolar, nPhiAngSv)
          PhiAngIn1g(:, 1) = 0;
          !transport xs and source setup for Subgroup Fixed Source Problem
          IF(RTMASTER) THEN
            CALL SetPlnLsigP(SigLamPot, xstr1g, lv, irc, Core, Fxr, ilv, iz, ig, PE)
            CALL SetSubGrpSrc1g(src1g, SigLamPot, xstr1g, Core, Fxr, iz, PE)
          ENDIF
          DO iter = 1, 100
            !Solving FSP
            IF(RtMaster) CALL CP_VA(phis1gd, phis1g, nfsr) !Save Previous data
            rt1 = nTracer_dclock(FALSE, FALSE)
            CALL RayTraceGM_OMP(RayInfo, Core, phis1g, PhiAngIn1g, xstr1g, src1g, jout, iz, .FALSE., nTracerCntl%FastMocLv)
            rt2 = nTracer_dclock(FALSE, FALSE)
            rtnet = rtnet + (rt2 - rt1)
            !Update Equivalence XS
            errmaxlv = SubGrpFspErr(phis1g, phis1gd, nfsr, PE)
            niter = niter + 1
            IF(RTMASTER) THEN
              CALL EquipXSGen(Phis1g, SigLamPot, xstr1g, ilv, irc, Core, Fxr, iz, ig, PE)
              CALL UpdtNDAF(id, ilv, irc, Core, Fxr, iz, ig, PE)
            ENDIF
            IF(errmaxlv .lt. epsm3) EXIT
            IF(RTMASTER) THEN
              CALL SetPlnLsigP(SigLamPot, xstr1g, lv, irc, Core, Fxr, ilv, iz, ig, PE)
              CALL SetSubGrpSrc1g(src1g, SigLamPot, xstr1g, Core, Fxr, iz, PE)
            ENDIF
          ENDDO ! DO iter = 1, 100
          errmax = max(errmax,errmaxlv)
          TimeChk%SubGrpFSPNum = TimeChk%SubGrpFSPNum + 1
        ENDDO ! DO ilv = 1, nlvflx
        niteravg = niteravg + niter
        TimeChk%SubGrpRTniter = TimeChk%SubGrpRTniter + niter

    ENDDO ! DO ig = ig1, ig2
    niteravg = niteravg / norg
    WRITE(mesg,'(6x,a,i3,a,i6,a,f7.2,1p,e13.3)') 'Pln', iz,'  Istp.', nid, '  In_Itr_Avg', niteravg, errmax
    IF(master.and..NOT.lSilent) CALL MESSAGE(io8, TRUE, TRUE, mesg)

  ENDDO ! DO irc = 1, nRes_Cor(iz)
ENDDO ! DO iz = myzb, myze

DEALLOCATE(xstr1g, SigLamPot, phis1g, phis1gd, PhiAngIn1g)

nTracerCntl%lSubGrpSweep = TRUE
#ifdef MPI_ENV
CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
#endif
Tend = nTracer_Dclock(FALSE, FALSE)
TimeChk%SubGrpTime = TimeChk%SubGrpTime + (Tend - Tbeg)
TimeChk%NetRTSubGrpTime = TimeChk%NetRTSubGrpTime + rtnet

CONTINUE
END SUBROUTINE
SUBROUTINE SubGrpFsp_MLG(Core, Fxr, THInfo, RayInfo,  GroupInfo, nTracerCntl, PE)  !!!!!!!!!!****************************************!!!!!!!!!!
USE PARAM
USE TYPEDEF,        ONLY : coreinfo_type,    Fxrinfo_type,   GroupInfo_Type,    RayInfo_Type,   PE_TYPE, THInfo_Type
USE CNTL,           ONLY : nTracerCntl_Type
#ifdef MPI_ENV
USE MPICOMM_MOD,    ONLY : MPI_SYNC,    MPI_MAX_REAL,   MPI_MAX_INT,     BCAST,     REDUCE
#endif
USE FILES,          ONLY : IO8
USE SUBGRP_MOD,     ONLY : SetPlnLsigP_MLG,  SetPlnLsigP_1gMLG,   SubGrpFspErr,  EquipXSGen_MLG,  &
                           EquipXsGen_1gMLG,  SetSubGrpSrc1g,  UpdtFnAdj, UpdtFtAdj
USE MOC_MOD,        ONLY : RayTraceGM_OMP
USE BasicOperation, ONLY : CP_CA,                CP_VA
USE IOUTIL,         ONLY : message
USE TIMER,          ONLY : nTracer_dclock,       TimeChk
USE OMP_LIB
USE XSLib_mod,      ONLY : mlgdata,mlgdata0
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type),POINTER :: Fxr(:, :)
TYPE(THInfo_Type) :: THInfo
TYPE(RayInfo_Type) :: RayInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE

INTEGER :: nfsr, nfxr, myzb, myze, iz, ilv, nlv, myitersum
INTEGER :: ig, ig1, ig2, nofg, norg, niter, iter, itermax, itersum, nPhiAngSv, nPolar

REAL :: lv, errmax, errmaxlv, Tbeg, Tend, rt1, rt2, rtnet
REAL, POINTER :: phis1g(:), phis1gd(:), PhiAngIn1g(:, :), jout(:,:,:)
REAL, POINTER :: siglamPot(:), xstr1g(:), src1g(:)

LOGICAL :: lCLD, lAIC, master, RTmaster, lDcpl, lSilent

nofg = GroupInfo%nofg; norg = GroupInfo%norg
ig1 = nofg + 1; ig2 = nofg + norg
nFxr = Core%nCoreFxr; nFsr = Core%nCoreFsr
myzb = PE%myzb; myze = PE%myze
nPolar = RayInfo%nPolarAngle; nPhiAngSv = RayInfo%nPhiAngSv
lSilent = .FALSE.
master = PE%master; RTmaster = PE%RTmaster

#ifdef MPI_ENV
lDcpl = nTracerCntl%lDcpl; lSilent = lDcpl
#endif

Tbeg = nTracer_Dclock(FALSE, FALSE)

itermax = 100
IF (any(Core%RadBC(1 : 4) .EQ. VoidCell)) itermax = 1

ALLOCATE(phis1g(nFsr)); ALLOCATE(phis1gd(nFsr))
ALLOCATE(xstr1g(nFsr)); ALLOCATE(src1g(nFsr))
ALLOCATE(PhiAngIn1g(nPolar, nPhiAngSv))
ALLOCATE(SigLamPot(nFxr))

write(mesg,'(a,f10.2,a)') "Reference Fuel Temperature",THInfo%RefFuelTemp(0)," C"
IF(master) call message(io8,TRUE,TRUE,mesg)
WRITE(mesg,'(a)') 'Solving Subgroup FSP (MLG)'
IF(master) CALL message(io8, TRUE, TRUE, mesg)

itersum = 0;
myitersum = 0;
DO iz = myzb, myze

  rtnet = 0
  nlv = mlgdata(iz)%f_nmaclv

  !Initialize the flux variable
  CALL CP_CA(phis1g, 1._8, nFsr); CALL CP_CA(PhiAngIn1g, 1._8, nPolar, nPhiAngSv)
  PhiAngIn1g(:, 1) = 0;

  IF(Core%lFuelPlane(iz)) THEN
  DO ig = ig1, ig2
      errmax = 0._8; niter = 0
      CALL UpdtFnAdj(Core, Fxr, ig, iz, PE)         ! NUMBER DENSITY CONSIDERATION FACTOR
      DO ilv = 1, nlv
        lv = mlgdata(iz)%f_maclv(ilv)
        !transport xs and source setup for Subgroup Fixed Source Problem
        IF(RTMASTER) THEN
          CALL SetPlnLsigP_MLG(SigLamPot, xstr1g, lv, Core, Fxr, ilv, iz, ig, PE)
          CALL SetSubGrpSrc1g(src1g, SigLamPot, xstr1g, Core, Fxr, iz, PE) ! SRC1G(ifsr) = SigLampot(ifxr) / xstr1g(ifsr)
        ENDIF
        DO iter = 1, itermax
          !Solving FSP
          IF(RtMaster) CALL CP_VA(phis1gd, phis1g, nfsr) !Save Previous data
          rt1 = nTracer_dclock(FALSE, FALSE)
          CALL RayTraceGM_OMP(RayInfo, Core, phis1g, PhiAngIn1g, xstr1g, src1g, jout, iz, .FALSE., nTracerCntl%FastMocLv)
          rt2 = nTracer_dclock(FALSE, FALSE)
          rtnet = rtnet + (rt2 - rt1)
          !Update Equivalence XS
          errmaxlv = SubGrpFspErr(phis1g, phis1gd, nfsr, PE)
          niter = niter + 1
          IF(RTMASTER) THEN
            CALL EquipXSGen_MLG(Phis1g, SigLamPot, xstr1g, ilv, Core, Fxr, iz, ig, PE) ! equivalence XS calculation with flux from FSP
            CALL UpdtFtAdj(Core, Fxr, ilv, ig, iz, PE) ! TCF update
          ENDIF
          IF(errmaxlv .lt. epsm3) EXIT
          IF(RTMASTER) THEN ! with updated TCF
            CALL SetPlnLsigP_MLG(SigLamPot, xstr1g, lv, Core, Fxr, ilv, iz, ig, PE)
            CALL SetSubGrpSrc1g(src1g, SigLamPot, xstr1g, Core, Fxr, iz, PE)
          ENDIF
        ENDDO ! DO iter = 1, 100
        errmax = max(errmax,errmaxlv)
        TimeChk%SubGrpFSPNum = TimeChk%SubGrpFSPNum + 1
      ENDDO ! DO ilv = 1, nlv
    WRITE(mesg,'(2x,2(a,i3),a,i5,1p,e13.3)') '[Fuel] Pln', iz, ' Group ',ig,'  In_itr', niter, errmax
    myitersum = myitersum + niter
    !TimeChk%SubGrpRTniter = TimeChk%SubGrpRTniter + niter
    IF(master .AND. .NOT. lSilent) CALL MESSAGE(io8, TRUE, TRUE, mesg)
  ENDDO ! ig
  ENDIF ! lfuelplane(iz)

  IF(Core%lcladPlane(iz)) THEN
    lCLD=.true.; lAIC=.false.
    nlv = mlgdata0%c_nmaclv1G

    errmax = 0._8; niter = 0
    DO ilv = 1, nlv
      lv = mlgdata0%c_maclv1G(ilv)
      !transport xs and source setup for Subgroup Fixed Source Problem
      IF(RTMASTER) THEN
        CALL SetPlnLsigP_1gMLG(SigLamPot, xstr1g, lv, Core, Fxr, iz, lCLD, lAIC, PE)
        CALL SetSubGrpSrc1g(src1g, SigLamPot, xstr1g, Core, Fxr, iz, PE)
      ENDIF
      DO iter = 1, itermax
        !Solving FSP
        IF(RtMaster) CALL CP_VA(phis1gd, phis1g, nfsr) !Save Previous data
        rt1 = nTracer_dclock(FALSE, FALSE)
        CALL RayTraceGM_OMP(RayInfo, Core, phis1g, PhiAngIn1g, xstr1g, src1g, jout, iz, .FALSE., nTracerCntl%FastMocLv)
        rt2 = nTracer_dclock(FALSE, FALSE)
        rtnet = rtnet + (rt2 - rt1)
        !Update Equivalence XS
        errmaxlv = SubGrpFspErr(phis1g, phis1gd, nfsr, PE)
        niter = niter + 1
        IF(errmaxlv .lt. epsm3) EXIT
      ENDDO ! DO iter = 1, 100
      IF(RTMASTER) CALL EquipXSGen_1gMLG(Phis1g, SigLamPot, xstr1g, ilv, Core, Fxr, iz, lCLD, lAIC, PE)
      errmax = max(errmax,errmaxlv)
      TimeChk%SubGrpFSPNum = TimeChk%SubGrpFSPNum + 1
    ENDDO ! DO ilv = 1, nlv
  WRITE(mesg,'(2x,a,i3,a,i5,1p,e13.3)') '[Clad] Pln', iz, '  In_itr', niter, errmax
  myitersum = myitersum + niter
  IF(master .AND. .NOT. lSilent) CALL MESSAGE(io8, TRUE, TRUE, mesg)
  ENDIF ! cladding

  IF (Core%lAICPlane(iz)) THEN
    nlv = mlgdata(iz)%f_nmaclv1G
    lCLD=.false.; lAIC=.true.
    errmax = 0._8; niter = 0
    DO ilv = 1, nlv
      lv = mlgdata(iz)%f_maclv1G(ilv)
      !transport xs and source setup for Subgroup Fixed Source Problem
      IF(RTMASTER) THEN
        CALL SetPlnLsigP_1gMLG(SigLamPot, xstr1g, lv, Core, Fxr, iz, lCLD, lAIC, PE)
        CALL SetSubGrpSrc1g(src1g, SigLamPot, xstr1g, Core, Fxr, iz, PE)
      ENDIF
      DO iter = 1, itermax
        !Solving FSP
        IF(RtMaster) CALL CP_VA(phis1gd, phis1g, nfsr) !Save Previous data
        rt1 = nTracer_dclock(FALSE, FALSE)
        CALL RayTraceGM_OMP(RayInfo, Core, phis1g, PhiAngIn1g, xstr1g, src1g, jout, iz, .FALSE., nTracerCntl%FastMocLv)
        rt2 = nTracer_dclock(FALSE, FALSE)
        rtnet = rtnet + (rt2 - rt1)
        !Update Equivalence XS
        errmaxlv = SubGrpFspErr(phis1g, phis1gd, nfsr, PE)
        niter = niter + 1
        IF(errmaxlv .lt. epsm3) EXIT
      ENDDO ! DO iter = 1, 100
      IF(RTMASTER) CALL EquipXSGen_1gMLG(Phis1g, SigLamPot, xstr1g, ilv, Core, Fxr, iz, lCLD, lAIC, PE)
      errmax = max(errmax,errmaxlv)
      TimeChk%SubGrpFSPNum = TimeChk%SubGrpFSPNum + 1
    ENDDO ! DO ilv = 1, nlv
  WRITE(mesg,'(2x,a,i3,a,i5,1p,e13.3)') '[AIC]  Pln', iz, '  In_itr', niter, errmax
  myitersum = myitersum + niter
  IF(master .AND. .NOT. lSilent) CALL MESSAGE(io8, TRUE, TRUE, mesg)
  ENDIF ! lAIC
ENDDO ! DO iz = myzb, myze

DEALLOCATE(xstr1g, SigLamPot, phis1g, phis1gd, PhiAngIn1g)

nTracerCntl%lSubGrpSweep = TRUE
itersum = myitersum;
#ifdef MPI_ENV
    CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
    CALL REDUCE(myitersum, itersum, PE%MPI_NTRACER_COMM, .FALSE.)
#endif
TimeChk%SubGrpRTniter = TimeChk%SubGrpRTniter + itersum
Tend = nTracer_Dclock(FALSE, FALSE)
TimeChk%SubGrpTime = TimeChk%SubGrpTime + (Tend - Tbeg)
TimeChk%NetRTSubGrpTime = TimeChk%NetRTSubGrpTime + rtnet

CONTINUE
END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CalcDancoff(Core, Fxr, RayInfo, THInfo, nTracerCntl, PE)  !!!!!!!!!!****************************************!!!!!!!!!!
USE PARAM
USE TYPEDEF,        ONLY : coreinfo_type,    Fxrinfo_type,   RayInfo_Type,  THInfo_Type,   PE_TYPE
USE CNTL,           ONLY : nTracerCntl_Type
#ifdef MPI_ENV
USE MPICOMM_MOD,    ONLY : MPI_SYNC,    MPI_MAX_REAL,   MPI_MAX_INT,     BCAST
#endif
USE FILES,          ONLY : IO8
USE SUBGRP_MOD,     ONLY : SetPlnLsigP_Dancoff, SetPlnLsigP_DancoffAIC, SetSubGrpSrc1g, Set_Dancoff, Set_DancoffAIC, SubGrpFspErr
USE MOC_MOD,        ONLY : RayTraceGM_OMP
USE BasicOperation, ONLY : CP_CA,                CP_VA
USE IOUTIL,         ONLY : message
USE TIMER,          ONLY : nTracer_dclock,       TimeChk
USE OMP_LIB
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type),POINTER :: Fxr(:, :)
TYPE(RayInfo_Type) :: RayInfo
TYPE(THInfo_Type) :: THInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE

INTEGER :: nfsr, nfxr, myzb, myze, iz, ilv, nlv, niterCP
INTEGER :: niter, iter, nPhiAngSv, nPolar

REAL :: errmax, Tbeg, Tend
REAL, POINTER :: phis1g(:), phis1gd(:), PhiAngIn1g(:, :), jout(:,:,:)
REAL, POINTER :: siglamPot(:), xstr1g(:), src1g(:)

LOGICAL :: master, RTmaster, lDcpl, lSilent

nFxr = Core%nCoreFxr; nFsr = Core%nCoreFsr
myzb = PE%myzb; myze = PE%myze
nPolar = RayInfo%nPolarAngle; nPhiAngSv = RayInfo%nPhiAngSv
lSilent = .FALSE.
master = PE%master; RTmaster = PE%RTmaster

#ifdef MPI_ENV
lDcpl = nTracerCntl%lDcpl; lSilent = lDcpl
#endif

Tbeg = nTracer_Dclock(FALSE, FALSE)

ALLOCATE(phis1g(nFsr)); ALLOCATE(phis1gd(nFsr))
ALLOCATE(xstr1g(nFsr)); ALLOCATE(src1g(nFsr))
ALLOCATE(PhiAngIn1g(nPolar, nPhiAngSv))
ALLOCATE(SigLamPot(nFxr))

write(mesg,'(a,f10.2,a)') "Reference Fuel Temperature",THInfo%RefFuelTemp(0)," C"
IF(master) call message(io8,TRUE,TRUE,mesg)
WRITE(mesg,'(a,i3)') 'Solving FSP for Pin-wise Dancoff Factor'
IF(master) CALL message(io8, TRUE, TRUE, mesg)

DO iz = myzb, myze
  IF(.NOT. Core%lFuelPlane(iz)) CYCLE

  !Initialize the flux variable
  CALL CP_CA(phis1g, 1._8, nFsr); CALL CP_CA(PhiAngIn1g, 1._8, nPolar, nPhiAngSv)
  PhiAngIn1g(:, 1) = 0;

  ! Enhanced Neutron Current Method to Calculate Dancoff Factor
  errmax = 0._8; niter = 0
  !transport xs and source setup for Subgroup Fixed Source Problem
  IF(RTMASTER) THEN
    CALL SetPlnLsigP_Dancoff(SigLamPot, xstr1g, Core, Fxr, iz, PE) ! set plane-wise XST, and XSP for Enhanced Neutron Current Method
    CALL SetSubGrpSrc1g(src1g, SigLamPot, xstr1g, Core, Fxr, iz, PE)
  ENDIF
  DO iter = 1, 100
    !Solving ENCM
    IF(RtMaster) CALL CP_VA(phis1gd, phis1g, nfsr) !Save Previous data
    CALL RayTraceGM_OMP(RayInfo, Core, phis1g, PhiAngIn1g, xstr1g, src1g, jout, iz, .FALSE., nTracerCntl%FastMocLv)
    !Update Equivalence XS
    errmax = SubGrpFspErr(phis1g, phis1gd, nfsr, PE)
    niter = niter + 1
    IF(errmax .lt. epsm3) EXIT
  ENDDO ! DO iter = 1, 100

  IF(RTMASTER) then
      CALL Set_Dancoff(SigLamPot, xstr1g, phis1g, Core, Fxr, iz, niterCP, PE)
      nTracerCntl%nCP_er = nTracerCntl%nCP_er + niterCP
  ENDIF

  WRITE(mesg,'(2x,a,i3,a,i5,1p,e13.3,a,i6)') '[Fuel] Pln', iz,'  In_itr',niter,errmax,'   # of CP : ',niterCP
  TimeChk%SubGrpRTniter = TimeChk%SubGrpRTniter + niter
  IF(master .AND. .NOT. lSilent) CALL MESSAGE(io8, TRUE, TRUE, mesg)

  IF (Core%lAICPlane(iz)) THEN

    !Initialize the flux variable
    CALL CP_CA(phis1g, 1._8, nFsr); CALL CP_CA(PhiAngIn1g, 1._8, nPolar, nPhiAngSv)
    PhiAngIn1g(:, 1) = 0;

    errmax = 0._8; niter = 0
    !transport xs and source setup for Subgroup Fixed Source Problem
    IF(RTMASTER) THEN
      CALL SetPlnLsigP_DancoffAIC(SigLamPot, xstr1g, Core, Fxr, iz, PE)
      CALL SetSubGrpSrc1g(src1g, SigLamPot, xstr1g, Core, Fxr, iz, PE)
    ENDIF
    DO iter = 1, 100
      !Solving FSP
      IF(RtMaster) CALL CP_VA(phis1gd, phis1g, nfsr) !Save Previous data
      CALL RayTraceGM_OMP(RayInfo, Core, phis1g, PhiAngIn1g, xstr1g, src1g, jout, iz, .FALSE., nTracerCntl%FastMocLv)
      !Update Equivalence XS
      errmax = SubGrpFspErr(phis1g, phis1gd, nfsr, PE)
      niter = niter + 1
      IF(errmax .lt. epsm3) EXIT
    ENDDO ! DO iter = 1, 100

    IF(RTMASTER) THEN
        CALL Set_DancoffAIC(SigLamPot, xstr1g, phis1g, Core, Fxr, iz, niterCP, PE)
        nTracerCntl%nCP_er = nTracerCntl%nCP_er + niterCP
    ENDIF

    WRITE(mesg,'(2x,a,i3,a,i5,1p,e13.3,a,i6)') '[AIC]  Pln', iz,'  In_itr',niter,errmax,'   # of CP : ',niterCP
    TimeChk%SubGrpRTniter = TimeChk%SubGrpRTniter + niter
    IF(master .AND. .NOT. lSilent) CALL MESSAGE(io8, TRUE, TRUE, mesg)

  ENDIF

ENDDO ! DO iz = myzb, myze

DEALLOCATE(xstr1g, SigLamPot, phis1g, phis1gd, PhiAngIn1g)

nTracerCntl%lSubGrpSweep = TRUE
#ifdef MPI_ENV
CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
#endif
Tend = nTracer_Dclock(FALSE, FALSE)
TimeChk%DancoffTime = TimeChk%DancoffTime + (Tend - Tbeg)

CONTINUE
END SUBROUTINE
SUBROUTINE CalcEscXSCP(Core, Fxr, GroupInfo, THInfo, nTracerCntl, PE)  !!!!!!!!!!****************************************!!!!!!!!!!
USE PARAM
USE TYPEDEF,        ONLY : coreinfo_type,    Fxrinfo_type, ResVarPin_Type,  Cell_TYPE, Pin_TYPE,  GroupInfo_Type, PE_TYPE, THInfo_Type
USE CNTL,           ONLY : nTracerCntl_Type
USE FILES,          ONLY : IO8
USE BasicOperation, ONLY : CP_CA,                CP_VA
USE IOUTIL,         ONLY : message
USE TIMER,          ONLY : nTracer_dclock,       TimeChk
USE OMP_LIB
USE XSLib_mod,      ONLY : mlgdata, mlgdata0
USE CP_mod
USE TH_Mod,         ONLY : GetPinFuelTemp
USE SUBGRP_MOD,     ONLY : SubGrpFspErr, GetSigtrlp_pin, UpdtFnAdj_pin, UpdtFtAdj_pin
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

integer :: nofg, norg, ig1, ig2, myzb, myze, nxy, iz, nitertot, nintertot_
integer :: nlvf, nlvc, ilv, FxrIdxSt, ipin, igt, ig, icel, nlocalFxr
integer :: iter, k, j, ifxr
integer,pointer,dimension(:) :: ridx,idx
real :: Tbeg, Tend, TempAvgsq, errmax, lv
real :: lvsum, lpsum, phisum, areasum, area
real,pointer,dimension(:) :: lvf,lvc,LocalSigtr,LocalSiglpC,LocalSigtrC,phi,phid,FtAdj,lvadj
real,pointer,dimension(:,:) :: FnAdj,LocalSiglp

logical :: master

nofg = GroupInfo%nofg; norg = GroupInfo%norg
ig1 = nofg + 1; ig2 = nofg + norg
myzb = PE%myzb; myze = PE%myze
master = PE%master

Tbeg = nTracer_Dclock(FALSE, FALSE)
WRITE(mesg,'(a,i3)') 'Solving Pin-wise Escape XS with 1D CP'
IF(master) CALL message(io8, TRUE, TRUE, mesg)

nxy = Core%nxy
Pin => Core%Pin; CellInfo => Core%CellInfo
ResVarPin => Core%ResVarPin

nintertot_=0
DO iz = myzb, myze
  IF(.NOT. Core%lFuelPlane(iz)) CYCLE

  nlvf = mlgdata(iz)%f_nmaclv
  nlvc = mlgdata0%c_nmaclv1G
  allocate(lvc(nlvc))
  do ilv = 1, nlvc
    lvc(ilv) = mlgdata0%c_maclv1G(ilv)
  enddo

  nitertot=0
  !$OMP PARALLEL DEFAULT(SHARED)      &
  !$OMP PRIVATE(ipin,FxrIdxSt,TempAvgsq,RP,igt,ilv,lvf,icel,nlocalFxr,FnAdj,ridx,idx,LocalSiglp,LocalSigtr, &
  !$OMP         LocalSiglpC,LocalSigtrC,phi,phid,FtAdj,lvadj,ig,iter,j,k,ifxr,myFxr,errmax,lvsum,lpsum,phisum, &
  !$OMP         areasum,area,lv)
  !$OMP DO REDUCTION(+:nitertot)
  DO ipin = 1, nxy
    if (.not.ResVarPin(ipin,iz)%lres) cycle
    FxrIdxSt = Pin(ipin)%FxrIdxSt
    TempAvgsq = dsqrt(GetPinFuelTemp(Core, Fxr, iz, ipin))
    RP => ResVarPin(ipin,iz)
    igt = RP%igt
    allocate(lvf(nlvf))
    do ilv = 1, nlvf
      lvf(ilv) = mlgdata(iz)%f_maclv_pin(ilv,igt)
    enddo

    icel = Pin(ipin)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
    allocate(FnAdj(nlocalFxr,ig1:ig2),ridx(nlocalFxr),idx(nlocalFxr))
    allocate(LocalSiglp(nlocalFxr,ig1:ig2),LocalSigtr(nlocalFxr))
    allocate(LocalSiglpC(nlocalFxr),LocalSigtrC(nlocalFxr))
    allocate(phi(nlocalFxr),phid(nlocalFxr),FtAdj(nlocalFxr),lvadj(nlocalFxr))
    call GetSigtrlp_pin(LocalSiglp,LocalSiglpC,Fxr,ig1,ig2,ridx,idx,FxrIdxSt,nlocalFxr,iz)

    if (.not.CellInfo(icel)%lAIC) then
      ! fuel
      call UpdtFnAdj_pin(FnAdj,Fxr,TempAvgsq,ig1,ig2,ridx,idx,nlocalFxr,iz)
      CALL CP_CA(phi, 1._8, nlocalFxr)
      do ig = ig1, ig2
        do ilv = 1, nlvf
          DO iter = 1, 100
            CALL CP_VA(phid, phi, nlocalFxr) !Save Previous data
            call UpdtFtAdj_pin(FtAdj,LocalSiglp(:,ig),Fxr,ilv,ig,ridx,idx,FxrIdxSt,nlocalFxr,TempAvgsq,iz)
            do j = 1, nLocalFxr
              k = ridx(j)
              ifxr = idx(j)
              myFxr => Fxr(ifxr, iz)
              !LocalSigtr(k) = LocalSiglp(k,ig)
              LocalSigtr(k) = LocalSiglpC(k)
              if (.not.myFxr%lres) cycle
              if (myFxr%lCLD) cycle
              lvadj(k) = FtAdj(k)*FnAdj(k,ig)*lvf(ilv)
              LocalSigtr(k) = LocalSigtr(k) + lvadj(k)
            enddo
            !call runCP(phi,RP%delr,RP%vol,RP%Qsurfvol,RP%invsurf4,RP%X,LocalSigtr,LocalSiglp(:,ig),nLocalFxr)
            call runCP(phi,RP%delr,RP%vol,RP%Qsurfvol,RP%invsurf4,RP%X,LocalSigtr,LocalSiglpC,nLocalFxr)
            do j = 1, nLocalFxr
              k = ridx(j)
              ifxr = idx(j)
              myFxr => Fxr(ifxr, iz)
              if (.not.myFxr%lres) cycle
              if (myFxr%lCLD) cycle
              !myFxr%XsEq_f_mg(ilv,ig) = lvadj(k) * phi(k) / (1._8 - phi(k)) - LocalSiglp(k,ig)
              myFxr%XsEq_f_mg(ilv,ig) = lvadj(k) * phi(k) / (1._8 - phi(k)) - LocalSiglpC(k)
            enddo
            errmax = SubGrpFspErr(phi, phid, nlocalFxr, PE)
            IF(errmax .lt. epsm3) EXIT
          ENDDO ! DO iter = 1, 100
          nitertot = nitertot + iter
          lvsum = 0._8; lpsum = 0._8; phisum = 0._8; areasum = 0._8
          do j = 1, nLocalFxr
            k = ridx(j)
            ifxr = idx(j)
            myFxr => Fxr(ifxr, iz)
            if (.not.myFxr%lres) cycle
            if (myFxr%lCLD) cycle
            area = myFxr%area
            lvsum = lvsum + lvadj(k) * phi(k) * area
            !lpsum = lpsum + LocalSiglp(k,ig) * phi(k) * area
            lpsum = lpsum + LocalSiglpC(k) * phi(k) * area
            phisum = phisum + phi(k) * area
            areasum = areasum + area
          enddo
          lvsum = lvsum / phisum
          lpsum = lpsum / phisum
          phisum = phisum / areasum
          RP%avgxseq_mg(ilv,ig) = lvsum * phisum / (1._8 - phisum) - lpsum
        enddo ! do ilv = 1, nlvf
      enddo ! do ig = ig1, ig2

    else

      ! AIC
      do ilv = 1, nlvf
        lv = lvf(ilv)
        do j = 1, nLocalFxr
          k = ridx(j)
          ifxr = idx(j)
          myFxr => Fxr(ifxr, iz)
          LocalSigtr(k) = LocalSiglpC(k)
          if (.not.myFxr%lres) cycle
          if (myFxr%lCLD) cycle
          LocalSigtr(k) = LocalSigtr(k) + lv
        enddo
        call runCP(phi,RP%delr,RP%vol,RP%Qsurfvol,RP%invsurf4,RP%X,LocalSigtr,LocalSiglpC,nLocalFxr)
        nitertot = nitertot + 1
        do j = 1, nLocalFxr
          k = ridx(j)
          ifxr = idx(j)
          myFxr => Fxr(ifxr, iz)
          if (.not.myFxr%lres) cycle
          if (myFxr%lCLD) cycle
          myFxr%XsEq_f_1g(ilv) = lv * phi(k) / (1._8 - phi(k)) - LocalSiglpC(k)
        enddo

        lpsum = 0._8; phisum = 0._8; areasum = 0._8
        do j = 1, nLocalFxr
          k = ridx(j)
          ifxr = idx(j)
          myFxr => Fxr(ifxr, iz)
          if (.not.myFxr%lres) cycle
          if (myFxr%lCLD) cycle
          area = myFxr%area
          lpsum = lpsum + LocalSiglpC(k) * phi(k) * area
          phisum = phisum + phi(k) * area
          areasum = areasum + area
        enddo
        lpsum = lpsum / phisum
        phisum = phisum / areasum
        RP%avgxseq_1g(ilv) = lv * phisum / (1._8 - phisum) - lpsum
      enddo ! do ilv = 1, nlvf

    endif

    ! Clad
    do ilv = 1, nlvc
      lv = lvc(ilv)
      do j = 1, nLocalFxr
        k = ridx(j)
        ifxr = idx(j)
        myFxr => Fxr(ifxr, iz)
        LocalSigtrC(k) = LocalSiglpC(k)
        if (.not.myFxr%lres) cycle
        if (.not.myFxr%lCLD) cycle
        LocalSigtrC(k) = LocalSigtrC(k) + lv
      enddo
      call runCP(phi,RP%delr,RP%vol,RP%Qsurfvol,RP%invsurf4,RP%X,LocalSigtrC,LocalSiglpC,nLocalFxr)
      nitertot = nitertot + 1
      do j = 1, nLocalFxr
        k = ridx(j)
        ifxr = idx(j)
        myFxr => Fxr(ifxr, iz)
        if (.not.myFxr%lres) cycle
        if (.not.myFxr%lCLD) cycle
        myFxr%XsEq_C_1g(ilv) = lv * phi(k) / (1._8 - phi(k)) - LocalSiglpC(k)
      enddo
    enddo ! do ilv = 1, nlvc

    deallocate(lvf,FnAdj,ridx,idx,LocalSiglp,LocalSigtr)
    deallocate(LocalSiglpC,LocalSigtrC,phi,phid,FtAdj,lvadj)

  ENDDO ! DO ipin = 1, nxy
  !$OMP END DO
  !$OMP END PARALLEL
  nintertot_ = nintertot_ + nitertot
  deallocate(lvc)

ENDDO ! DO iz = myzb, myze
#ifdef MPI_ENV
CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
#endif
#ifdef MPI_ENV
CALL REDUCE(nintertot_, nitertot, PE%MPI_NTRACER_COMM, .TRUE.)
#endif
nTracerCntl%nCP = nTracerCntl%nCP + nitertot
WRITE(mesg,'(a,i8,a)') 'Finished the Calculation of Pin-wise Escape XSs with ',nitertot,' 1D CP Calc.'
IF(master) CALL MESSAGE(io8, TRUE, TRUE, mesg)

nTracerCntl%lSubGrpSweep = TRUE
Tend = nTracer_Dclock(FALSE, FALSE)
TimeChk%CPTime = TimeChk%CPTime + (Tend - Tbeg)

CONTINUE
END SUBROUTINE
SUBROUTINE UpdtFnAdj_pin(FnAdj,Fxr,TempAvgsq,ig1,ig2,ridx,idx,nlocalFxr,iz)
! NUMBER DENSITY Correction factor for EDC model (NDCF compared to pin RI average)
USE TYPEDEF,        ONLY : Fxrinfo_type
USE XSLib_mod,      ONLY : libdata, mapnucl, ldiso
USE XsUtil_mod,     ONLY : LineIntPol
IMPLICIT NONE
TYPE(FxrInfo_Type),POINTER :: Fxr(:, :),myfxr
INTEGER :: nlocalFxr
INTEGER :: ridx(nlocalFxr),idx(nlocalFxr),ig1,ig2,iz
INTEGER :: j,k,ifxr,niso,iso,id,npot,nRiTemp,ig
REAL :: FnAdj(nlocalFxr,ig1:ig2),var(200),ria
REAL :: Npin(ig1:ig2),N(ig1:ig2),pinareasum,TempAvgsq
INTEGER, POINTER :: idiso(:)
REAL, POINTER :: pnum(:)
TYPE(libdata), POINTER :: isodata

Npin(ig1:ig2)=0._8; pinareasum=0._8
DO j = 1, nLocalFxr
  k = ridx(j)
  FnAdj(k,ig1:ig2) = 0._8
  ifxr = idx(j)
  myFxr => Fxr(ifxr, iz)
  if (.not.myFxr%lres) cycle
  if (myFxr%lCLD) cycle
  niso = myFxr%niso; idiso => myFxr%idiso; pnum => myFxr%pnum
  N(ig1:ig2) = 0._8
  DO iso = 1, niso
    id = mapnucl(idiso(iso)); isodata => ldiso(id)
    if (.not.isodata%lreso) cycle
    npot = isodata%nsig0; nRiTemp = isodata%nrtemp
    do ig = ig1, ig2
      var(1:nRiTemp) = isodata%ri_a(npot,ig, 1:nRiTemp)
      ria = LineIntpol(TempAvgsq, nRiTemp, isodata%rtempsq(1:nRiTemp), var(1:nRiTemp))
      N(ig) = N(ig) + pnum(iso) * ria
    enddo
  ENDDO
  pinareasum = pinareasum + myFxr%area
  do ig = ig1, ig2
    Npin(ig) = Npin(ig) + N(ig) * myFxr%area
    FnAdj(k,ig) = N(ig)
  enddo
ENDDO ! DO j = 1, nLocalFxr
!do ig = ig1, ig2
!  Npin(ig) = Npin(ig)/pinareasum
!enddo
Npin(ig1:ig2) = Npin(ig1:ig2)/pinareasum
DO j = 1, nLocalFxr
  k = ridx(j)
  ifxr = idx(j)
  myFxr => Fxr(ifxr, iz)
  if (.not.myFxr%lres) cycle
  if (myFxr%lCLD) cycle
  do ig = ig1, ig2
    FnAdj(k,ig) = FnAdj(k,ig)/Npin(ig)
    !FnAdj(k,ig) = 1.
  enddo
ENDDO
END SUBROUTINE
SUBROUTINE GetSigtrlp_pin(LocalSiglp,LocalSiglpC,Fxr,ig1,ig2,ridx,idx,FxrIdxSt,nlocalFxr,iz)
USE TYPEDEF,        ONLY : Fxrinfo_type
USE XSLib_mod,      ONLY : libdata, mapnucl, ldiso
IMPLICIT NONE
TYPE(FxrInfo_Type),POINTER :: Fxr(:, :),myfxr
INTEGER :: nlocalFxr
INTEGER :: ridx(nlocalFxr),idx(nlocalFxr),ig1,ig2,iz
INTEGER :: k,j,ifxr,niso,iso,id,ig,FxrIdxSt
REAL :: LocalSiglp(nlocalFxr,ig1:ig2)
REAL :: LocalSiglpC(nlocalFxr)
INTEGER, POINTER :: idiso(:)
REAL, POINTER :: pnum(:)
TYPE(libdata), POINTER :: isodata
DO j = 1, nLocalFxr
  k = nLocalFxr - j + 1
  ridx(j) = k
  ifxr = FxrIdxSt + j - 1
  idx(j) = ifxr
  myFxr => Fxr(ifxr, iz)
  niso = myFxr%niso; pnum => myFxr%pnum; idiso => myFxr%idiso
  LocalSiglp(k,ig1:ig2)=0._8; LocalSiglpC(k)=0._8
  DO iso = 1, niso
    id = mapnucl(idiso(iso));   isodata => ldiso(id)
    DO ig = ig1, ig2
      LocalSiglp(k,ig) = LocalSiglp(k,ig) + pnum(iso) * isodata%lamsigp(ig)
    ENDDO
    LocalSiglpC(k) = LocalSiglpC(k) + pnum(iso) * isodata%sigp
  ENDDO
ENDDO
END SUBROUTINE
SUBROUTINE UpdtFtAdj_pin(FtAdj,LocalSiglp,Fxr,ilv,ig,ridx,idx,FxrIdxSt,nlocalFxr,TempAvgsq,iz)
USE TYPEDEF,        ONLY : Fxrinfo_type
USE XSLib_mod,      ONLY : libdata, mapnucl, mapnuclres, ldiso
USE XsUtil_mod,     ONLY : LineIntPol, LineIntPol2
IMPLICIT NONE
TYPE(FxrInfo_Type),POINTER :: Fxr(:, :),myfxr
INTEGER :: nlocalFxr
INTEGER :: ridx(nlocalFxr),idx(nlocalFxr),ilv,FxrIdxSt,iz
INTEGER :: i,j,k,ifxr,niso,iso,id,npot,nRiTemp,ig,idres,jso,jd,jdres
REAL :: LocalSiglp(nlocalFxr),FtAdj(nlocalFxr),var(200),ria,xdat(20),ydat(20)
REAL :: Tempsq,Navg,Nreg,ind,micsigb,micsig0,I_avg,I_reg,sigbsq,TempAvgsq
INTEGER, POINTER :: idiso(:)
REAL, POINTER :: pnum(:)
TYPE(libdata), POINTER :: isodata, jsodata
DO j = 1, nLocalFxr
  k = ridx(j)
  ifxr = idx(j)
  myFxr => Fxr(ifxr, iz)
  FtAdj(k) = 1.
  if (.not.myFxr%lres) cycle
  if (myFxr%lCLD) cycle
  if (myFxr%xseq_f_mg(ilv,ig).eq.0._8) cycle
  Tempsq = dsqrt(myFxr%temp)
  niso = myFxr%niso; idiso => myFxr%idiso; pnum => myFxr%pnum
  Nreg = 0._8; Navg = 0._8
  DO iso = 1, niso
    id = mapnucl(idiso(iso));   isodata => ldiso(id)
    if (.NOT.isodata%lreso) CYCLE
    idres = mapnuclres(idiso(iso),iz)
    ind = 0._8
    DO jso = 1, niso
      jd = mapnucl(idiso(jso));   jsodata => ldiso(jd)
      jdres = mapnuclres(idiso(jso),iz)
      if (jdres.eq.idres) ind = ind + pnum(jso)
    ENDDO
    micsigb = (myFxr%XsEq_f_mg(ilv,ig) + LocalSiglp(k)) / ind / myFxr%FtAdj(ilv,ig)
    micsig0 = micsigb - isodata%lamsigp(ig)
    if (micsig0.le.0._8) micsig0 = 1E-10
    sigbsq = dsqrt(micsig0)
    npot = isodata%nsig0
    nRiTemp = isodata%nrtemp
    DO i = 1, nRiTemp
      xdat(i) = isodata%rtempsq(i)
      var(1:npot) = isodata%ri_a(1:npot, ig, i)
      ydat(i) = LineIntPol2(sigbsq ,npot, isodata%sig0sq(1:npot), var(1:npot))
    ENDDO
    I_reg = LineIntpol(tempsq, nRiTemp, xdat(1:nRiTemp), ydat(1:nRiTemp))
    Nreg = Nreg + pnum(iso) * I_reg
    I_avg = LineIntpol(TempAvgsq, nRiTemp, xdat(1:nRiTemp), ydat(1:nRiTemp))
    Navg = Navg + pnum(iso) * I_avg
  ENDDO ! DO iso = 1, niso
  FtAdj(k) = Nreg / Navg
  myFxr%FtAdj(ilv,ig) = FtAdj(k)
ENDDO ! DO j = 1, nLocalFxr
END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SetPlnLsigP(Siglp, Sigtr, lv, irc, Core, Fxr, ilv, iz, ig, PE)
USE TYPEDEF,     ONLY : CoreInfo_Type, FxrInfo_Type, pin_type, cell_type, PE_Type
USE XsLib_Mod,   ONLY : libdata, ldiso, mapnucl, mapnuclres
IMPLICIT NONE
TYPE(PE_Type) :: PE
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type),POINTER :: Fxr(:, :), myFxr
TYPE(Pin_TYPE), POINTER :: PIN(:)
TYPE(Cell_TYPE), POINTER :: CellInfo(:)
TYPE(libdata), POINTER :: isodata

INTEGER,INTENT(IN) :: irc, ilv, iz, ig
REAL,INTENT(IN) :: lv
REAL,POINTER :: Siglp(:), Sigtr(:)

INTEGER :: FxrIdxSt, FsrIdxSt, icel, nlocalFxr, nxy, nFsrInFxr, ipin, ifxr, ifsr, j, i
INTEGER :: niso, idres, id, iso, nThread
INTEGER, POINTER :: idiso(:)
REAL :: Localsiglp, LocalSigtr
REAL, POINTER :: pnum(:)

nxy = Core%nxy
Pin => Core%Pin; CellInfo => Core%CellInfo

nThread = PE%nThread
CALL OMP_SET_NUM_THREADS(nThread)

!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(ipin, FxrIdxSt, FsrIdxSt, icel, nlocalFxr, j, ifxr, myFxr, niso, pnum, idiso, LocalSiglp, LocalSigtr, &
!$OMP         iso, id, isodata, idres, nFsrInFxr, i, ifsr )
!$OMP DO
DO ipin = 1, nxy
  FxrIdxSt = Pin(ipin)%FxrIdxSt; FsrIdxSt = Pin(ipin)%FsrIdxSt
  icel = Pin(ipin)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j -1
    myFxr => Fxr(ifxr, iz)
    niso = myFxr%niso; pnum => myFxr%pnum; idiso => myFxr%idiso
    LocalSiglp = 0._8; LocalSigtr = 0._8
    if (myFxr%lres) then
      DO iso = 1, niso
        id = mapnucl(idiso(iso));   isodata => ldiso(id)
        idres = mapnuclres(idiso(iso),iz)
          if (isodata%sigp == 0.d0) then
        LocalSiglp = LocalSiglp + pnum(iso) * isodata%lamsigp(ig)
          else
            LocalSiglp = LocalSiglp + pnum(iso) * isodata%sigp
          end if
        if (idres.eq.irc) LocalSigtr = LocalSigtr + myFxr%NDAF(ilv,irc,ig) * pnum(iso) * lv
      ENDDO
    LocalSigtr = LocalSigtr + LocalSiglp
    else
      DO iso = 1, niso
        id = mapnucl(idiso(iso));   isodata => ldiso(id)
        idres = mapnuclres(idiso(iso),iz)
        if (isodata%sigp == 0.d0) then
          LocalSiglp = LocalSiglp + pnum(iso) * isodata%lamsigp(ig)
        else
          LocalSiglp = LocalSiglp + pnum(iso) * isodata%sigp
        end if
      ENDDO
      LocalSigtr = LocalSiglp
    endif

    Siglp(ifxr) = LocalSiglp
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      Sigtr(ifsr) = LocalSigtr
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL
NULLIFY(myFxr,pnum,idiso,isodata,PIN,CellInfo)
END SUBROUTINE
SUBROUTINE SetPlnLsigP_CAT(Siglp, Sigtr, lv, icat, Core, Fxr, ilv, iz, ig, PE)
USE TYPEDEF,     ONLY : CoreInfo_Type, FxrInfo_Type, pin_type, cell_type, PE_Type
USE XsLib_Mod,   ONLY : libdata, ldiso, mapnucl, ResoCat
USE XsUtil_mod,  ONLY : LineIntPol
USE TH_Mod,      ONLY : GetPinFuelTemp
IMPLICIT NONE
TYPE(PE_Type) :: PE
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type),POINTER :: Fxr(:, :), myFxr
TYPE(Pin_TYPE), POINTER :: PIN(:)
TYPE(Cell_TYPE), POINTER :: CellInfo(:)
TYPE(libdata), POINTER :: isodata, repisodata

INTEGER,INTENT(IN) :: icat, ilv, iz, ig
REAL,INTENT(IN) :: lv
REAL,POINTER  :: Siglp(:), Sigtr(:)

INTEGER :: FxrIdxSt, FsrIdxSt, icel, nlocalFxr, nxy, nFsrInFxr, ipin, ifxr, ifsr
INTEGER :: niso, j, i, id, iso, repnid, repid, npot, nRiTemp, nThread
INTEGER, POINTER :: idiso(:)
REAL :: tempavgsq, Localsiglp, LocalSigtr, ria, repria, effnum, var(200)
REAL, POINTER :: pnum(:)

nxy = Core%nxy
Pin => Core%Pin; CellInfo => Core%CellInfo

repnid = ResoCat(icat)%repid; repid = mapnucl(repnid)
repisodata => ldiso(repid)

nThread = PE%nThread
CALL OMP_SET_NUM_THREADS(nThread)

!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(ipin, FxrIdxSt, FsrIdxSt, icel, nlocalFxr, tempavgsq, j, ifxr, myFxr, niso, pnum, idiso, LocalSiglp, LocalSigtr, &
!$OMP         npot, nRiTemp, repria, effnum, iso, id, isodata, ria, nFsrInFxr, i, ifsr, var )
!$OMP DO
DO ipin = 1, nxy
  FxrIdxSt = Pin(ipin)%FxrIdxSt; FsrIdxSt = Pin(ipin)%FsrIdxSt
  icel = Pin(ipin)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
  tempavgsq = dsqrt(GetPinFuelTemp(Core, Fxr, iz, ipin))
  npot = repisodata%nsig0; nRiTemp = repisodata%nrtemp
  var(1:nRiTemp) = repisodata%ri_a(npot,ig, 1:nRiTemp)
  repria = LineIntpol(tempavgsq, nRiTemp, repisodata%rtempsq(1:nRiTemp), var(1:nRiTemp))
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j -1
    myFxr => Fxr(ifxr, iz)
    niso = myFxr%niso; pnum => myFxr%pnum; idiso => myFxr%idiso
    LocalSiglp = 0._8
    DO iso = 1, niso
      id = mapnucl(idiso(iso));   isodata => ldiso(id)
      !LocalSiglp = LocalSiglp + pnum(iso) * isodata%lamsigp(ig)
      if (isodata%sigp == 0.d0) then
      LocalSiglp = LocalSiglp + pnum(iso) * isodata%lamsigp(ig)
      else
        LocalSiglp = LocalSiglp + pnum(iso) * isodata%sigp
      end if
    ENDDO
    LocalSigtr = LocalSiglp

    if (myFxr%lres) then
    effnum = 0._8
    DO iso = 1, niso
      id = mapnucl(idiso(iso));   isodata => ldiso(id)
      if (isodata%icat.ne.icat) cycle
      npot = isodata%nsig0; nRiTemp = isodata%nrtemp
      var(1:nRiTemp) = isodata%ri_a(npot,ig, 1:nRiTemp)
      ria = LineIntpol(tempavgsq, nRiTemp, isodata%rtempsq(1:nRiTemp), var(1:nRiTemp))
      effnum = effnum + pnum(iso)*ria
    ENDDO
    effnum = effnum / repria
    if (effnum.gt.0._8) LocalSigtr = LocalSigtr + lv * effnum * myFxr%NDAF(ilv,icat,ig)
    endif

    Siglp(ifxr) = LocalSiglp
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      Sigtr(ifsr) = LocalSigtr
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL
NULLIFY(myFxr,pnum,idiso,isodata,PIN,CellInfo)
END SUBROUTINE
SUBROUTINE SetPlnLsigP_1gMLG(Siglp, Sigtr, lv, Core, Fxr, iz, lCLD, lAIC, PE)
USE TYPEDEF,     ONLY : CoreInfo_Type, FxrInfo_Type, pin_type, cell_type, PE_Type
USE XsLib_Mod,   ONLY : libdata, ldiso, mapnucl
IMPLICIT NONE
TYPE(PE_TYPE) :: PE
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type),POINTER :: Fxr(:, :), myFxr
TYPE(Pin_TYPE), POINTER :: PIN(:)
TYPE(Cell_TYPE), POINTER :: CellInfo(:)
TYPE(libdata), POINTER :: isodata

INTEGER,INTENT(IN) :: iz
REAL,INTENT(IN) :: lv
REAL,POINTER :: Siglp(:), Sigtr(:)
LOGICAL,INTENT(IN) :: lCLD, lAIC

INTEGER :: ipin, ifxr, ifsr, j, i,FxrIdxSt, FsrIdxSt, icel, nlocalFxr, nxy, nFsrInFxr
INTEGER :: niso, id, iso, nThread
INTEGER, POINTER :: idiso(:)
REAL :: Localsiglp, LocalSigtr
REAL, POINTER :: pnum(:)

nxy = Core%nxy
Pin => Core%Pin; CellInfo => Core%CellInfo

nThread = PE%nThread
CALL OMP_SET_NUM_THREADS(nThread)

!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(ipin, FxrIdxSt, FsrIdxSt, icel, nlocalFxr, j, ifxr, myFxr, niso, pnum, idiso, LocalSiglp, LocalSigtr, &
!$OMP         iso, id, isodata, nFsrInFxr, i, ifsr )
!$OMP DO
DO ipin = 1, nxy
  FxrIdxSt = Pin(ipin)%FxrIdxSt; FsrIdxSt = Pin(ipin)%FsrIdxSt
  icel = Pin(ipin)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j -1
    myFxr => Fxr(ifxr, iz)
    niso = myFxr%niso; pnum => myFxr%pnum; idiso => myFxr%idiso
    LocalSiglp = 0._8
    DO iso = 1, niso
      id = mapnucl(idiso(iso));   isodata => ldiso(id)
      !LocalSiglp = LocalSiglp + pnum(iso) * isodata%lamsigp1G
      LocalSiglp = LocalSiglp + pnum(iso) * isodata%sigp
    ENDDO
    LocalSigtr = LocalSiglp
    if (lCLD) then
        if (Fxr(ifxr, iz)%lres.AND.Fxr(ifxr, iz)%lCLD) LocalSigtr = LocalSigtr + lv
    elseif (lAIC) then
        if (Fxr(ifxr, iz)%lres.AND.Fxr(ifxr, iz)%lAIC) LocalSigtr = LocalSigtr + lv
    endif
    Siglp(ifxr) = LocalSiglp
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      Sigtr(ifsr) = LocalSigtr
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL
NULLIFY(myFxr,pnum,idiso,isodata,PIN,CellInfo)
END SUBROUTINE
SUBROUTINE SetPlnLsigP_MLG(Siglp, Sigtr, lv, Core, Fxr, ilv, iz, ig, PE)
USE TYPEDEF,     ONLY : CoreInfo_Type, FxrInfo_Type, pin_type, cell_type , PE_Type
USE XsLib_Mod,   ONLY : libdata, ldiso, mapnucl
IMPLICIT NONE
TYPE(PE_TYPE) :: PE
TYPE(Pin_TYPE), POINTER :: PIN(:)
TYPE(Cell_TYPE), POINTER :: CellInfo(:)
TYPE(libdata), POINTER :: isodata
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type),POINTER :: Fxr(:, :), myFxr

INTEGER,INTENT(IN) :: ilv, iz, ig
REAL,INTENT(IN) :: lv
REAL,POINTER :: Siglp(:), Sigtr(:)

INTEGER :: ipin, nxy, ifxr, ifsr, j, i, FxrIdxSt, FsrIdxSt, icel, nlocalFxr, nFsrInFxr
INTEGER :: id, iso, niso, nThread
INTEGER, POINTER :: idiso(:)
REAL :: Localsiglp, LocalSigtr
REAL,POINTER :: pnum(:)

nxy = Core%nxy
Pin => Core%Pin; CellInfo => Core%CellInfo

nThread = PE%nThread
CALL OMP_SET_NUM_THREADS(nThread)

!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(ipin, FxrIdxSt, FsrIdxSt, icel, nlocalFxr, j, ifxr, myFxr, niso, pnum, idiso, LocalSiglp, LocalSigtr, &
!$OMP         iso, id, isodata, nFsrInFxr, i, ifsr )
!$OMP DO
DO ipin = 1, nxy
  FxrIdxSt = Pin(ipin)%FxrIdxSt; FsrIdxSt = Pin(ipin)%FsrIdxSt
  icel = Pin(ipin)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j -1
    myFxr => Fxr(ifxr, iz)
    niso = myFxr%niso; pnum => myFxr%pnum; idiso => myFxr%idiso
    LocalSiglp = 0._8
    DO iso = 1, niso
      id = mapnucl(idiso(iso));   isodata => ldiso(id)
      !LocalSiglp = LocalSiglp + pnum(iso) * isodata%lamsigp(ig)
      !LocalSiglp = LocalSiglp + pnum(iso) * isodata%sigp
      if (isodata%sigp == 0.d0) then
        LocalSiglp = LocalSiglp + pnum(iso) * isodata%lamsigp(ig)
      else
        LocalSiglp = LocalSiglp + pnum(iso) * isodata%sigp
      end if
    ENDDO
    LocalSigtr = LocalSiglp
    if (Fxr(ifxr, iz)%lres) then
        if((.not.Fxr(ifxr, iz)%lCLD).AND.(.not.Fxr(ifxr, iz)%lAIC)) LocalSigtr = LocalSigtr + lv * myFxr%FnAdj(ig) * myFxr%FtAdj(ilv,ig)
    endif
    Siglp(ifxr) = LocalSiglp
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      Sigtr(ifsr) = LocalSigtr
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL
NULLIFY(myFxr,pnum,idiso,isodata,PIN,CellInfo)
END SUBROUTINE
SUBROUTINE SetPlnLsigP_Dancoff(Siglp, Sigtr, Core, Fxr, iz, PE)
USE TYPEDEF,     ONLY : CoreInfo_Type, FxrInfo_Type, ResVarPin_Type, pin_type, cell_type, PE_Type
USE XsLib_Mod,   ONLY : libdata, ldiso, mapnucl
USE OMP_LIB
IMPLICIT NONE
TYPE(PE_TYPE) :: PE
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type),POINTER :: Fxr(:, :), myFxr
TYPE(Pin_TYPE), POINTER :: PIN(:)
TYPE(ResVarPin_Type), POINTER :: ResVarPin(:,:)
TYPE(Cell_TYPE), POINTER :: CellInfo(:)
TYPE(libdata), POINTER :: isodata

INTEGER,INTENT(IN) :: iz
REAL,POINTER :: Siglp(:), Sigtr(:)

INTEGER :: ipin, ifxr, ifsr, j, i,FxrIdxSt, FsrIdxSt, icel, nlocalFxr, nxy, nFsrInFxr
INTEGER :: niso, id, iso, nThread
INTEGER, POINTER :: idiso(:)
REAL :: Localsiglp, LocalSigtr
REAL, POINTER :: pnum(:)

REAL,PARAMETER :: InfSigt=1E+5

nxy = Core%nxy
Pin => Core%Pin; CellInfo => Core%CellInfo
ResVarPin => Core%ResVarPin

nThread = PE%nThread
CALL OMP_SET_NUM_THREADS(nThread)
!$OMP PARALLEL
!$OMP DO PRIVATE(ipin, FxrIdxSt, FsrIdxSt, icel, nlocalFxr, j, ifxr, myfxr, niso, pnum, idiso, LocalSiglp, iso, id, isodata, LocalSigtr, nFsrInFxr, i, ifsr)
DO ipin = 1, nxy
  FxrIdxSt = Pin(ipin)%FxrIdxSt; FsrIdxSt = Pin(ipin)%FsrIdxSt
  icel = Pin(ipin)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j -1
    myFxr => Fxr(ifxr, iz)
    niso = myFxr%niso;
    pnum => myFxr%pnum;
    idiso => myFxr%idiso
    LocalSiglp = 0._8
    DO iso = 1, niso
      id = mapnucl(idiso(iso));   isodata => ldiso(id)
      LocalSiglp = LocalSiglp + pnum(iso) * isodata%sigp
    ENDDO
    LocalSigtr = LocalSiglp
    Siglp(ifxr) = LocalSiglp
    if (myFxr%lres) then
      if (myFxr%lFUEL) LocalSigtr = InfSigt
    endif
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      Sigtr(ifsr) = LocalSigtr
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL
NULLIFY(myFxr,pnum,idiso,isodata,PIN,CellInfo)
END SUBROUTINE

SUBROUTINE SetPlnLsigP_DancoffAIC(Siglp, Sigtr, Core, Fxr, iz, PE)
USE TYPEDEF,     ONLY : CoreInfo_Type, FxrInfo_Type, ResVarPin_Type, pin_type, cell_type, PE_Type
USE XsLib_Mod,   ONLY : libdata, ldiso, mapnucl
IMPLICIT NONE
TYPE(PE_TYPE) :: PE
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type),POINTER :: Fxr(:, :), myFxr
TYPE(Pin_TYPE), POINTER :: PIN(:)
TYPE(ResVarPin_Type), POINTER :: ResVarPin(:,:)
TYPE(Cell_TYPE), POINTER :: CellInfo(:)
TYPE(libdata), POINTER :: isodata

INTEGER,INTENT(IN) :: iz
REAL,POINTER :: Siglp(:), Sigtr(:)

INTEGER :: ipin, ifxr, ifsr, j, i,FxrIdxSt, FsrIdxSt, icel, nlocalFxr, nxy, nFsrInFxr
INTEGER :: niso, id, iso, nThread
INTEGER, POINTER :: idiso(:)
REAL :: Localsiglp, LocalSigtr
REAL, POINTER :: pnum(:)

REAL,PARAMETER :: InfSigt=1E+5

nxy = Core%nxy
Pin => Core%Pin; CellInfo => Core%CellInfo
ResVarPin => Core%ResVarPin

nThread = PE%nThread
CALL OMP_SET_NUM_THREADS(nThread)

!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(ipin, FxrIdxSt, FsrIdxSt, icel, nlocalFxr, j, ifxr, myFxr, niso, pnum, idiso, LocalSiglp, LocalSigtr, &
!$OMP         iso, id, isodata, nFsrInFxr, i, ifsr )
!$OMP DO
DO ipin = 1, nxy
  FxrIdxSt = Pin(ipin)%FxrIdxSt; FsrIdxSt = Pin(ipin)%FsrIdxSt
  icel = Pin(ipin)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j -1
    myFxr => Fxr(ifxr, iz)
    niso = myFxr%niso; pnum => myFxr%pnum; idiso => myFxr%idiso
    LocalSiglp = 0._8
    DO iso = 1, niso
      id = mapnucl(idiso(iso));   isodata => ldiso(id)
      LocalSiglp = LocalSiglp + pnum(iso) * isodata%sigp
    ENDDO
    LocalSigtr = LocalSiglp
    Siglp(ifxr) = LocalSiglp
    if (myFxr%lres) then
      if (myFxr%lAIC) LocalSigtr = InfSigt
    endif
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      Sigtr(ifsr) = LocalSigtr
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL
NULLIFY(myFxr,pnum,idiso,isodata,PIN,CellInfo)
END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SetSubGrpSrc1g(src, SigLamPot, xstr1g, Core, Fxr, iz, PE)
USE TYPEDEF, ONLY : CoreInfo_Type, FxrInfo_Type, pin_type, cell_type, PE_Type
IMPLICIT NONE
TYPE(PE_TYPE) :: PE
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER  :: Fxr(:, :)
TYPE(Cell_TYPE), POINTER :: CellInfo(:)
TYPE(Pin_TYPE), POINTER :: PIN(:)

INTEGER,INTENT(IN) :: iz
REAL,POINTER,INTENT(IN) :: xstr1g(:), SigLamPot(:)
REAL,POINTER :: src(:)

INTEGER :: FsrIdxSt, FxrIdxSt, icel, nlocalFxr, nFsrInFxr, nxy, ipin, ifxr, ifsr, j, i, nThread
REAL :: xstrinv, localsrc

nxy = Core%nxy
Pin => Core%Pin; CellInfo => Core%CellInfo

nThread = PE%nThread
CALL OMP_SET_NUM_THREADS(nThread)

!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(ipin, FsrIdxSt, FxrIdxSt, icel, nlocalFxr, j, ifxr, nFsrInFxr, ifsr, xstrinv, i, localsrc)
!$OMP DO
DO ipin = 1, nxy
  FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
  icel = Pin(ipin)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j -1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    ifsr = Fxr(ifxr, iz)%FsrIdxSt; xstrinv = 1._8/xstr1g(ifsr)
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      localsrc = SigLamPot(ifxr)
      src(ifsr) = localsrc * xstrinv
    ENDDO   !Fsr in Fxr Sweep
  ENDDO  ! Fxr Sweep
ENDDO !Pin Sweep
!$OMP END DO
!$OMP END PARALLEL
Nullify(Pin,CellInfo)
END SUBROUTINE
SUBROUTINE Set_Dancoff(Siglp, Sigtrfsr, phi, Core, Fxr, iz, nitertot, PE)
USE TYPEDEF,     ONLY : CoreInfo_Type, FxrInfo_Type, pin_type, cell_type, ResVarPin_Type, PE_Type
USE SUBGRP_MOD,  ONLY : Reset1Dpingeom
IMPLICIT NONE
TYPE(PE_TYPE) :: PE
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type),POINTER :: Fxr(:, :), myFxr
TYPE(Pin_TYPE), POINTER :: PIN(:)
TYPE(Cell_TYPE), POINTER :: CellInfo(:)
TYPE(ResVarPin_Type), POINTER :: ResVarPin(:,:),RP

INTEGER,INTENT(IN) :: iz
REAL,POINTER :: Siglp(:), Sigtrfsr(:), phi(:)

INTEGER :: ipin,ifxr,ifsr,j,i,k,l,FxrIdxSt,FsrIdxSt,icel,nlocalFxr,nxy,nFsrInFxr,niter,nitertot
INTEGER :: nThread
REAL :: phisumf,phifxr,siglpsum,rrsum,vol,volsum,phivol,rrsumlocal
REAL :: Dancoff,eR,Sigtrfxr(100),siglpfxr(100)

nxy = Core%nxy
Pin => Core%Pin; CellInfo => Core%CellInfo
ResVarPin => Core%ResVarPin


nThread = PE%nThread
CALL OMP_SET_NUM_THREADS(nThread)

nitertot=0
!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(ipin,FxrIdxSt,FsrIdxSt,RP,icel,nlocalFxr,Sigtrfxr,siglpfxr,rrsum,volsum,siglpsum,phisumf, &
!$OMP         j,k,ifxr,nFsrInFxr,myFxr,phifxr,rrsumlocal,i,ifsr,l,vol,phivol,Dancoff,eR,niter)
!$OMP DO REDUCTION(+:nitertot)
DO ipin = 1, nxy
  FxrIdxSt = Pin(ipin)%FxrIdxSt; FsrIdxSt = Pin(ipin)%FsrIdxSt
  if (.not.ResVarPin(ipin,iz)%lres) cycle
  RP => ResVarPin(ipin,iz)
  icel = Pin(ipin)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
  if (CellInfo(icel)%lAIC) cycle
  Sigtrfxr(1:nlocalFxr)=0._8; siglpfxr(1:nlocalFxr)=0._8
  rrsum=0._8; volsum = 0._8; siglpsum = 0._8; phisumf=0._8
  DO j = 1, nLocalFxr
    k = nLocalFxr - j + 1
    ifxr = FxrIdxSt + j -1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    myFxr => Fxr(ifxr, iz)
    if (myFxr%lfuel) then
      phifxr = 0._8; rrsumlocal=0._8
      DO i = 1, nFsrInFxr
        l = Cellinfo(icel)%MapFxr2FsrIdx(i, j); ifsr = FsrIdxSt + l - 1
        vol = CellInfo(icel)%vol(l)
        phivol = phi(ifsr) * vol
        rrsumlocal = rrsumlocal + Sigtrfsr(ifsr) * phivol
        phifxr = phifxr + phivol
      ENDDO
      rrsum = rrsum + rrsumlocal
      Sigtrfxr(k) = rrsumlocal/phifxr
      phisumf = phisumf + phifxr
      siglpsum = siglpsum + Siglp(ifxr) * phifxr
      volsum = volsum + myFxr%area
    else
      l = Cellinfo(icel)%MapFxr2FsrIdx(1, j); ifsr = FsrIdxSt + l - 1
      Sigtrfxr(k)=Sigtrfsr(ifsr)
    endif
    siglpfxr(k)=Siglp(ifxr)
  ENDDO
  if (volsum.gt.0._8) then
    siglpsum = siglpsum / phisumf
    rrsum = rrsum / volsum
    Dancoff = RP%lbar * (rrsum - siglpsum)
    RP%Dancoff = Dancoff
  endif
  call EquivRsearch(Core,icel,eR,niter,Dancoff,RP%rad_cp(1:nlocalFxr),Sigtrfxr(1:nlocalFxr),siglpfxr(1:nlocalFxr),nlocalFxr,RP%lbar)
  nitertot = nitertot + niter
  call Reset1Dpingeom(RP,eR,Core,icel,nlocalfxr)
ENDDO
!$OMP END DO
!$OMP END PARALLEL

NULLIFY(myFxr,ResVarPin,PIN,CellInfo)
END SUBROUTINE
SUBROUTINE Set_DancoffAIC(Siglp, Sigtrfsr, phi, Core, Fxr, iz, nitertot, PE)
USE TYPEDEF,     ONLY : CoreInfo_Type, FxrInfo_Type, pin_type, cell_type, ResVarPin_Type, PE_Type
USE SUBGRP_MOD,  ONLY : Reset1Dpingeom
IMPLICIT NONE
TYPE(PE_TYPE) :: PE
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type),POINTER :: Fxr(:, :), myFxr
TYPE(Pin_TYPE), POINTER :: PIN(:)
TYPE(Cell_TYPE), POINTER :: CellInfo(:)
TYPE(ResVarPin_Type), POINTER :: ResVarPin(:,:),RP

INTEGER,INTENT(IN) :: iz
REAL,POINTER :: Siglp(:), Sigtrfsr(:), phi(:)

INTEGER :: ipin,ifxr,ifsr,j,i,k,l,FxrIdxSt,FsrIdxSt,icel,nlocalFxr,nxy,nFsrInFxr,niter,nitertot
INTEGER :: nThread
REAL :: phisumf,phifxr,siglpsum,rrsum,vol,volsum,phivol,rrsumlocal
REAL :: Dancoff,eR,Sigtrfxr(100),siglpfxr(100)

nxy = Core%nxy
Pin => Core%Pin; CellInfo => Core%CellInfo
ResVarPin => Core%ResVarPin


nThread = PE%nThread
CALL OMP_SET_NUM_THREADS(nThread)

nitertot=0
!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(ipin,FxrIdxSt,FsrIdxSt,RP,icel,nlocalFxr,Sigtrfxr,siglpfxr,rrsum,volsum,siglpsum,phisumf, &
!$OMP         j,k,ifxr,nFsrInFxr,myFxr,phifxr,rrsumlocal,i,ifsr,l,vol,phivol,Dancoff,eR,niter)
!$OMP DO REDUCTION(+:nitertot)
DO ipin = 1, nxy
  FxrIdxSt = Pin(ipin)%FxrIdxSt; FsrIdxSt = Pin(ipin)%FsrIdxSt
  if (.not.ResVarPin(ipin,iz)%lres) cycle
  RP => ResVarPin(ipin,iz)
  icel = Pin(ipin)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
  if (.not.CellInfo(icel)%lAIC) cycle
  Sigtrfxr(1:nlocalFxr)=0._8; siglpfxr(1:nlocalFxr)=0._8
  rrsum=0._8; volsum = 0._8; siglpsum = 0._8; phisumf=0._8
  DO j = 1, nLocalFxr
    k = nLocalFxr - j + 1
    ifxr = FxrIdxSt + j -1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    myFxr => Fxr(ifxr, iz)
    if (myFxr%lAIC) then
      phifxr = 0._8; rrsumlocal=0._8
      DO i = 1, nFsrInFxr
        l = Cellinfo(icel)%MapFxr2FsrIdx(i, j); ifsr = FsrIdxSt + l - 1
        vol = CellInfo(icel)%vol(l)
        phivol = phi(ifsr) * vol
        rrsumlocal = rrsumlocal + Sigtrfsr(ifsr) * phivol
        phifxr = phifxr + phivol
      ENDDO
      rrsum = rrsum + rrsumlocal
      Sigtrfxr(k) = rrsumlocal/phifxr
      phisumf = phisumf + phifxr
      siglpsum = siglpsum + Siglp(ifxr) * phifxr
      volsum = volsum + myFxr%area
    else
      l = Cellinfo(icel)%MapFxr2FsrIdx(1, j); ifsr = FsrIdxSt + l - 1
      Sigtrfxr(k)=Sigtrfsr(ifsr)
    endif
    siglpfxr(k)=Siglp(ifxr)
  ENDDO
  if (volsum.gt.0._8) then
    siglpsum = siglpsum / phisumf
    rrsum = rrsum / volsum
    Dancoff = RP%lbar * (rrsum - siglpsum)
    RP%Dancoff = Dancoff
  endif
  call EquivRsearch(Core,icel,eR,niter,Dancoff,RP%rad_cp(1:nlocalFxr),Sigtrfxr(1:nlocalFxr),siglpfxr(1:nlocalFxr),nlocalFxr,RP%lbar)
  nitertot = nitertot + niter
  call Reset1Dpingeom(RP,eR,Core,icel,nlocalfxr)
ENDDO
!$OMP END DO
!$OMP END PARALLEL

NULLIFY(myFxr,ResVarPin,PIN,CellInfo)
END SUBROUTINE
SUBROUTINE Reset1Dpingeom(RP,eR,Core,icel,nlocalfxr)
  USE PARAM,    ONLY : INVPI
  USE TYPEDEF,  ONLY : CoreInfo_Type,ResVarPin_Type, Cell_TYPE
  USE CP_mod
  IMPLICIT NONE
  TYPE(CoreInfo_Type) :: Core
  integer :: nlocalfxr,icel
  real :: eR
  TYPE(ResVarPin_Type), POINTER :: RP
  TYPE(Cell_TYPE), POINTER :: CellInfo(:)

  integer :: i, j, idx, m, k, n
  real :: SUMX, Y, Y22, Qsurf, vol, subvol, delr(nlocalfxr), r2(nlocalfxr), X(nlocalfxr,nlocalfxr,ngauss)

  CellInfo => Core%CellInfo
  nlocalfxr = CellInfo(icel)%nFxr
  RP%EquiRad = eR
  Qsurf = hpi * eR
  RP%invsurf4 = 1._8/Qsurf
  vol = CellInfo(icel)%fuelgapcldvol
  subvol = (PI*eR*eR-vol)*CellInfo(icel)%invnmodfxr
  idx = CellInfo(icel)%cldfxridx
  RP%vol(1:idx)=CellInfo(icel)%fxrvol(1:idx)
  do i = 1, CellInfo(icel)%nmodfxr
    vol = vol + subvol
    idx = idx + 1
    RP%rad_cp(idx) = dsqrt(vol*INVPI)
    RP%vol(idx) = subvol
  enddo
  RP%Qsurfvol(1) = RP%vol(1) * Qsurf
  delr(1) = RP%rad_cp(1)
  r2(1) = RP%rad_cp(1) * RP%rad_cp(1)
  do i = 2, nlocalFxr
    RP%Qsurfvol(i) = RP%vol(i) * Qsurf
    delr(i) = RP%rad_cp(i) - RP%rad_cp(i-1)
    r2(i) = RP%rad_cp(i) * RP%rad_cp(i)
  enddo

  DO m=1,ngauss
    DO k=1,nlocalFxr
      X(:,k,m) = 0._8
      Y=RP%rad_cp(k)-delr(k)*Quad(m); Y22=Y*Y
      X(k,k,m)=dsqrt(r2(k)-Y22)
      SUMX=X(k,k,m)
      DO n=k+1,nlocalFxr
        X(n,k,m)=dsqrt(r2(n)-Y22)-SUMX
        SUMX=SUMX+X(n,k,m)
      ENDDO
    ENDDO
  ENDDO

  do i = 1, nlocalFxr
    RP%delr(i) = delr(i) * 4._8
  enddo
  do m=1,ngauss
    do j = 1, nlocalFxr
      do i = 1, nlocalFxr
        RP%X(i,j,m)=X(i,j,m)
      enddo
    enddo
  enddo
END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION SubGrpFspErr(phis1g, phis1gd, nfsr, PE)
USE TYPEDEF, only : PE_Type
USE PARAM
IMPLICIT NONE
TYPE(PE_TYPE) :: PE
REAL :: SubGrpFspErr
REAL, POINTER :: phis1g(:), phis1gd(:)
INTEGER :: nfsr
INTEGER :: i, nThread
REAL :: errmax, err

nThread = PE%nThread
CALL OMP_SET_NUM_THREADS(nThread)

errmax = 0._8
!$OMP PARALLEL DO REDUCTION(max:errmax)      &
!$OMP PRIVATE(i, err)
DO i = 1, nfsr
  IF(phis1g(i) .lt. 0._8) CYCLE
  err = abs((phis1g(i)-phis1gd(i))/phis1gd(i))
  errmax = max(err,errmax)
ENDDO
!$OMP END PARALLEL DO
SubGrpFspErr = errmax
END FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EquipXSGen(phis1g, SigLamPot, xstr1g, ilv, irc, Core, Fxr, iz, ig, PE)
USE Typedef,     ONLY : CoreInfo_Type, FxrInfo_Type, pin_type, ResVarPin_Type, cell_type, PE_Type
IMPLICIT NONE
TYPE(PE_TYPE) :: PE
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :), myFxr
TYPE(Cell_TYpe), POINTER :: CellInfo(:)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(ResVarPin_Type), POINTER :: ResVarPin(:,:)

INTEGER,INTENT(IN) :: iz, ilv, irc, ig
REAL, POINTER,INTENT(IN) :: phis1g(:), SigLamPot(:), xstr1g(:)
INTEGER :: FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr, nxy, ipin, ifxr, ifsr, icel, i, j, l, nThread
REAL :: vol, volsum, phisum, xseq, maclp, maclv, avgphisum, avgvolsum, maclpavg, maclvavg

nxy = Core%nxy
Pin => Core%Pin; CellInfo => Core%CellInfo
ResVarPin => Core%ResVarPin

nThread = PE%nThread
CALL OMP_SET_NUM_THREADS(nThread)

!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(ipin, FsrIdxSt, FxrIdxSt, icel, nlocalFxr, avgphisum, avgvolsum, maclvavg, maclpavg, j, ifxr, nFsrInFxr, &
!$OMP         myFxr, ifsr, maclp, maclv, phisum, volsum, i, l, vol, xseq )
!$OMP DO
DO ipin = 1, nxy
  FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
  icel = Pin(ipin)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
  if (.not.CellInfo(icel)%lres) CYCLE
  avgphisum = 0._8; avgvolsum = 0._8
  maclvavg = 0._8; maclpavg = 0._8
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j -1; nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    myFxr => Fxr(ifxr, iz)
    if (.NOT. myFxr%lres) CYCLE
    ifsr = myFxr%FsrIdxSt
    maclp = SigLamPot(ifxr)
    maclv = xstr1g(ifsr) - maclp
    if (maclv.eq.0._8) cycle
    phisum = 0; volsum = 0
    DO i = 1, nFsrInFxr
      l = Cellinfo(icel)%MapFxr2FsrIdx(i, j); vol = CellInfo(icel)%vol(l)
      ifsr = FsrIdxSt + l - 1
      phisum = phisum + phis1g(ifsr) * vol; volsum = volsum + vol
    ENDDO
    maclvavg = maclvavg + maclv * phisum
    maclpavg = maclpavg + maclp * phisum
    avgphisum = avgphisum + phisum
    avgvolsum = avgvolsum + volsum
    phisum = phisum / volsum
    if (abs(phisum-1._8).lt.1E-10) then
      xseq = 1E+10
    else
      xseq = - maclp + maclv * phisum/(1._8 - phisum)
    endif
    myFxr%xseq(ilv, irc, ig) = xseq
  ENDDO
  if (maclpavg.eq.0._8) cycle
  if (.not.ResVarPin(ipin,iz)%lres) CYCLE
  maclvavg = maclvavg / avgphisum
  maclpavg = maclpavg / avgphisum
  avgphisum = avgphisum / avgvolsum
  if (abs(avgphisum-1._8).lt.1E-10) then
    xseq = 1E+10
  else
    xseq = - maclpavg + maclvavg * avgphisum/(1._8 - avgphisum)
  endif
  ResVarPin(ipin,iz)%avgxseq(ilv, irc, ig) = xseq
ENDDO
!$OMP END DO
!$OMP END PARALLEL
NULLIFY(PIN,ResVarPin,CellInfo,myFxr)
END SUBROUTINE
SUBROUTINE EquipXSGen_1gMLG(phis1g, SigLamPot, xstr1g, ilv, Core, Fxr, iz, lCLD, lAIC, PE)
USE Typedef,     ONLY : CoreInfo_Type, FxrInfo_Type, pin_type, ResVarPin_Type, cell_type, PE_Type
IMPLICIT NONE
TYPE(PE_TYPE) :: PE
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :), myFxr
TYPE(Cell_TYpe), POINTER :: CellInfo(:)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(ResVarPin_Type), POINTER :: ResVarPin(:,:)

INTEGER,INTENT(IN) :: iz, ilv
REAL,POINTER,INTENT(IN) :: phis1g(:), SigLamPot(:), xstr1g(:)
LOGICAL,INTENT(IN) :: lCLD, lAIC
INTEGER :: ipin, ifxr, ifsr, icel, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr, nxy, i, j, l, nThread
REAL :: vol, volsum, phisum, xseq, maclp, maclv, avgphisum, avgvolsum, maclpavg, maclvavg

nxy = Core%nxy
Pin => Core%Pin; CellInfo => Core%CellInfo
ResVarPin => Core%ResVarPin

nThread = PE%nThread
CALL OMP_SET_NUM_THREADS(nThread)

!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(ipin, FsrIdxSt, FxrIdxSt, icel, nlocalFxr, j, ifxr, nFsrInFxr, myFxr, ifsr, maclp, maclv, &
!$OMP         phisum, volsum, i, l, vol, xseq, avgphisum, avgvolsum, maclpavg, maclvavg )
!$OMP DO
DO ipin = 1, nxy
  FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
  icel = Pin(ipin)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
  if (.not.CellInfo(icel)%lres) CYCLE
  IF (lCLD) THEN
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j -1; nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
      myFxr => Fxr(ifxr, iz); ifsr = myFxr%FsrIdxSt
      IF(.NOT.myFxr%lres) CYCLE
      ifsr = myFxr%FsrIdxSt
      maclp = SigLamPot(ifxr)
      maclv = xstr1g(ifsr) - maclp
      if (maclv.eq.0._8) cycle
      phisum = 0; volsum = 0
      DO i = 1, nFsrInFxr
        l = Cellinfo(icel)%MapFxr2FsrIdx(i, j); vol = CellInfo(icel)%vol(l)
        ifsr = FsrIdxSt + l - 1
        phisum = phisum + phis1g(ifsr) * vol; volsum = volsum + vol
      ENDDO
      phisum = phisum / volsum
      if (abs(phisum-1._8).lt.1E-10) then
        xseq = 1E+10
      else
        xseq = - maclp + maclv * phisum/(1._8 - phisum)
      endif
      myFxr%xseq_c_1g(ilv) = xseq
    ENDDO
  ELSEIF (lAIC) THEN
    if (.not.CellInfo(icel)%lAIC) CYCLE
    maclpavg = 0._8; maclvavg = 0._8; avgphisum = 0._8; avgvolsum = 0._8
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j -1; nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
      myFxr => Fxr(ifxr, iz)
      IF(.NOT.myFxr%lres) CYCLE
      ifsr = myFxr%FsrIdxSt
      maclp = SigLamPot(ifxr)
      maclv = xstr1g(ifsr) - maclp
      if (maclv.eq.0._8) cycle
      phisum = 0; volsum = 0
      DO i = 1, nFsrInFxr
        l = Cellinfo(icel)%MapFxr2FsrIdx(i, j); vol = CellInfo(icel)%vol(l)
        ifsr = FsrIdxSt + l - 1
        phisum = phisum + phis1g(ifsr) * vol; volsum = volsum + vol
      ENDDO
      maclpavg = maclpavg + maclp * phisum
      maclvavg = maclvavg + maclv * phisum
      avgphisum = avgphisum + phisum
      avgvolsum = avgvolsum + volsum
      phisum = phisum / volsum
      if (abs(phisum-1._8).lt.1E-10) then
        xseq = 1E+10
      else
        xseq = - maclp + maclv * phisum/(1._8 - phisum)
      endif
      myFxr%xseq_f_1g(ilv) = xseq
    ENDDO
    if (maclpavg.eq.0._8) cycle
    maclpavg = maclpavg / avgphisum
    maclvavg = maclvavg / avgphisum
    avgphisum = avgphisum / avgvolsum
    if (abs(avgphisum-1._8).lt.1E-10) then
      xseq = 1E+10
    else
      xseq = - maclpavg + maclvavg * avgphisum/(1._8 - avgphisum)
    endif
    ResVarPin(ipin,iz)%avgxseq_1g(ilv) = xseq
  ENDIF
ENDDO
!$OMP END DO
!$OMP END PARALLEL
NULLIFY(PIN,ResVarPin,CellInfo,myFxr)
END SUBROUTINE
SUBROUTINE EquipXSGen_MLG(phis1g, SigLamPot, xstr1g, ilv, Core, Fxr, iz, ig, PE)
USE Typedef,     ONLY : CoreInfo_Type, FxrInfo_Type, pin_type, ResVarPin_Type, cell_type, PE_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :), myFxr
TYPE(PE_TYPE) :: PE
TYPE(Cell_TYpe), POINTER :: CellInfo(:)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(ResVarPin_Type), POINTER :: ResVarPin(:,:)

INTEGER,INTENT(IN) :: iz, ilv, ig
REAL, POINTER,INTENT(IN) :: phis1g(:), SigLamPot(:), xstr1g(:)
INTEGER :: ipin, ifxr, ifsr, icel, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr, nxy, i, j, l, nThread
REAL :: vol, volsum, phisum, xseq, maclp, maclv, avgphisum, avgvolsum, maclpavg, maclvavg

nxy = Core%nxy
Pin => Core%Pin; CellInfo => Core%CellInfo
ResVarPin => Core%ResVarPin

nThread = PE%nThread
CALL OMP_SET_NUM_THREADS(nThread)

!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(ipin, FsrIdxSt, FxrIdxSt, icel, nlocalFxr, avgphisum, avgvolsum, maclvavg, maclpavg, j, ifxr, nFsrInFxr, &
!$OMP         myFxr, ifsr, maclp, maclv, phisum, volsum, i, l, vol, xseq )
!$OMP DO
DO ipin = 1, nxy
  FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
  icel = Pin(ipin)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
  if (.not.CellInfo(icel)%lres) CYCLE
  avgphisum = 0._8; avgvolsum = 0._8
  maclvavg = 0._8; maclpavg = 0._8
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j -1; nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    myFxr => Fxr(ifxr, iz); ifsr = myFxr%FsrIdxSt
    IF(.NOT.myFxr%lres) CYCLE
    ifsr = myFxr%FsrIdxSt
    maclp = SigLamPot(ifxr)
    maclv = xstr1g(ifsr) - maclp
    if (maclv.eq.0._8) cycle
    phisum = 0; volsum = 0
    DO i = 1, nFsrInFxr
      l = Cellinfo(icel)%MapFxr2FsrIdx(i, j); vol = CellInfo(icel)%vol(l)
      ifsr = FsrIdxSt + l - 1
      phisum = phisum + phis1g(ifsr) * vol; volsum = volsum + vol
    ENDDO
    maclpavg = maclpavg + maclp * phisum
    maclvavg = maclvavg + maclv * phisum
    avgphisum = avgphisum + phisum
    avgvolsum = avgvolsum + volsum
    phisum = phisum / volsum
    if (abs(phisum-1._8).lt.1E-10) then
      xseq = 1E+10
    else
      xseq = - maclp + maclv * phisum/(1._8 - phisum)
    endif
    !xseq = - maclp + maclv * phisum/(1._8 - phisum)
    myFxr%xseq_f_mg(ilv,ig) = xseq
  ENDDO
  if (maclpavg.eq.0._8) cycle
  maclpavg = maclpavg / avgphisum
  maclvavg = maclvavg / avgphisum
  avgphisum = avgphisum / avgvolsum
  if (abs(avgphisum-1._8).lt.1E-10) then
    xseq = 1E+10
  else
    xseq = - maclpavg + maclvavg * avgphisum/(1._8 - avgphisum)
  endif
!  xseq = - maclpavg + maclvavg * avgphisum/(1._8 - avgphisum)
  ResVarPin(ipin,iz)%avgxseq_mg(ilv,ig) = xseq
ENDDO
!$OMP END DO
!$OMP END PARALLEL
NULLIFY(PIN,ResVarPin,CellInfo,myFxr)
END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE UpdtNDAF(id, ilv, irc, Core, Fxr, iz, ig, PE)
USE Typedef,      ONLY : CoreInfo_Type, FxrInfo_Type, Cell_TYpe, Pin_Type, PE_Type
USE XsLib_Mod,    ONLY : libdata, ldiso, mapnucl, mapnuclres
USE XsUtil_mod,   ONLY : LineIntPol, LineIntPol2
USE TH_Mod,       ONLY : GetPinFuelTemp
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(PE_TYPE) :: PE
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :), myFxr
TYPE(Cell_TYpe), POINTER :: CellInfo(:)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(libdata), POINTER :: isodata,jsodata

INTEGER,INTENT(IN) :: id, ilv, irc, iz, ig
INTEGER :: ipin, nxy, ifxr, FxrIdxSt, nlocalFxr, icel, j, niso, jd, iso, jdres, i, npot, nRiTemp, nThread
INTEGER, POINTER :: idiso(:)
INTEGER, PARAMETER :: nmaxRITEMP = 20, nMaxVar = 200
REAL :: TempAvgsq, Tempsq, siglp, ind, micsigb, micsig0, sigbsq, xdat(nmaxRITEMP), ydat(nmaxRITEMP), I_reg, I_avg, var(nMaxVar)
REAL, POINTER :: pnum(:)

nxy = Core%nxy
Pin => Core%Pin; CellInfo => Core%CellInfo
isodata => ldiso(id)

nThread = PE%nThread
CALL OMP_SET_NUM_THREADS(nThread)

!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(ipin, FxrIdxSt, icel, nlocalFxr, TempAvgsq, j, ifxr, myFxr, Tempsq, niso, idiso, pnum, siglp, ind, iso, jd, &
!$OMP         jsodata, jdres, micsigb, micsig0, sigbsq, npot, nRiTemp, i, xdat, ydat, I_reg, I_avg, var )
!$OMP DO
DO ipin = 1, nxy
  FxrIdxSt = Pin(ipin)%FxrIdxSt
  icel = Pin(ipin)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
  if (.not.CellInfo(icel)%lres) CYCLE
  TempAvgsq = dsqrt(GetPinFuelTemp(Core, Fxr, iz, ipin))
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j - 1
    myFxr => Fxr(ifxr, iz)
    if (.not.myFxr%lres) cycle
    if (myFxr%lCLD) cycle
    if (myFxr%xseq(ilv,irc,ig).eq.0._8) cycle

    Tempsq = dsqrt(myFxr%temp)
    niso = myFxr%niso; idiso => myFxr%idiso; pnum => myFxr%pnum

    siglp = 0._8; ind = 0._8
    DO iso = 1, niso
      jd = mapnucl(idiso(iso));   jsodata => ldiso(jd)
      jdres = mapnuclres(idiso(iso),iz)
      siglp = siglp + pnum(iso) * jsodata%lamsigp(ig)
      if (jdres.eq.irc) ind = ind + pnum(iso)
    ENDDO
    micsigb = (myFxr%xseq(ilv,irc,ig) + siglp) / ind / myFxr%NDAF(ilv,irc,ig)
    micsig0 = micsigb - isodata%lamsigp(ig)
    if (micsig0.lt.1E-6) cycle
    sigbsq = dsqrt(micsig0)
    npot = isodata%nsig0; nRiTemp = isodata%nrtemp
    DO i = 1, nRiTemp
      xdat(i) = isodata%rtempsq(i)
      var(1:npot) = isodata%ri_a(1:npot, ig, i)
      ydat(i) = LineIntPol2(sigbsq ,npot, isodata%sig0sq(1:npot), var(1:npot))
    ENDDO
    I_reg = LineIntpol(tempsq, nRiTemp, xdat(1:nRiTemp), ydat(1:nRiTemp))
    I_avg = LineIntpol(TempAvgsq, nRiTemp, xdat(1:nRiTemp), ydat(1:nRiTemp))

    myFxr%NDAF(ilv,irc,ig) = I_reg/I_avg

  ENDDO ! DO j = 1, nLocalFxr
ENDDO ! DO ipin = 1, nxy
!$OMP END DO
!$OMP END PARALLEL
NULLIFY(PIN,CellInfo,myFxr,isodata,jsodata,pnum,idiso)
END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE UpdtNDAF_CAT(ilv, icat, Core, Fxr, iz, ig, PE)
USE Typedef,      ONLY : CoreInfo_Type, FxrInfo_Type, Cell_TYpe, Pin_Type, PE_Type
USE XsLib_Mod,    ONLY : libdata, ldiso, mapnucl, ResoCat
USE XsUtil_mod,   ONLY : LineIntPol, LineIntPol2
USE TH_Mod,       ONLY : GetPinFuelTemp
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(PE_TYPE) :: PE
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :), myFxr
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(libdata), POINTER :: isodata, repisodata

INTEGER,INTENT(IN) :: ilv, icat, iz, ig
INTEGER :: ipin, nxy, ifxr, FxrIdxSt, nlocalFxr, icel, j
INTEGER :: repnid, repid, niso, iso, i, id, npot, nRiTemp, nThread
INTEGER, POINTER :: idiso(:)
REAL :: TempAvgsq, Tempsq, siglp, micsigb, micsig0, sigbsq, repria, ria, effnum, xdat(20), ydat(20), I_reg, I_avg, var(200)
REAL, POINTER :: pnum(:)

nxy = Core%nxy
Pin => Core%Pin; CellInfo => Core%CellInfo

nThread = PE%nThread
CALL OMP_SET_NUM_THREADS(nThread)

repnid = ResoCat(icat)%repid; repid = mapnucl(repnid)
repisodata => ldiso(repid)

!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(ipin, FxrIdxSt, icel, nlocalFxr, TempAvgsq, npot, nRiTemp, repria, j, ifxr, myFxr, Tempsq, niso, idiso, pnum, &
!$OMP         effnum, siglp, iso, id, isodata, ria, micsigb, micsig0, sigbsq, i, xdat, ydat, I_reg, I_avg, var )
!$OMP DO
DO ipin = 1, nxy
  FxrIdxSt = Pin(ipin)%FxrIdxSt
  icel = Pin(ipin)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
  if (.not.CellInfo(icel)%lres) CYCLE
  TempAvgsq = dsqrt(GetPinFuelTemp(Core, Fxr, iz, ipin))
  npot = repisodata%nsig0; nRiTemp = repisodata%nrtemp
  var(1:nRiTemp) = repisodata%ri_a(npot,ig, 1:nRiTemp)
  repria = LineIntpol(TempAvgsq, nRiTemp, repisodata%rtempsq(1:nRiTemp), var(1:nRiTemp))

  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j - 1
    myFxr => Fxr(ifxr, iz)
    if (.not.myFxr%lres) cycle
    if (myFxr%lCLD) cycle
    if (myFxr%xseq(ilv,icat,ig).eq.0._8) cycle

    Tempsq = dsqrt(myFxr%temp)
    niso = myFxr%niso; idiso => myFxr%idiso; pnum => myFxr%pnum

    effnum = 0._8; siglp = 0._8
    DO iso = 1, niso
      id = mapnucl(idiso(iso));   isodata => ldiso(id)
      siglp = siglp + pnum(iso) * isodata%lamsigp(ig)
      if (icat.ne.isodata%icat) cycle
      npot = isodata%nsig0; nRiTemp = isodata%nrtemp
      var(1:nRiTemp) = isodata%ri_a(npot,ig, 1:nRiTemp)
      ria = LineIntpol(TempAvgsq, nRiTemp, isodata%rtempsq(1:nRiTemp), var(1:nRiTemp) )
      effnum = effnum + pnum(iso) * ria
    ENDDO
    effnum = effnum / repria

    micsigb = (myFxr%xseq(ilv,icat,ig) + siglp) / effnum / myFxr%NDAF(ilv,icat,ig)
    micsig0 = micsigb - repisodata%lamsigp(ig)
    if (micsig0.lt.1E-6) cycle
    sigbsq = dsqrt(micsig0)
    npot = repisodata%nsig0; nRiTemp = repisodata%nrtemp
    DO i = 1, nRiTemp
      xdat(i) = repisodata%rtempsq(i)
      var(1:npot) = repisodata%ri_a(1:npot, ig, i)
      ydat(i) = LineIntPol2(sigbsq ,npot, repisodata%sig0sq(1:npot), var(1:npot))
    ENDDO
    I_reg = LineIntpol(tempsq, nRiTemp, xdat(1:nRiTemp), ydat(1:nRiTemp))
    I_avg = LineIntpol(TempAvgsq, nRiTemp, xdat(1:nRiTemp), ydat(1:nRiTemp))

    myFxr%NDAF(ilv,icat,ig) = I_reg / I_avg

  ENDDO ! DO j = 1, nLocalFxr
ENDDO ! DO ipin = 1, nxy
!$OMP END DO
!$OMP END PARALLEL
NULLIFY(PIN,CellInfo,myFxr,isodata,repisodata,pnum,idiso)
END SUBROUTINE
SUBROUTINE UpdtFnAdj(Core, Fxr, ig, iz, PE)
! NUMBER DENSITY Correction factor for 2D MOC MLG (NDCF compared to core RI avergae)
USE Typedef,      ONLY : CoreInfo_Type, FxrInfo_Type, Cell_TYpe, Pin_Type, ResVarPin_Type, PE_Type
USE XsLib_Mod,    ONLY : libdata, ldiso, mapnucl
USE XsUtil_mod,   ONLY : LineIntpol
USE TH_Mod,       ONLY : GetPinFuelTemp
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(PE_TYPE) :: PE
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :), myFxr
TYPE(Cell_TYpe), POINTER :: CellInfo(:)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(ResVarPin_Type), POINTER :: ResVarPin(:,:)
TYPE(libdata), POINTER :: isodata

INTEGER,INTENT(IN) :: iz, ig
INTEGER :: nxy, ipin, FxrIdxSt, icel, nlocalFxr, j, ifxr, niso, iso, jd, nRiTemp, npot, nThread
INTEGER, POINTER :: idiso(:)
REAL :: TempAvgsq, N, Nsum, Npin, areasum, pinareasum, ria, var(200)
REAL, POINTER :: pnum(:), NresPinAvg(:)

nxy = Core%nxy
Pin => Core%Pin; CellInfo => Core%CellInfo
ResVarPin => Core%ResVarPin
allocate(NresPinAvg(nxy))

nThread = PE%nThread
CALL OMP_SET_NUM_THREADS(nThread)

Nsum=0._8 ! core average macro RI
areasum=0._8
!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(ipin,  FxrIdxSt, icel, nlocalFxr, TempAvgsq, Npin, pinareasum, j, ifxr, myFxr, &
!$OMP         niso, idiso, pnum, N, iso, jd, isodata, npot, nRiTemp, ria, var )
!$OMP DO REDUCTION(+:Nsum, areasum)
DO ipin = 1, nxy
  FxrIdxSt = Pin(ipin)%FxrIdxSt
  icel = Pin(ipin)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
  if (.not.CellInfo(icel)%lres) CYCLE
  if (CellInfo(icel)%lAIC) CYCLE
  TempAvgsq = dsqrt(GetPinFuelTemp(Core, Fxr, iz, ipin))
  Npin = 0._8 ! RI*pnum*area summation in Pin
  pinareasum = 0._8
  NresPinAvg(ipin) = 0._8 ! Pin average macro RI
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j - 1
    myFxr => Fxr(ifxr, iz)
    if (.not.myFxr%lres) cycle
    if (myFxr%lCLD) cycle
    niso = myFxr%niso; idiso => myFxr%idiso; pnum => myFxr%pnum
    N = 0._8! macro RI(RI*pnum) summation in FXR
    DO iso = 1, niso
      jd = mapnucl(idiso(iso)); isodata => ldiso(jd)
      if (.not.isodata%lreso) cycle
      npot = isodata%nsig0; nRiTemp = isodata%nrtemp
      var(1:nRiTemp) = isodata%ri_a(npot,ig, 1:nRiTemp) ! infinite RI with Temperatures
      ria = LineIntpol(TempAvgsq, nRiTemp, isodata%rtempsq(1:nRiTemp), var(1:nRiTemp))
      N = N + pnum(iso) * ria
    ENDDO
    Npin = Npin + N * myFxr%area
    pinareasum = pinareasum + myFxr%area
    myFxr%FnAdj(ig) = N
  ENDDO ! DO j = 1, nLocalFxr
  if (pinareasum.gt.0._8) NresPinAvg(ipin) = Npin / pinareasum
  areasum = areasum + pinareasum
  Nsum = Nsum + Npin
ENDDO ! DO ipin = 1, nxy
!$OMP END DO
!$OMP END PARALLEL
Nsum = Nsum / areasum

!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(ipin, FxrIdxSt, icel, nlocalFxr, j, ifxr, myFxr)
!$OMP DO
DO ipin = 1, nxy
  FxrIdxSt = Pin(ipin)%FxrIdxSt
  icel = Pin(ipin)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
  if (.not.ResVarPin(ipin,iz)%lres) CYCLE
  if (.not.CellInfo(icel)%lres) CYCLE
  if (CellInfo(icel)%lAIC) CYCLE
  ResVarPin(ipin,iz)%FnAdj(ig) = NresPinAvg(ipin) / Nsum
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j - 1
    myFxr => Fxr(ifxr, iz)
    if (.not.myFxr%lres) cycle
    if (myFxr%lCLD) cycle
    myFxr%FnAdj(ig) = myFxr%FnAdj(ig) / Nsum
  ENDDO ! DO j = 1, nLocalFxr
ENDDO ! DO ipin = 1, nxy
!$OMP END DO
!$OMP END PARALLEL
deallocate(NresPinAvg)
NULLIFY(PIN,ResVarPin,CellInfo,myFxr,isodata,pnum,idiso)
END SUBROUTINE
SUBROUTINE UpdtFtAdj(Core, Fxr, ilv, ig, iz, PE)
USE Typedef,      ONLY : CoreInfo_Type, FxrInfo_Type, Cell_TYpe, Pin_Type, PE_Type
USE XsLib_Mod,    ONLY : libdata, ldiso, mapnucl, mapnuclres
USE XsUtil_mod,   ONLY : LineIntPol, LineIntPol2
USE TH_Mod,       ONLY : GetPinFuelTemp
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(PE_TYPE) :: PE
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :), myFxr
TYPE(Cell_TYpe), POINTER :: CellInfo(:)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(libdata), POINTER :: isodata, jsodata

INTEGER,INTENT(IN) :: ilv, ig, iz
INTEGER :: nxy, ipin, FxrIdxSt, icel, nlocalFxr, j, ifxr
INTEGER :: niso, iso, jso, jd, id, idres, jdres, nThread, npot, nRiTemp, i
INTEGER, POINTER :: idiso(:)
REAL :: ind, siglp, TempAvgsq, Tempsq, micsigb, micsig0, sigbsq, xdat(20), ydat(20), Nreg, Navg, I_avg, I_reg, var(200)
REAL, POINTER :: pnum(:)

nxy = Core%nxy
Pin => Core%Pin; CellInfo => Core%CellInfo

!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(ipin, FxrIdxSt, icel, nlocalFxr, TempAvgsq, j, ifxr, myFxr, Tempsq, niso, idiso, pnum, siglp, jso, jd, &
!$OMP         jsodata, iso, id, isodata, idres, jdres, ind, micsigb, micsig0, sigbsq, npot, nRiTemp, &
!$OMP         i, xdat, ydat, I_avg, I_reg, Navg, Nreg, var )
!$OMP DO
DO ipin = 1, nxy
  FxrIdxSt = Pin(ipin)%FxrIdxSt
  icel = Pin(ipin)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
  if (.not.CellInfo(icel)%lres) CYCLE
  if (CellInfo(icel)%lAIC) CYCLE
  TempAvgsq = dsqrt(GetPinFuelTemp(Core, Fxr, iz, ipin))
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j - 1
    myFxr => Fxr(ifxr, iz)
    if (.not.myFxr%lres) cycle
    if (myFxr%lCLD) cycle
    if (myFxr%xseq_f_mg(ilv,ig).eq.0._8) cycle

    Tempsq = dsqrt(myFxr%temp)
    niso = myFxr%niso; idiso => myFxr%idiso; pnum => myFxr%pnum

    siglp = 0._8
    DO jso = 1, niso
      jd = mapnucl(idiso(jso));   jsodata => ldiso(jd)
      siglp = siglp + pnum(jso) * jsodata%lamsigp(ig)
    ENDDO

    Nreg = 0._8; Navg = 0._8
    DO iso = 1, niso
      id = mapnucl(idiso(iso));   isodata => ldiso(id)
      if (.NOT.isodata%lreso) CYCLE
      idres = mapnuclres(idiso(iso),iz)
      ind = 0._8
      DO jso = 1, niso
        jd = mapnucl(idiso(jso));   jsodata => ldiso(jd)
        jdres = mapnuclres(idiso(jso),iz)
        if (jdres.eq.idres) ind = ind + pnum(jso)
      ENDDO

      micsigb = (myFxr%XsEq_f_mg(ilv,ig) + siglp) / ind / myFxr%FtAdj(ilv,ig)
      micsig0 = micsigb - isodata%lamsigp(ig)
      if (micsig0.le.0._8) micsig0 = 1E-10
      sigbsq = dsqrt(micsig0)
      npot = isodata%nsig0
      nRiTemp = isodata%nrtemp
      DO i = 1, nRiTemp
        xdat(i) = isodata%rtempsq(i)
        var(1:npot) = isodata%ri_a(1:npot, ig, i)
        ydat(i) = LineIntPol2(sigbsq ,npot, isodata%sig0sq(1:npot), var(1:npot))
      ENDDO
      I_reg = LineIntpol(tempsq, nRiTemp, xdat(1:nRiTemp), ydat(1:nRiTemp))
      Nreg = Nreg + pnum(iso) * I_reg
      I_avg = LineIntpol(TempAvgsq, nRiTemp, xdat(1:nRiTemp), ydat(1:nRiTemp))
      Navg = Navg + pnum(iso) * I_avg

    ENDDO ! DO iso = 1, niso

    myFxr%FtAdj(ilv,ig) = Nreg / Navg

  ENDDO ! DO j = 1, nLocalFxr
ENDDO ! DO ipin = 1, nxy
!$OMP END DO
!$OMP END PARALLEL
NULLIFY(PIN,CellInfo,myFxr,isodata,jsodata,pnum,idiso)
END SUBROUTINE

SUBROUTINE UpdtCoreIsoInfo(Core, Fxr, PE)
USE TYPEDEF,      ONLY : CoreInfo_Type,    FxrInfo_Type,   Pin_TYPE, ResVarPin_TYPE, Cell_TYPE, PE_Type
USE XSLIB_MOD,    ONLY : nelthel, CoreResIsoUpdate, mapnuclres
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type),POINTER :: Fxr(:, :), myFxr
TYPE(PE_TYPE) :: PE
TYPE(Pin_TYPE), POINTER :: PIN(:)
TYPE(ResVarPin_TYPE), POINTER :: ResVarPIN(:,:)
TYPE(Cell_TYPE), POINTER :: CellInfo(:)

INTEGER :: nxy, myzb, myze, xyb, xye, iz, ipin, FxrIdxSt, icel, nlocalFxr, j, ifxr
INTEGER :: niso, nisotot, iso, nid, idx, idisotot(nelthel), jso, idres, jdres
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
      IF (.NOT.ResVarPin(ipin,iz)%lresA) CYCLE
      idisotot = 0; pnumtot = 0._8; nisotot = 0; volsum = 0._8
      DO j =1, nLocalFxr
        ifxr = FxrIdxSt + j -1
        myFxr => Fxr(ifxr, iz)
        IF( .NOT. myFxr%lres ) CYCLE
        niso = myFxr%niso
        pnum => myFxr%pnum; idiso => myFxr%idiso
        CALL CoreResIsoUpdate(idiso(1:niso), niso, iz)
        IF ( myFxr%lCLD ) CYCLE   ! fuel or AIC

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
        ENDDO
        volsum = volsum + vol
      ENDDO
      DO iso = 1, nisotot
        pnumtot(iso) = pnumtot(iso) / volsum
      ENDDO
      ResVarPin(ipin,iz)%niso = nisotot
      ResVarPin(ipin,iz)%idiso(1:nisotot)=idisotot(1:nisotot)
      ResVarPin(ipin,iz)%pnum(1:nisotot)=pnumtot(1:nisotot)
  ENDDO
ENDDO


DO iz = myzb, myze
  DO ipin = xyb, xye
      icel = Pin(ipin)%Cell(iz)
      nlocalFxr = CellInfo(icel)%nFxr
      FxrIdxSt = Pin(ipin)%FxrIdxSt
      IF (.NOT.ResVarPin(ipin,iz)%lresA) CYCLE
      idisotot = 0; pnumtot = 0._8; nisotot = 0; volsum = 0._8
      DO j =1, nLocalFxr
        ifxr = FxrIdxSt + j -1
        myFxr => Fxr(ifxr, iz)
        IF( .NOT. myFxr%lres ) CYCLE
!        IF (myFxr%lGD) print*, ifxr
        niso = myFxr%niso
        pnum => myFxr%pnum; idiso => myFxr%idiso
        myFxr%niso_Res = 0
        myFXR%pnum_Res = 0.
        myFXR%idx_Res = 0
        myFXR%idiso_Res = 0
        DO iso = 1, niso
          idres = mapnuclres(idiso(iso),iz);
!          IF (ifxr.EQ.152) WRITE(*, "(3I ES10.3)") iso, idres, idiso(iso), pnum(iso)
          IF (idres.EQ.0) CYCLE
          DO jso = 1, myFxr%niso_Res
            jdres = mapnuclres(myFXR%idiso_Res(jso) ,iz)
            IF (idres .EQ. jdres) EXIT
          END DO
          IF (jso .LE. myFxr%niso_Res .AND. myFxr%niso_Res .GT. 0) THEN
            myFXR%pnum_Res(jso) = myFXR%pnum_Res(jso) + pnum(iso);
            myFXR%idx_Res(iso) = jso
!            IF (ifxr .EQ. 152) print*, "BOOM!"
          ELSE
            myFxr%niso_Res = myFxr%niso_Res + 1
            myFXR%pnum_Res(myFxr%niso_Res) = pnum(iso);
            myFXR%idx_Res(iso) = myFxr%niso_Res
            myFXR%idiso_Res(myFxr%niso_Res) = idiso(iso)
!            IF (ifxr .EQ. 152) print*, "YUUM!"
          END IF
!          IF (ifxr.EQ.152) WRITE(*, "(3I ES10.2)") myFxr%idx_Res(iso), myFxr%niso_Res, jso, myFxr%pnum_Res(min(jso,myFxr%niso_Res))
        END DO
      ENDDO
      nisotot = ResVarPIN(ipin,iz)%niso

      IF (.NOT. ResVarPin(ipin,iz)%lres) CYCLE
      ResVarPin(ipin,iz)%niso_Res = 0
      ResVarPin(ipin,iz)%pnum_Res = 0
      ResVarPin(ipin,iz)%idx_Res = 0
      ResVarPin(ipin,iz)%idiso_Res = 0
      DO iso = 1, nisotot
        idres = mapnuclres(ResVarPin(ipin,iz)%idiso(iso),iz);
!        IF (ipin.EQ.13) WRITE(*, "(3I ES10.3)") iso, idres, ResVarPin(ipin,iz)%idiso(iso), ResVarPin(ipin,iz)%pnum(iso)
        IF (idres.EQ.0) CYCLE
        DO jso = 1, ResVarPin(ipin,iz)%niso_Res
          jdres = mapnuclres(ResVarPin(ipin,iz)%idiso_Res(jso),iz);
          IF (idres .EQ. jdres) EXIT
        END DO
        IF (jso .LE. ResVarPin(ipin,iz)%niso_Res .AND. ResVarPin(ipin,iz)%niso_Res .GT. 0) THEN
          ResVarPin(ipin,iz)%pnum_Res(jso) = ResVarPin(ipin,iz)%pnum_Res(jso) + ResVarPin(ipin,iz)%pnum(iso);
          ResVarPin(ipin,iz)%idx_Res(iso) = jso
!          IF (ipin.EQ.13) print*, "BOOM!"
        ELSE
          ResVarPin(ipin,iz)%niso_Res = ResVarPin(ipin,iz)%niso_Res + 1
          ResVarPin(ipin,iz)%pnum_Res(ResVarPin(ipin,iz)%niso_Res) = ResVarPin(ipin,iz)%pnum(iso);
          ResVarPin(ipin,iz)%idx_Res(iso) = ResVarPin(ipin,iz)%niso_Res
          ResVarPin(ipin,iz)%idiso_Res(ResVarPin(ipin,iz)%niso_Res) = ResVarPin(ipin,iz)%idiso(iso)
!          IF (ipin.EQ.13) print*, "YUUM!"
        END IF
!        IF (ipin.EQ.13) WRITE(*, "(3I ES10.2)") ResVarPin(ipin,iz)%idx_Res(iso), ResVarPin(ipin,iz)%niso_Res, jso, ResVarPin(ipin,iz)%pnum_Res(min(jso,ResVarPin(ipin,iz)%niso_Res))
      END DO
  ENDDO
ENDDO

NULLIFY(Pin,ResVarPin,CellInfo,myFxr,idiso,pnum)
END SUBROUTINE
SUBROUTINE UpdtResoCat(PE)
USE TYPEDEF,      ONLY : PE_Type
USE XSLIB_MOD,    ONLY : nRes_Cor, idres_cor, nCat, ResoCat, mapnuclRes, nActiveCat, ResoCatUse
IMPLICIT NONE
TYPE(PE_TYPE) :: PE
INTEGER :: i, id, ic, iiso, iz, myzb, myze

myzb = PE%myzb; myze = PE%myze
DO iz = myzb, myze
    do ic=1,nCat
        ResoCatUse(ic,iz)=.false.
    enddo
    do i = 1, nRes_Cor(iz)
        id=idres_cor(i,iz)
        do ic = 1, nCat
            if (ResoCatUse(ic,iz)) cycle
            do iiso = 1, ResoCat(ic)%niso
                if (id.eq.ResoCat(ic)%idiso(iiso)) then
                    ResoCatUse(ic,iz)=.true.
                    exit
                endif
            enddo
        enddo
    enddo

    nActiveCat(iz) = 0
    do ic = 1, nCat
        if (.not.ResoCatUse(ic,iz)) cycle
        nActiveCat(iz) = nActiveCat(iz) + 1
        id=ResoCat(ic)%repid
        if (mapnuclRes(id,iz).eq.0) then
            nRes_Cor(iz) = nRes_Cor(iz) + 1
            idres_cor(nRes_Cor,iz) = id
            mapnuclRes(id,iz) = nRes_Cor(iz)
            mapnuclRes(id+500,iz) = nRes_Cor(iz)
        endif
    enddo
ENDDO
END SUBROUTINE
subroutine EquivRsearch(Core,icel,eR,niter,D,rad,sigt,q,nr,lb)
! subroutine to search equivalent 1D cell radius
USE TYPEDEF,  ONLY : CoreInfo_Type,ResVarPin_Type, Cell_TYPE
use CP_mod
implicit none
TYPE(CoreInfo_Type) :: Core
TYPE(Cell_TYPE), POINTER :: CellInfo(:)
integer,intent(in) :: nr,icel
real,intent(out) :: eR
real,intent(in) :: D,rad(nr),sigt(nr),q(nr),lb
integer,intent(out) :: niter
integer :: nfueldiv,cldifxr
real :: phi(nr),r,volexceptmod
real :: D1,D2,r1(nr),r2(nr),f1,f2,fp,err

CellInfo => Core%CellInfo
nfueldiv = CellInfo(icel)%nfueldiv
volexceptmod = CellInfo(icel)%fuelgapcldvol
cldifxr = CellInfo(icel)%cldfxridx

r1=rad; r2=rad
r=r1(nr-1)*1.01
call adjustgeom(r1,r,volexceptmod,cldifxr,nr)
call runCP_(phi,r1,sigt,q,nr); D1=(sum(phi(1:nfueldiv)*sigt(1:nfueldiv))/nfueldiv-q(1))*lb; f1=(D1/D-1._8)
call runCP_(phi,r2,sigt,q,nr); D2=(sum(phi(1:nfueldiv)*sigt(1:nfueldiv))/nfueldiv-q(1))*lb; f2=(D2/D-1._8)

fp=f1*f2
if (fp.eq.0._8) then
    if (f1.eq.0._8) then
        eR=r1(nr)
    else
        eR=r2(nr)
    endif
    return
else
    err=1._8; niter=2
    do while (err.gt.1E-10)
        if (fp.gt.0._8) then
            if (abs(f1).gt.abs(f2)) then
                r=r2(nr)*r2(nr)/r1(nr)
                call adjustgeom(r1,r,volexceptmod,cldifxr,nr)
                call runCP_(phi,r1,sigt,q,nr); D1=(sum(phi(1:nfueldiv)*sigt(1:nfueldiv))/nfueldiv-q(1))*lb; f1=(D1/D-1._8)
            else
                r=r1(nr)*r1(nr)/r2(nr)
                call adjustgeom(r2,r,volexceptmod,cldifxr,nr)
                call runCP_(phi,r2,sigt,q,nr); D2=(sum(phi(1:nfueldiv)*sigt(1:nfueldiv))/nfueldiv-q(1))*lb; f2=(D2/D-1._8)
            endif
        else
            if (abs(f1).lt.abs(f2)) then
                r=dsqrt(r1(nr)*r2(nr))
                call adjustgeom(r2,r,volexceptmod,cldifxr,nr)
                call runCP_(phi,r2,sigt,q,nr); D2=(sum(phi(1:nfueldiv)*sigt(1:nfueldiv))/nfueldiv-q(1))*lb; f2=(D2/D-1._8)
            else
                r=dsqrt(r1(nr)*r2(nr))
                call adjustgeom(r1,r,volexceptmod,cldifxr,nr)
                call runCP_(phi,r1,sigt,q,nr); D1=(sum(phi(1:nfueldiv)*sigt(1:nfueldiv))/nfueldiv-q(1))*lb; f1=(D1/D-1._8)
            endif
        endif
        fp=f1*f2
        err=dabs(fp)
        niter=niter+1
        if (niter.gt.1000) exit
    enddo
endif
eR=dsqrt(r1(nr)*r2(nr))

contains

subroutine adjustgeom(r1,r,volexceptmod,cldifxr,nr)
use param, only : pi
implicit none
integer,intent(in) :: nr,cldifxr
real,intent(in) :: r,volexceptmod
real,intent(inout) :: r1(nr)
integer :: ir
real :: subvol,volsum,ndiv

ndiv=nr-cldifxr
subvol=(pi*r**2-volexceptmod)/ndiv
volsum=volexceptmod
do ir = cldifxr+1,nr
    volsum = volsum + subvol
    r1(ir)=dsqrt(volsum/pi)
enddo
end subroutine

END SUBROUTINE
