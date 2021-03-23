#include <defines.h>
SUBROUTINE SubGrpEffXsGen(Core, Fxr, THInfo, eigv, GroupInfo, nTracerCntl, PE)
!Effective Xs Generation
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,       FxrInfo_Type,          GroupInfo_Type,     &
                           PE_Type,             THInfo_Type,           pin_type,           &
                           ResVarpin_type,      cell_type,             XsMac_Type
USE CNTL,           ONLY : nTracerCntl_Type
USE MacXsLib_mod,   ONLY : EffMacXS,   EffRIFPin,       MacXsBase
USE BasicOperation, ONLY : CP_VA
USE IOUTIL,         ONLY : message
USE FILES,          ONLY : io8
USE TH_Mod,         ONLY : GetPinFuelTemp
USE SUBGRP_MOD,     ONLY : UpdtCoreIsoInfo
USE MPIComm_Mod,    ONLY : MPI_SYNC
USE TIMER,          ONLY : nTracer_dclock,       TimeChk
USE OMP_LIB
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(GroupInfo_TYPE) :: GroupInfo
TYPE(THInfo_Type) :: THInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
TYPE(Pin_Type), POINTER :: PIN(:)
TYPE(ResVarPin_Type), POINTER :: ResVarPin(:,:),mypin
TYPE(Cell_Type), POINTER :: CellInfo(:), myCell
TYPE(FxrInfo_Type), POINTER :: myFxr
TYPE(XsMac_Type), SAVE :: XsMac(nThreadMax)

REAL,INTENT(IN) :: eigv
INTEGER :: nxy, myzb, myze, iz, ipin, xyb, xye, FxrIdxSt, nlocalFxr, j, ifxr, icel  !--- (xyb, xye)CNJ Edit : Domain Decomposition + MPI
INTEGER :: ng, ig, iresoGrp1, iresoGrp2, niso, iso, tid
REAL :: PinFuelTempAvgsq, XsMacAold, XsMacFold, XsMacNFold, XsMackFold, XsMacSold, XsMacStrold
REAL :: isoXsMacAold(700),isoXsMacFold(700),isoXsMacNFold(700),isoXsMackFold(700),isoXsMacSold(700),isoXsMacS1old(700),isoXsMacSSold(700),isoxsmaccapold(700)
LOGICAL :: lxslib, lSubGrpSweep, lAIC, lmcxs
REAL :: Tbeg, Tend

INTEGER :: iso_mc,niso_mc,imc,iz_mc,ixa_mc,iya_mc,ixya_mc,iasytype_mc,ixc_mc,iyc_mc,ixyc_mc,ixy_mc,ifxr1_mc,ifxr2_mc

IF (.not.nTracerCntl%lrestrmt) RETURN
lxslib = nTracerCntl%lxslib
lSubGrpSweep = nTracerCntl%lSubGrpSweep
IF(.NOT. lxslib) RETURN

WRITE(mesg, '(A)') 'Update Subgroup Effective XSs...'
IF(PE%Master) CALL message(io8, TRUE, TRUE, mesg)
Tbeg = nTracer_Dclock(FALSE, FALSE)

Pin => Core%Pin; CellInfo => Core%CellInfo
ResVarPin => Core%ResVarPin;
nxy = Core%nxy;  myzb = PE%myzb; myze = PE%myze

!--- CNJ Edit : Domain Decomposition + MPI
xyb = PE%myPinBeg; xye = PE%myPinEnd
IF (PE%RTMASTER) THEN
  xyb = 1; xye = nxy
ENDIF

ng = GroupInfo%ng
iresoGrp1 = GroupInfo%nofg + 1; iresoGrp2 = GroupInfo%nofg + GroupInfo%norg

CALL UpdtCoreIsoInfo(Core, Fxr, PE)
if (.not.nTracerCntl%lMLG) then
    if (nTracerCntl%lCAT) CALL UpdtResoCat(PE)
endif

CALL omp_set_num_threads(PE%nThread)

DO iz = myzb, myze
  IF (.NOT. Core%lFuelPlane(iz)) CYCLE
  !! The routine 'EffMacXs' is not parallelized.
  !! Instead, pin-wise parallelization is adopted. (2018-03-24 by PHS)
  !$OMP PARALLEL DEFAULT(SHARED)      &
  !$OMP PRIVATE(ipin, mypin, FxrIdxSt, icel, myCell, lAIC, nlocalFxr, PinFuelTempAvgsq, j, ifxr, myFxr, niso, &
  !$OMP         tid, ig, XsMacAold, XsMacSold, XsMacFold, XsMacNFold, XsMackFold, XsMacStrold, isoXsMacFold, isoXsMacAold, isoXsMacSold, isoXsMacNFold, isoXsMackFold, &
  !$OMP         isoXsMacS1old, isoXsMacSSold, IsoXsMacCapOld, iso)
  !$  tid = omp_get_thread_num()+1
  !$OMP DO
  DO ipin = xyb, xye
    mypin => ResVarPin(ipin,iz)
    FxrIdxSt = Pin(ipin)%FxrIdxSt
    icel = Pin(ipin)%Cell(iz)
    myCell => CellInfo(icel)
    lAIC = myCell%lAIC
    nlocalFxr = CellInfo(icel)%nFxr

    lmcxs=.false.
    IF (nTracerCntl%lMCXS) THEN
      do imc=1,nTracerCntl%nmcxs
        iz_mc = nTracerCntl%mcxs(imc)%iz
        ixa_mc = nTracerCntl%mcxs(imc)%iasyx; iya_mc = nTracerCntl%mcxs(imc)%iasyy
        ixya_mc = Core%CoreIdx(ixa_mc, iya_mc); iasytype_mc = Core%CoreMap(ixya_mc)
        ixc_mc = nTracerCntl%mcxs(imc)%ix; iyc_mc = nTracerCntl%mcxs(imc)%iy
        ixyc_mc = Core%AsyInfo(iasyType_mc)%pin2DIdx(ixc_mc, iyc_mc)
        ixy_mc = Core%Asy(ixya_mc)%GlobalPinIdx(ixyc_mc)
        if (ipin.eq.ixy_mc) exit
      enddo
      if (imc.le.nTracerCntl%nmcxs) then
          lmcxs=.true.
          niso_mc=nTracerCntl%mcxs(imc)%niso
          ifxr1_mc=nTracerCntl%mcxs(imc)%ifxr1
          ifxr2_mc=nTracerCntl%mcxs(imc)%ifxr2
      endif
    ENDIF

    PinFuelTempAvgsq = dsqrt(GetPinFuelTemp(Core, Fxr, iz, ipin))

    IF (mypin%lres.AND.nTracerCntl%lRIF.AND..NOT.nTRACERCntl%lRIFFXR) CALL EffRIFPin(mypin, PinFuelTempAvgsq, iz, lAIC, PE)

    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j -1
      myFxr => Fxr(ifxr, iz)
      IF( .NOT. myFxr%lres ) CYCLE

      niso = myFxr%niso

      do ig = iResoGrp1, iResoGrp2
        myFxr%fresoa(ig) = 1._8
        myFxr%fresof(ig) = 1._8
        myFxr%fresoNf(ig) = 1._8
        myFxr%fresokf(ig) = 1._8
        do iso = 1, niso
          myFXR%fresoAIso(iso,ig) = 1._8
          myFXR%fresoFIso(iso,ig) = 1._8
        enddo
      enddo
      IF (nTracerCntl%lRST) THEN
        do ig = iResoGrp1, iResoGrp2
          myFxr%fresos(ig) = 1._8
          myFxr%fresostr(ig) = 1._8
          do iso = 1, niso
            myFXR%fresoSIso(iso,ig) = 1._8
            myFXR%fresoSSIso(iso,ig) = 1._8
            myFXR%fresoS1Iso(iso,ig) = 1._8
          enddo
        enddo
      ENDIF
      IF(.NOT. lSubGrpSweep) CYCLE

      CALL MacXsBase(XsMac(tid), myFxr, iResoGrp1, iResoGrp2, ng, eigv, FALSE, TRUE)

      DO ig = iResoGrp1, iResoGrp2
        XsMacAold = XsMac(tid)%XsMacA(ig); XsMacFold = XsMac(tid)%XsMacF(ig);
        XsMacNFold = XsMac(tid)%XsMacNF(ig); XsMackFold = XsMac(tid)%XsMackF(ig)
        XsMacSold = XsMac(tid)%XsMacS(ig); XsMacStrold = XsMac(tid)%XsMacStr(ig);
        CALL CP_VA(isoXsMacFold(1:niso),XsMac(tid)%IsoXsMacF(1:niso,ig),niso)
        CALL CP_VA(isoXsMacNFold(1:niso),XsMac(tid)%IsoXsMacNF(1:niso,ig),niso)
        CALL CP_VA(isoXsMackFold(1:niso),XsMac(tid)%IsoXsMackF(1:niso,ig),niso)
        CALL CP_VA(isoXsMacAold(1:niso),XsMac(tid)%IsoXsMacA(1:niso,ig),niso)
        IF (nTracerCntl%lRST) THEN
          CALL CP_VA(isoXsMacSold(1:niso),XsMac(tid)%IsoXsMacS0(1:niso,ig),niso)
          CALL CP_VA(isoXsMacS1old(1:niso),XsMac(tid)%IsoXsMacS1(1:niso,ig),niso)
          CALL CP_VA(isoXsMacSSold(1:niso),XsMac(tid)%IsoXsMacSS(1:niso,ig),niso)
        ENDIF
        CALL CP_VA(IsoXsMacCapOld(1:niso),XsMac(tid)%IsoXsRadCap(1:niso,ig),niso)     !-- JSU EDIT 20170727
        CALL EffMacXs(XsMac(tid), mypin, myFxr, PinFuelTempAvgsq, niso, ig, ng, .TRUE., iz, PE, ifxr)

        if (lmcxs.and.((j.ge.ifxr1_mc).and.(j.le.ifxr2_mc))) then
          do iso=1,niso
            do iso_mc=1,niso_mc
              if (nTracerCntl%mcxs(imc)%idiso(iso_mc).eq.myFxr%idiso(iso)) exit
            enddo
            if (iso_mc.gt.niso_mc) cycle
            XsMac(tid)%XsMacA(ig) = XsMac(tid)%XsMacA(ig) - XsMac(tid)%IsoXsMacA(iso,ig) + nTracerCntl%mcxs(imc)%isoxsa(ig,j,iso_mc)
            XsMac(tid)%IsoXsMacA(iso,ig) = nTracerCntl%mcxs(imc)%isoxsa(ig,j,iso_mc)
            XsMac(tid)%XsMacNF(ig) = XsMac(tid)%XsMacNF(ig) - XsMac(tid)%IsoXsMacNF(iso,ig) + nTracerCntl%mcxs(imc)%isoxsNF(ig,j,iso_mc)
            XsMac(tid)%IsoXsMacNF(iso,ig)= nTracerCntl%mcxs(imc)%isoxsnf(ig,j,iso_mc)
          enddo
          myFxr%fresoa(ig) = XsMac(tid)%XsMacA(ig) / XsMacAold
          IF(XsMacNFold .gt. epsm8) THEN
            myFxr%fresof(ig) = XsMac(tid)%XsMacf(ig) / XsMacFold
            myFxr%fresoNf(ig) = XsMac(tid)%XsMacNf(ig) / XsMacNFold
            myFxr%fresokf(ig) = XsMac(tid)%XsMackf(ig) / XsMackFold
          END IF
          DO iso = 1, niso
            myFXR%fresoAIso(iso,ig) = XsMac(tid)%IsoXsMacA(iso,ig) / isoXsMacAold(iso)
            IF (isoXsMacNFold(iso) .GT. epsm8) myFXR%fresoFIso(iso,ig) = XsMac(tid)%IsoXsMacNF(iso,ig) / isoXsMacNFold(iso)
          ENDDO
          IF (nTRACERCntl%lGamma) THEN
            IF(XsMacNFold .gt. epsm8) THEN ! Radioactive Capture resonance fraction...
              DO iso = 1, niso
                myFxr%fresocapIso(iso, ig) = (XsMac(tid)%IsoXsMacA(iso,ig) - XsMac(tid)%IsoXsMacF(iso,ig)) / (isoXsMacAold(iso) - isoXsMacFold(iso))
              END DO
            ELSE
              DO iso = 1, niso
                myFxr%fresocapIso(iso, ig) = XsMac(tid)%IsoXsMacA(iso,ig) / isoXsMacAold(iso)
              END DO
            END IF
          END IF
        else
          myFxr%fresoa(ig) = XsMac(tid)%XsMacA(ig) / XsMacAold
          IF(XsMacFold .gt. epsm8) THEN
            myFxr%fresof(ig) = XsMac(tid)%XsMacf(ig) / XsMacFold
            myFxr%fresoNf(ig) = XsMac(tid)%XsMacNf(ig) / XsMacNFold
            myFxr%fresokf(ig) = XsMac(tid)%XsMackf(ig) / XsMackFold
!            myFxr%fresokf(ig) = myFxr%fresof(ig)
          END IF

          DO iso = 1, niso
            myFXR%fresoAIso(iso,ig) = XsMac(tid)%IsoXsMacA(iso,ig) / isoXsMacAold(iso)
            IF (isoXsMacFold(iso) .GT. epsm8) myFXR%fresoFIso(iso,ig) = XsMac(tid)%IsoXsMacF(iso,ig) / isoXsMacFold(iso)
          ENDDO
          IF (nTRACERCntl%lGamma) THEN
          IF(isoXsMacFold(iso) .gt. epsm8) THEN ! Radioactive Capture resonance fraction...
            DO iso = 1, niso
              myFxr%fresocapIso(iso, ig) = (XsMac(tid)%IsoXsMacA(iso,ig) - XsMac(tid)%IsoXsMacF(iso,ig)) / (isoXsMacAold(iso) - isoXsMacFold(iso))
            END DO
          ELSE
            DO iso = 1, niso
              myFxr%fresocapIso(iso, ig) = XsMac(tid)%IsoXsMacA(iso,ig) / isoXsMacAold(iso)
            END DO
          END IF
          END IF
        endif


        IF (nTracerCntl%lRST) THEN
          myFxr%fresos(ig) = XsMac(tid)%XsMacS(ig) / XsMacSold
          myFxr%fresostr(ig) = XsMac(tid)%XsMacStr(ig) / XsMacStrold
          DO iso = 1, niso
            if (abs(isoXsMacSold(iso)).GT. epsm8) myFXR%fresoSIso(iso,ig) = XsMac(tid)%IsoXsMacS0(iso,ig)/isoXsMacSold(iso)
            if (abs(isoXsMacS1old(iso)).GT. epsm8) myFXR%fresoS1Iso(iso,ig) = XsMac(tid)%IsoXsMacS1(iso,ig)/isoXsMacS1old(iso)
            if (abs(isoXsMacSSold(iso)).GT. epsm8) myFXR%fresoSSIso(iso,ig) = XsMac(tid)%IsoXsMacSS(iso,ig)/isoXsMacSSold(iso)
          ENDDO
        ENDIF

      ENDDO

    ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
ENDDO

CALL MPI_SYNC(PE%MPI_NTRACER_COMM)

NULLIFY(PIN, ResVarPin, CellINfo, myFxr, myPin, myCell)
Tend = nTracer_Dclock(FALSE, FALSE)
TimeChk%SubGrpGenEFFXSTime = TimeChk%SubGrpGenEFFXSTime + (Tend - Tbeg)

END SUBROUTINE

SUBROUTINE FxrChiGen(Core, Fxr, FmInfo, GroupInfo, PE, nTracerCntl)
USE TYPEDEF,        ONLY :   CoreInfo_Type,     FxrInfo_Type,  FmInfo_Type,            &
                             GroupInfo_Type,    PE_Type,                               &
                             Cell_Type,         Pin_Type
USE CNTL,           ONLY : nTracerCntl_type
USE MacXsLib_mod,   ONLY : GetMacChi
USE BasicOperation, ONLY : CP_CA
IMPLICIT NONE
TYPE(CoreInfo_Type),INTENT(IN) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(FmInfo_Type),INTENT(IN) :: FmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(PE_Type) :: PE
TYPE(nTracerCntl_Type) :: nTracerCntl

!Pointing Variables
TYPE(FxrInfo_Type), POINTER :: myFxr
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(Pin_Type), POINTER :: Pin(:)
REAL, POINTER :: phis(:, :, :)

INTEGER :: myzb, myze, nxy, nfxr, nfsr, ng, nchi
INTEGER :: nlocalFxr, nFsrInFxr, FsrIdxSt, FxrIdxSt
INTEGER :: ipin ,icel, iz, ifsr, ifxr, ig
INTEGER :: i, j, k, l

REAL, POINTER :: Spectrum(:)
REAL :: vol, volsum, phisum
LOGICAL :: lxslib

lxsLib = nTracerCntl%lxslib
IF(.NOT. lxsLib) RETURN

CellInfo => Core%CellInfo; Pin => Core%Pin
phis => FmInfo%phis

myzb = PE%myzb; myze = PE%myze
ng = GroupInfo%ng; nchi = GroupInfo%nchi
nxy = COre%nxy; nFsr = Core%nCoreFsr; nFxr = Core%nCoreFxr
!ALLOCATE(Spectrum(ng))
DO iz = myzb, myze
  !$OMP PARALLEL DEFAULT(SHARED)&
  !$OMP PRIVATE(ipin,FxrIdxSt,FsrIdxSt,icel,nlocalFxr,j,ifxr,myFxr,nFsrInFxr,ig,i,l,vol,ifsr,volsum,phisum,Spectrum)
  ALLOCATE(Spectrum(ng))
  !$OMP DO SCHEDULE(GUIDED)
  DO ipin = 1, nxy
    FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
    icel = Pin(ipin)%Cell(iz)
    nlocalFxr = CellInfo(icel)%nFxr

    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j -1;  myFxr => Fxr(ifxr, iz)
      IF(.NOT. myFXR%lfUEL) CYCLE
      nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
      !CALL CP_CA(Spectrum, 0._8, ng)
      Spectrum = 0._8;
      !Obtain a spectrum of Flat Xs Region
      DO ig = 1, ng
        volsum = 0; phisum = 0
        DO i = 1, nFsrInFxr
          l = Cellinfo(icel)%MapFxr2FsrIdx(i, j); vol = CellInfo(icel)%vol(l)
          ifsr = FsrIdxSt + l - 1; volsum = volsum + vol
          phisum = phisum + phis(ifsr, iz, ig)* vol
        ENDDO
        Spectrum(ig) = phisum/volsum
      ENDDO

      CALL GetMacChi(myFXR, Spectrum, 1, nchi, nchi, ng, .FALSE.)
    ENDDO
  ENDDO
  !$OMP END DO
  DEALLOCATE(Spectrum)
  !$OMP END PARALLEL
ENDDO

!Free the variables
!Deallocate(Spectrum)
IF(ASSOCIATED(CellInfo)) NULLIFY(CellInfo)
IF(ASSOCIATED(Pin)) NULLIFY(Pin)
IF(ASSOCIATED(Phis)) NULLIFY(Phis)
!NULLIFY(MyFxr)

END SUBROUTINE
