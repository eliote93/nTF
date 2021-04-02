#include <defines.h>
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetRtMacXsNM(core, Fxr, xstnm, iz, ng, lxslib, lTrCorrection, lRST, lssph, lssphreg, PE)

USE PARAM
USE OMP_LIB
USE TYPEDEF,      ONLY : coreinfo_type, Fxrinfo_type, Cell_Type, pin_Type, PE_TYPE, XsMac_Type
USE Core_mod,     ONLY : GroupInfo
USE BenchXs,      ONLY : GetXstrBen, GetXsTrDynBen, GetXsTBen, GetXsTDynBen
USE MacXsLib_mod, ONLY : MacXsBase,  BaseMacSTr
USE SPH_mod,      ONLY : ssphfnm
USE XSLIB_MOD,    ONLY : igresb,igrese
USE TRAN_MOD,     ONLY : TranInfo, TranCntl

IMPLICIT NONE

TYPE(coreinfo_type) :: Core
TYPE(Fxrinfo_type), DIMENSION(:) :: Fxr
TYPE(PE_TYPE) :: PE

INTEGER :: iz, ng
LOGICAL :: lxslib, lTrCorrection, lRST, lssph, lssphreg
! ----------------------------------------------------
TYPE(Pin_Type),  POINTER, DIMENSION(:) :: Pin
TYPE(Cell_Type), POINTER, DIMENSION(:) :: CellInfo

INTEGER :: nCoreFxr, nCoreFsr, nxy, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr, nofg, norg, xyb, xye
INTEGER :: icel, ipin, ifxr, ifsr, ifsrlocal, itype, ig, tid
INTEGER :: i, j, k

REAL :: macstr

REAL, DIMENSION(ng)     :: xsmactr
REAL, DIMENSION(ng, 40) :: SPHfac

LOGICAL :: lresogrp(ng), lres, lress

REAL, POINTER, DIMENSION(:,:) :: xstnm

TYPE(XsMac_Type), POINTER, DIMENSION(:) :: XsMac
! ----------------------------------------------------

Pin      => Core%Pin
CellInfo => Core%CellInfo
nCoreFsr  = Core%nCoreFsr
nCoreFxr  = Core%nCoreFxr
nxy       = Core%nxy

xyb = PE%myPinBeg
xye = PE%myPinEnd ! Domain Dcmp + MPI

ALLOCATE (XsMac (PE%nThread))

IF(lxslib) THEN
  nofg = GroupInfo%nofg
  norg = GroupInfo%norg
  
  DO ig = 1, ng
    lresogrp(ig) = FALSE
    IF (ig.GT.nofg .AND. ig.LE.(nofg + norg)) lresogrp(ig) = TRUE
  ENDDO
ENDIF

SPHfac = ONE
! ----------------------------------------------------
!$OMP PARALLEL DEFAULT(SHARED)                                                                                      &
!$OMP PRIVATE(i, j, k, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr,                                                    &
!$OMP         icel, ipin, ifxr, ifsr, ifsrlocal, itype, ig, tid, xsmactr, lres, lress, SPHfac)
tid = omp_get_thread_num() + 1
!$OMP DO SCHEDULE(DYNAMIC)
DO ipin = xyb, xye
  FsrIdxSt  = Pin(ipin)%FsrIdxSt
  FxrIdxSt  = Pin(ipin)%FxrIdxSt
  icel      = Pin(ipin)%Cell(iz)
  nlocalFxr = CellInfo(icel)%nFxr 
  
  DO j = 1, nLocalFxr
    ifxr      = FxrIdxSt + j -1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    ! ------------------------------------------------
    IF (lxslib) THEN
      CALL MacXsBase(XSMac(tid), Fxr(ifxr), 1, ng, ng, ONE, FALSE, TRUE, TRUE)
      
      DO ig = 1, ng
        lres  = lresogrp(ig) .AND. Fxr(ifxr)%lres
        lress = lres .AND. lRST
        
        IF(lres) XsMac(tid)%XsMacA(ig) = XsMac(tid)%XsMacA(ig) * Fxr(ifxr)%fresoa(ig)
        
        IF(lTrCorrection) THEN
          CALL BaseMacSTr(XsMac(tid), Fxr(ifxr), ig, ng, TRUE)  
          
          IF(lress) XsMac(tid)%XsMacStr(ig) = XsMac(tid)%XsMacStr(ig) * Fxr(ifxr)%fresostr(ig)  
          
          xsmactr(ig) = XsMac(tid)%XsMacA(ig) + XsMac(tid)%XsMacstr(ig)
        ELSE
          IF (lress) XsMac(tid)%XsMacS(ig) = XsMac(tid)%XsMacS(ig) * Fxr(ifxr)%fresos(ig)
          
          xsmactr(ig) = XsMac(tid)%XsMacA(ig) + XsMac(tid)%XsMacs(ig)
        END IF
#ifdef inflow
        xsmactr(ig) = xsmactr(ig) + fxr(ifxr)%Delinflow(ig)
#endif
      END DO
      
      IF (lsSPH .AND. CellInfo(icel)%lsSPH) THEN
        DO ig = igresb, igrese
          IF (lsSPHreg) THEN
            SPHfac(ig, j) = Fxr(ifxr)%SPHfactor(ig)
          ELSE
            SPHfac(ig, j) = CellInfo(icel)%SPHfactor(j, ig)
          END IF
          
          xsmactr(ig) = xsmactr(ig) * SPHfac(ig, j)
        END DO
      END IF
    ! ----------------------------------------------------
    ELSE
      ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(1,j)
      itype     = Fxr(ifxr)%imix
      
      IF (lTrCorrection) THEN
        IF (TranCntl%lDynamicBen) THEN
          xsmactr = GetXsTrDynBen(itype, TranInfo%fuelTemp(ipin, iz), 1, ng)
        ELSE
          xsmactr = GetXsTrBen(itype, 1, ng)
        END IF
      ELSE
        IF (TranCntl%lDynamicBen) THEN
          xsmactr = GetXsTDynBen(itype, TranInfo%fuelTemp(ipin, iz), 1, ng)
        ELSE
          xsmactr = GetXsTBen(itype, 1, ng)
        END IF
      END IF
    END IF
    ! ----------------------------------------------------
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      
      xstnm(:, ifsr) = xsmactr
      
      IF (lsSPH .AND. CellInfo(icel)%lsSPH) ssphfnm(igresb:igrese, ifsr, iz) = SPHfac(igresb:igrese, j)
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL
! ----------------------------------------------------
DO i = 1, PE%nThread
  IF (XsMac(i)%lIsoAlloc) THEN
    DEALLOCATE (XsMac(i)%IsoXsMacA,  XsMac(i)%IsoXsMacf,  XsMac(i)%ISoXsMacnf, XsMac(i)%IsoXsMacSS)
    DEALLOCATE (XsMac(i)%IsoXsMacS1, XsMac(i)%IsoXsMacS0, XsMac(i)%IsoXsMackf)
    DEALLOCATE (XsMac(i)%IsoXsMacT,  XsMac(i)%IsoXsMacTr, XsMac(i)%IsoXsRadCap)
  END IF
  
  IF (XsMac(i)%lalloc) THEN
    DEALLOCATE (XsMac(i)%xsmaca,    XsMac(i)%xsmacf,    XsMac(i)%xsmacnf, XsMac(i)%xsmackf)
    DEALLOCATE (XsMac(i)%xsmacs,    XsMac(i)%xsmacsm,   XsMac(i)%xsmacstr)
    DEALLOCATE (XsMac(i)%xsmacp1sm, XsMac(i)%xsmacp2sm, XsMac(i)%xsmacp3sm)
    DEALLOCATE (XsMac(i)%xsmact,    XsMac(i)%xsmactr,   XsMac(i)%Chi)
  END IF 
END DO

DEALLOCATE (XsMac)
! ----------------------------------------------------
NULLIFY (Pin)
NULLIFY (CellInfo)
! ----------------------------------------------------

END SUBROUTINE SetRtMacXsNM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE AddBucklingNM(Core, Fxr, xstnm, bsq, iz, ng, lxslib, lRST)

USE PARAM
USE TYPEDEF,      ONLY : coreinfo_type,       Cell_Type,      pin_Type,      &
                         Fxrinfo_type,        XsMac_Type
USE Core_mod,     ONLY : GroupInfo
USE BenchXs,      ONLY : GetXstrBen,          GetXstrDynBen
USE MacXsLib_mod, ONLY : MacXsBase,           BaseMacSTr
USE PE_MOD,       ONLY : PE
USE TRAN_MOD,     ONLY : TranInfo,            TranCntl
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(Fxrinfo_type), POINTER :: Fxr(:,:)
REAL, POINTER :: xstnm(:, :)
INTEGER :: iz, ng
REAL :: Bsq
LOGICAL :: lxslib, lRST

TYPE(FxrInfo_Type), POINTER, SAVE :: myFxr
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(XsMac_Type), SAVE :: XsMac 
REAL :: XsD
REAL :: xsmactr(ng)
INTEGER :: nCoreFsr, nCoreFxr, nxy
INTEGER :: FsrIdxSt, ipin, icel, ifsr, ifxr, itype, ifsrlocal, ig
INTEGER :: FxrIdxSt, nLocalFxr, nFsrInFxr
INTEGER :: i, j
INTEGER :: nofg, norg
INTEGER :: xyb, xye   !--- CNJ Edit : Domain Decomposition + MPI
LOGICAL :: lresogrp(ng), lres

Pin => Core%Pin
CellInfo => Core%CellInfo
nCoreFsr = Core%nCoreFsr
nCoreFxr = Core%nCoreFxr
nxy = Core%nxy
xyb = PE%myPinBeg; xye = PE%myPinEnd   !--- CNJ Edit : Domain Decomposition + MPI

IF(lxslib) THEN
  nofg = GroupInfo%nofg; norg = GroupInfo%norg
  DO ig = 1, ng
    lresogrp(ig) = FALSE
    if(ig .gt. nofg .and. ig .le. (nofg + norg)) lresogrp(ig) =TRUE
  ENDDO
ENDIF

DO ipin = xyb, xye
  FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
  icel = Pin(ipin)%Cell(iz)
  nlocalFxr = CellInfo(icel)%nFxr 
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j -1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    myFxr => Fxr(ifxr, iz)
    IF(lxslib) THEN
      CALL MacXsBase(XSMac, myFxr, 1, ng, ng, 1._8, FALSE, TRUE, TRUE)
      DO ig = 1, ng
        CALL BaseMacSTr(XsMac, myFxr, ig, ng, TRUE)  
        lres = lresogrp(ig) .and. myFxr%lres
        IF(lres) THEN 
          XsMac%XsMacA(ig) = XsMac%XsMacA(ig) * myFxr%fresoa(ig)
          IF (lRST) XsMac%XsMacStr(ig) = XsMac%XsMacStr(ig) * myFxr%fresostr(ig)
        ENDIF
        xsmactr(ig) = XsMac%XsMacA(ig) + XsMac%XsMacstr(ig)
      ENDDO
    ELSE
      ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(1, j)
      itype = myFxr%imix
      IF(TranCntl%lDynamicBen) THEN
        xsmactr = GetXsTrDynBen(itype, TranInfo%fuelTemp(ipin, iz), 1, ng)
      ELSE
        xsmactr = GetXsTrBen(itype, 1, ng)
      END IF
    ENDIF
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      DO ig = 1, ng
        IF((xstnm(ig, ifsr)) .LT. epsm3) CYCLE
        XsD = 1._8 / (3._8 * xsmactr(ig)); XsD = XsD * Bsq
        xstnm(ig, ifsr) = xstnm(ig, ifsr) + XsD
      ENDDO
    ENDDO
  ENDDO
ENDDO

NULLIFY(Pin)
NULLIFY(CellInfo)
NULLIFY(myFxr)

END SUBROUTINE AddBucklingNM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetRtMacXsNM_Cusping(core, FmInfo, Fxr, xstnm, iz, ng, lxslib, lTrCorrection, lRST, lssph, lssphreg, PE)

USE PARAM
USE TYPEDEF,  ONLY : coreinfo_type,    Fxrinfo_type,      Cell_Type, pin_Type,   &
                     PE_TYPE,          FmInfo_type,                              &
                     XsMac_Type
USE Core_mod, ONLY : GroupInfo

USE BenchXs,  ONLY : GetXstrBen,       MacXsBen,          GetXsTrBen_Cusping,    GetXsTrDynBen,         GetXsTrDynBen_Cusping, &
                     GetXstBen,        GetXstBen_Cusping, GetXstDynBen,          GetXstDynBen_Cusping,  DynMacXsBen
USE MacXsLib_mod,   ONLY : MacXsBase,  BaseMacSTr,        MacXsBase_Cusping,     BaseMacSTr_Cusping
USE OMP_LIB
USE SPH_mod,  ONLY : ssphfnm
USE XSLIB_MOD,ONLY : igresb,igrese
USE TRAN_MOD, ONLY : TranInfo,         TranCntl
IMPLICIT NONE
TYPE(coreinfo_type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(Fxrinfo_type) :: Fxr(:)
TYPE(PE_TYPE) :: PE
REAL, POINTER :: xstnm(:, :)
INTEGER :: iz, ng
logical :: lxslib, lTrCorrection, lRST, lssph, lssphreg

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
!TYPE(FxrInfo_Type), POINTER, SAVE :: myFxr
INTEGER :: nCoreFxr, nCoreFsr, nxy, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr
INTEGER :: icel, ipin, ifxr, ifsr, ifsrlocal, itype, ig, tid
INTEGER :: i, j, k
INTEGER :: nofg, norg
INTEGER :: xyb, xye   !--- CNJ Edit : Domain Decomposition + MPI
REAL, ALLOCATABLE :: phiz(:, :), philz(:, :), phiuz(:, :)
REAL :: vol, volsum
REAL :: wt, wtbar, wtg, wtgbar, uflux, lflux, hz, hzl, hzu
REAL :: xsmactr(ng)
REAL :: SPHfac(ng,40)

LOGICAL :: lresogrp(ng), lres, lress
LOGICAL :: lCusping

REAL :: macstr
TYPE(XsMac_Type), POINTER :: XsMac(:)
TYPE(XsMac_Type), SAVE :: XsMac0(nThreadMax), XsMac1(nThreadMax)

Pin => Core%Pin
CellInfo => Core%CellInfo
nCoreFsr = Core%nCoreFsr
nCoreFxr = Core%nCoreFxr
nxy = Core%nxy
xyb = PE%myPinBeg; xye = PE%myPinEnd   !--- CNJ Edit : Domain Decomposition + MPI

ALLOCATE(XsMac(PE%nThread))
ALLOCATE(phiz(ng, PE%nThread), philz(ng, PE%nThread), phiuz(ng, PE%nThread))

IF(lxslib) THEN
  nofg = GroupInfo%nofg; norg = GroupInfo%norg
  DO ig = 1, ng
    lresogrp(ig) = FALSE
    IF(ig .gt. nofg .and. ig .le. (nofg + norg)) lresogrp(ig) =TRUE
  ENDDO
ENDIF
SPHfac=1._8

!$OMP PARALLEL DEFAULT(SHARED)                                                                                      &
!$OMP PRIVATE(i, j, k, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr, vol, volsum,                                       &
!$OMP         icel, ipin, ifxr, ifsr, ifsrlocal, itype, ig, tid, xsmactr, lres, lress, SPHfac)
tid = omp_get_thread_num() + 1
!$OMP DO SCHEDULE(DYNAMIC)
DO ipin = xyb, xye
  FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
  icel = Pin(ipin)%Cell(iz)
  nlocalFxr = CellInfo(icel)%nFxr 
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j -1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    IF(lxslib) THEN
      IF(Fxr(ifxr)%lCusping) THEN
        phiz(:,tid) = 0.
        philz(:,tid) = 0.
        phiuz(:,tid) = 0.
        volsum = 0.
        DO ig = 1, ng
          DO i = 1, nFsrInFxr
            ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(i, j)
            ifsr = FsrIdxSt + ifsrlocal - 1
            vol = CellInfo(icel)%vol(ifsrlocal)
            IF(ig .EQ. 1) volsum = volsum + vol
            phiz(ig, tid) = phiz(ig, tid) + fmInfo%phis(ifsr, iz, ig) * vol
            IF(iz .EQ. PE%myzb) THEN
              philz(ig, tid) = philz(ig, tid) + FmInfo%neighphis(ifsr, ig, BOTTOM) * vol
            ELSE
              philz(ig, tid) = philz(ig, tid) + FmInfo%phis(ifsr, iz-1, ig) * vol
            END IF
            IF(iz .EQ. PE%myze) THEN
              phiuz(ig, tid) = phiuz(ig, tid) + FmInfo%neighphis(ifsr, ig, TOP) * vol
            ELSE
              phiuz(ig, tid) = phiuz(ig, tid) + FmInfo%phis(ifsr, iz+1, ig) * vol
            END IF
          END DO
          phiz(ig, tid) = phiz(ig, tid) / volsum
          philz(ig, tid) = philz(ig, tid) / volsum
          phiuz(ig, tid) = phiuz(ig, tid) / volsum
        END DO
        CALL MacXsBase_Cusping(XSMac0(tid), Fxr(ifxr)%iso0, Fxr(ifxr), 1, ng, ng, 1., FALSE, TRUE)
        CALL MacXsBase_Cusping(XSMac1(tid), Fxr(ifxr)%iso1, Fxr(ifxr), 1, ng, ng, 1., FALSE, TRUE)
        CALL MacXsBase(XSMac(tid), Fxr(ifxr), 1, ng, ng, 1._8, FALSE, TRUE, TRUE)
        IF (lTrCorrection) THEN
          DO ig = 1, ng
            CALL BaseMacStr_Cusping(XsMac0(tid), Fxr(ifxr)%iso0, Fxr(ifxr), ig, ng)
            CALL BaseMacStr_Cusping(XsMac1(tid), Fxr(ifxr)%iso1, Fxr(ifxr), ig, ng)
            CALL BaseMacStr(XsMac(tid), Fxr(ifxr), ig, ng, .TRUE.)
          ENDDO
        ENDIF
        wt = Fxr(ifxr)%wt
        wtbar = 1._8 - wt
        hz = Core%hzfm(iz); hzl = Core%hzfm(iz-1); hzu = Core%hzfm(iz+1)
        DO ig = 1, ng
          uflux = (hzu * phiuz(ig, tid) + wt * hz * phiz(ig, tid)) / (hzu + wt * hz)
          lflux = (hzl * philz(ig, tid) + wtbar * hz * phiz(ig, tid)) / (hzl + wtbar * hz)
          wtg = wt * uflux / (wt * uflux + wtbar * lflux)
          wtgbar = 1. - wtg

          lres = lresogrp(ig) .and. Fxr(ifxr)%lres
          lress = lres .and. lRST

          XsMac(tid)%XsMacA(ig) = XsMac0(tid)%XsMacA(ig) * wtgbar + XsMac1(tid)%XsMacA(ig) * wtg
          IF(lres) XsMac(tid)%XsMacA(ig) = XsMac(tid)%XsMacA(ig) * Fxr(ifxr)%fresoa(ig)
          IF(lTrCorrection) THEN
            XsMac(tid)%XsMacStr(ig) = XsMac0(tid)%XsMacStr(ig) * wtgbar + XsMac1(tid)%XsMacStr(ig) * wtg
            IF(lress) XsMac(tid)%XsMacStr(ig) = XsMac(tid)%XsMacStr(ig) * Fxr(ifxr)%fresostr(ig)
            xsmactr(ig) = XsMac(tid)%XsMacA(ig) + XsMac(tid)%XsMacstr(ig)
          ELSE
            XsMac(tid)%XsMacS(ig) = XsMac0(tid)%XsMacS(ig) * wtgbar + XsMac1(tid)%XsMacS(ig) * wtg
            IF(lress) XsMac(tid)%XsMacS(ig) = XsMac(tid)%XsMacS(ig) * Fxr(ifxr)%fresos(ig)
            xsmactr(ig) = XsMac(tid)%XsMacA(ig) + XsMac(tid)%XsMacs(ig)
          END IF
#ifdef inflow      
          xsmactr(ig) = xsmactr(ig) + fxr(ifxr)%Delinflow(ig)
#endif
        END DO
      ELSE
        CALL MacXsBase(XSMac(tid), Fxr(ifxr), 1, ng, ng, 1._8, FALSE, TRUE, TRUE)
        DO ig = 1, ng
          lres = lresogrp(ig) .and. Fxr(ifxr)%lres
          lress = lres .and. lRST
          IF(lres) XsMac(tid)%XsMacA(ig) = XsMac(tid)%XsMacA(ig) * Fxr(ifxr)%fresoa(ig)
          IF(lTrCorrection) THEN
            CALL BaseMacSTr(XsMac(tid), Fxr(ifxr), ig, ng, TRUE)
            IF(lress) XsMac(tid)%XsMacStr(ig) = XsMac(tid)%XsMacStr(ig) * Fxr(ifxr)%fresostr(ig)
            xsmactr(ig) = XsMac(tid)%XsMacA(ig) + XsMac(tid)%XsMacstr(ig)
          ELSE
            IF(lress) XsMac(tid)%XsMacS(ig) = XsMac(tid)%XsMacS(ig) * Fxr(ifxr)%fresos(ig)
            xsmactr(ig) = XsMac(tid)%XsMacA(ig) + XsMac(tid)%XsMacs(ig)
          ENDIF
#ifdef inflow      
          xsmactr(ig) = xsmactr(ig) + fxr(ifxr)%Delinflow(ig)
#endif
        ENDDO
      END IF

      IF (lsSPH) THEN
        DO ig = igresb, igrese
          IF (CellInfo(icel)%lsSPH) THEN
            IF (lsSPHreg) THEN
              SPHfac(ig,j)=Fxr(ifxr)%SPHfactor(ig)
            ELSE
              SPHfac(ig,j)=CellInfo(icel)%SPHfactor(j,ig)
            ENDIF
            xsmactr(ig)=xsmactr(ig)*SPHfac(ig,j)
          ENDIF
        ENDDO
      ENDIF
    ELSE
      itype = Fxr(ifxr)%imix
      IF(TranCntl%lDynamicBen) THEN
        lCusping = DynMacXsBen(itype)%lCusping
      ELSE
        lCusping = MacXsBen(itype)%lCusping
      END IF
      IF(lCusping) THEN
        phiz(:, tid) = 0.
        philz(:, tid) = 0.
        phiuz(:, tid) = 0.
        volsum = 0.
        DO ig = 1, ng
          DO i = 1, nFsrInFxr
            ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(i,j)
            ifsr = FsrIdxSt + ifsrlocal - 1
            vol = CellInfo(icel)%Vol(ifsrlocal)
            IF(ig .EQ. 1) volsum = volsum + vol
            phiz(ig, tid) = phiz(ig, tid) + FmInfo%phis(ifsr, iz, ig) * vol
            IF(iz .EQ. PE%myzb) THEN
              philz(ig, tid) = philz(ig, tid) + FmInfo%neighphis(ifsr, ig, BOTTOM) * vol
            ELSE
              philz(ig, tid) = philz(ig, tid) + FmInfo%phis(ifsr, iz-1, ig) * vol 
            END IF
            IF(iz .EQ. PE%myze) THEN
              phiuz(ig, tid) = phiuz(ig, tid) + FmInfo%neighphis(ifsr, ig, TOP) * vol
            ELSE
              phiuz(ig, tid) = phiuz(ig, tid) + FmInfo%phis(ifsr, iz+1, ig) * vol
            END IF
          END DO 
          phiz(ig, tid) = phiz(ig, tid) / volsum
          philz(ig, tid) = philz(ig, tid) / volsum
          phiuz(ig, tid) = phiuz(ig, tid) / volsum
        END DO
        IF(lTrCorrection) THEN
          IF(TranCntl%lDynamicBen) THEN
            xsmactr = GetXsTrDynBen_Cusping(itype, TranInfo%fuelTemp(ipin, iz), 1, ng, phiz(:, tid), philz(:, tid), phiuz(:, tid), &
              Core%hzfm(iz), Core%hzfm(iz-1), Core%hzfm(iz+1))
          ELSE
            xsmactr = GetXsTrBen_Cusping(itype, 1, ng, phiz(:, tid), philz(:, tid), phiuz(:, tid), &
              Core%hzfm(iz), Core%hzfm(iz-1), Core%hzfm(iz+1))
          END IF
        ELSE
          IF(TranCntl%lDynamicBen) THEN
            xsmactr = GetXsTDynBen_Cusping(itype, TranInfo%fuelTemp(ipin, iz), 1, ng, phiz(:, tid), philz(:, tid), phiuz(:, tid), &
              Core%hzfm(iz), Core%hzfm(iz-1), Core%hzfm(iz+1))
          ELSE
            xsmactr = GetXsTBen_Cusping(itype, 1, ng, phiz(:, tid), philz(:, tid), phiuz(:, tid), &
              Core%hzfm(iz), Core%hzfm(iz-1), Core%hzfm(iz+1))
          END IF
        END IF
      ELSE
        IF(lTrCorrection) THEN
          IF(TranCntl%lDynamicBen) THEN
            xsmactr = GetXsTrDynBen(itype, TranInfo%fuelTemp(ipin, iz), 1, ng)
          ELSE
            xsmactr = GetXsTrBen(itype, 1, ng)
          END IF
        ELSE
          IF(TranCntl%lDynamicBen) THEN
            xsmactr = GetXsTDynBen(itype, TranInfo%fuelTemp(ipin, iz), 1, ng)
          ELSE
            xsmactr = GetXsTBen(itype, 1, ng)
          END IF
        END IF
      END IF
    ENDIF
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      xstnm(:, ifsr) = xsmactr
      IF (lsSPH) THEN
        IF (CellInfo(icel)%lsSPH) ssphfnm(igresb:igrese,ifsr,iz) = SPHfac(igresb:igrese,j)
      ENDIF
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

DEALLOCATE(XsMac)
DEALLOCATE(phiz, philz, phiuz)

NULLIFY(Pin)
NULLIFY(CellInfo)

END SUBROUTINE SetRtMacXsNM_Cusping
! ------------------------------------------------------------------------------------------------------------