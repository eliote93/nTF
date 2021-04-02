#include <defines.h>
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetRtMacXsGM(core, Fxr, xstr, iz, ig, ng, lxslib, lTrCorrection, lRST, lsSPH, lsSPHreg, PE)

USE PARAM
USE OMP_LIB
USE TYPEDEF,      ONLY : coreinfo_type, Fxrinfo_type, Cell_Type, pin_Type, PE_TYPE, XsMac_Type
USE Core_mod,     ONLY : GroupInfo
USE BenchXs,      ONLY : GetXstrBen, GetXsTBen, GetXsTrDynBen,  GetXsTDynBen
USE MacXsLib_mod, ONLY : MacXsBase, BaseMacSTr
USE SPH_mod,      ONLY : ssphf
USE XSLIB_MOD,    ONLY : igresb,igrese
USE TRAN_MOD,     ONLY : TranInfo, TranCntl

IMPLICIT NONE

TYPE(coreinfo_type) :: Core
TYPE(Fxrinfo_type), DIMENSION(:) :: Fxr
TYPE(PE_TYPE) :: PE
REAL, POINTER, DIMENSION(:) :: xstr

INTEGER :: iz, ig, ng
LOGICAL :: lxslib, lTrCorrection, lRST, lsSPH, lsSPHreg
! ----------------------------------------------------
INTEGER :: nCoreFxr, nCoreFsr, nxy, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr, nofg, norg
INTEGER :: icel, ipin, ifxr, ifsr, ifsrlocal, itype, tid
INTEGER :: i, j, k

LOGICAL :: lresogrp, lres, lress

REAL :: macstr

REAL, DIMENSION(40) :: SPHfac
REAL, DIMENSION(ng) :: xsmactr

TYPE(XsMac_Type), SAVE, DIMENSION(nThreadMax) :: XsMac

TYPE(Pin_Type),  POINTER, DIMENSION(:) :: Pin
TYPE(Cell_Type), POINTER, DIMENSION(:) :: CellInfo
! ----------------------------------------------------

Pin      => Core%Pin
CellInfo => Core%CellInfo
nCoreFsr  = Core%nCoreFsr
nCoreFxr  = Core%nCoreFxr
nxy       = Core%nxy

tid = 1

IF (lxslib) THEN
  nofg = GroupInfo%nofg
  norg = GroupInfo%norg
  
  lresogrp = FALSE
  IF (ig.GT.nofg .AND. ig.LE.(nofg+norg)) lresogrp = TRUE
END IF

SPHfac = ONE
! ----------------------------------------------------
!$  call omp_set_dynamic(FALSE)
!$  call omp_set_num_threads(PE%nThread) 

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr, icel, ipin, ifxr, ifsr, ifsrlocal, itype, tid, xsmactr, lres, lress, SPHfac)
!$  tid = omp_get_thread_num()+1
!$OMP DO
DO ipin = 1, nxy
  FsrIdxSt  = Pin(ipin)%FsrIdxSt
  FxrIdxSt  = Pin(ipin)%FxrIdxSt
  icel      = Pin(ipin)%Cell(iz)
  nlocalFxr = CellInfo(icel)%nFxr 
  
  DO j = 1, nLocalFxr
    ifxr      = FxrIdxSt + j -1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    ! ----------------------------------------------------
    IF (lxslib) THEN
      !Obtain Absoprtion Xs
      CALL MacXsBase(XSMac(tid), Fxr(ifxr), ig, ig, ng, ONE, FALSE, TRUE, TRUE)  
      
      lres  = lresogrp .AND. Fxr(ifxr)%lres
      lress = lres .AND. lRST
      
      ! Self-Sheilding Effect of Absorption/Scattering XS      
      IF (lres) XsMac(tid)%XsMacA(ig) = XsMac(tid)%XsMacA(ig) * Fxr(ifxr)%fresoA(ig)
      
      ! Obtain Outgoing total XS(trasport Corrected)
      IF (lTrCorrection) THEN
        CALL BaseMacSTr(XsMac(tid), Fxr(ifxr), ig, ng, TRUE)  
        
        IF (lress) XsMac(tid)%XsMacStr(ig) = XsMac(tid)%XsMacStr(ig) * Fxr(ifxr)%fresoStr(ig)
        
        xsmactr(ig) = XsMac(tid)%XsMacA(ig) + XsMac(tid)%XsMacstr(ig)
      ELSE
        IF (lress) XsMac(tid)%XsMacS(ig) = XsMac(tid)%XsMacS(ig) * Fxr(ifxr)%fresoS(ig)
        
        xsmactr(ig) = XsMac(tid)%XsMacA(ig) + XsMac(tid)%XsMacs(ig)
      ENDIF
#ifdef inflow
      xsmactr(ig) = xsmactr(ig) + fxr(ifxr)%Delinflow(ig)
#endif
      IF (lsSPH .AND. ig.GE.igresb .AND. ig.LE.igrese .AND. CellInfo(icel)%lsSPH) THEN
        IF (lsSPHreg) THEN
          SPHfac(j) = Fxr(ifxr)%SPHfactor(ig)
        ELSE
          SPHfac(j) = CellInfo(icel)%SPHfactor(j, ig)
        END IF
        
        xsmactr(ig) = xsmactr(ig) * SPHfac(j)
      ENDIF
    ! ----------------------------------------------------
    ELSE
      ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(1,j)
      itype     = Fxr(ifxr)%imix
      
      IF (lTrCorrection) THEN
        IF (TranCntl%lDynamicBen) THEN
          xsmactr(ig:ig) = GetXsTrDynBen(itype, TranInfo%fuelTemp(ipin, iz), ig, ig)
        ELSE
          xsmactr(ig:ig) = GetXsTrBen(itype, ig, ig)
        END IF
      ELSE
        IF(TranCntl%lDynamicBen) THEN
          xsmactr(ig:ig) = GetXsTDynBen(itype, TranInfo%fuelTemp(ipin, iz), ig, ig)
        ELSE
          xsmactr(ig:ig) = GetXsTBen(itype, ig, ig)
        END IF
      END IF
    END IF
    ! ------------------------------------------------
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      
      xstr(ifsr) = xsmactr(ig)
      
      IF (lsSPH .AND. ig.GE.igresb .AND. ig.LE.igrese .AND. CellInfo(icel)%lsSPH) ssphf(ifsr,iz,ig) = SPHfac(j)
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL
! ----------------------------------------------------
NULLIFY (Pin)
NULLIFY (CellInfo)
! ----------------------------------------------------

END SUBROUTINE SetRtMacXsGM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE AddBucklingGM(Core, Fxr, xstr1g, bsq, iz, ig, ng, lxslib, lRST)

USE PARAM
USE TYPEDEF,      ONLY : coreinfo_type,       Cell_Type,      pin_Type,      &
                         Fxrinfo_type,        XsMac_Type
USE Core_mod,     ONLY : GroupInfo
USE BenchXs,      ONLY : GetXstrBen,          GetXstrDynBen  
USE MacXsLib_mod,   ONLY : MacXsBase,  BaseMacSTr
USE TRAN_MOD,     ONLY : TranInfo,            TranCntl
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(Fxrinfo_type), POINTER :: Fxr(:,:)
REAL, POINTER :: xstr1g(:)
INTEGER :: iz, ig, ng
REAL :: Bsq
LOGICAL :: lxslib, lRST
TYPE(FxrInfo_Type), POINTER, SAVE :: myFxr
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(XsMac_Type), SAVE :: XsMac 
REAL :: XsD
REAL :: xsmactr(ng)
INTEGER :: nxy, nCoreFsr, nCoreFxr
INTEGER :: FsrIdxSt, ipin, icel, ifsr, ifxr, itype, ifsrlocal
INTEGER :: FxrIdxSt, nLocalFxr, nFsrInFxr
INTEGER :: i, j
INTEGER :: nofg, norg
LOGICAL :: lresogrp, lres

Pin => Core%Pin
CellInfo => Core%CellInfo
nCoreFsr = Core%nCoreFsr
nCoreFxr = Core%nCoreFxr
nxy = Core%nxy

IF(lxslib) THEN
  nofg = GroupInfo%nofg; norg = GroupInfo%norg
  lresogrp = FALSE
  if(ig.gt.nofg.and.ig.le.(nofg+norg)) lresogrp =TRUE
ENDIF

DO ipin = 1, nxy
  FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
  icel = Pin(ipin)%Cell(iz)
  nlocalFxr = CellInfo(icel)%nFxr 
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j -1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    !Get Transport XS
    myFxr => Fxr(ifxr, iz)
    IF(lxslib) THEN
      
      !Obtain Absoprtion Xs
      CALL MacXsBase(XSMac, myFxr, ig, ig, ng, 1._8, FALSE, TRUE, TRUE)  
      CALL BaseMacSTr(XsMac, myFxr, ig, ng, TRUE)  
      !Self-Sheilding Effect of Absorption XS
      lres = lresogrp .and. myFxr%lres
      IF(lres) THEN 
        XsMac%XsMacA(ig) = XsMac%XsMacA(ig) * myFxr%fresoa(ig)
        IF (lRST) XsMac%XsMacStr(ig) = XsMac%XsMacStr(ig) * myFxr%fresostr(ig)
      ENDIF
      !Obtain Outgoing total XS(trasport Corrected)      
      xsmactr(ig) = XsMac%XsMacA(ig) + XsMac%XsMacstr(ig)
    ELSE
      ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(1,j)
      !itype = CellInfo(icel)%iReg(ifsrlocal)
      itype = myFxr%imix
      IF(TranCntl%lDynamicBen) THEN
        xsmactr(ig:ig) = GetXsTrDynBen(itype, TranInfo%fuelTemp(ipin, iz), ig, ig)
      ELSE
        xsmactr(ig:ig) = GetXsTrBen(itype, ig, ig)
      END IF
      !CALL xstrben(itype, ig, ig, xsmactr)
    ENDIF
    !Assgin  
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      IF((xstr1g(ifsr)) .LT. epsm3) CYCLE
      XsD = 1._8 / (3._8 * xsmactr(ig)); XsD = XsD*Bsq
      xstr1g(ifsr) = xstr1g(ifsr) + XsD
      !xstr(ifsr, iz) = xstr(ifsr, iz)
    ENDDO
  ENDDO !End of Fxr Sweep
ENDDO !End of Pin Sweep

NULLIFY(Pin)
NULLIFY(CellInfo)
NULLIFY(myFxr)

END SUBROUTINE AddBucklingGM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetRtMacXsGM_Cusping(core, FmInfo, Fxr, xstr, phis, iz, ig, ng, lxslib, lTrCorrection, lRST, lsSPH, lsSPHreg, PE)

USE PARAM
USE TYPEDEF,  ONLY : coreinfo_type,    Fxrinfo_type,      Cell_Type, pin_Type,   &
                     PE_TYPE,          FmInfo_Type,                              &
                     XsMac_Type
USE Core_mod, ONLY : GroupInfo

USE BenchXs,  ONLY : GetXstrBen, GetXsTBen, GetXstrBen_Cusping, GetXstBen_Cusping, MacXsBen, &
                     GetXstrDynBen, GetXsTDynBen, GetXstrDynBen_Cusping, GetXstDynBen_Cusping, DynMacXsBen
USE MacXsLib_mod,   ONLY : MacXsBase,  BaseMacSTr
USE OMP_LIB
USE SPH_mod,  ONLY : ssphf
USE XSLIB_MOD,ONLY : igresb,igrese
USE TRAN_MOD, ONLY : TranInfo,  TranCntl
IMPLICIT NONE
TYPE(coreinfo_type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(Fxrinfo_type) :: Fxr(:)
TYPE(PE_TYPE) :: PE
REAL, POINTER :: xstr(:)
REAL, POINTER :: phis(:, :, :)
REAL :: SPHfac(40)
INTEGER :: iz, ig, ng
logical :: lxslib, lTrCorrection, lRST, lsSPH, lsSPHreg


TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
!TYPE(FxrInfo_Type), POINTER, SAVE :: myFxr
INTEGER :: nCoreFxr, nCoreFsr, nxy, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr
INTEGER :: icel, ipin, ifxr, ifsr, ifsrlocal, itype, tid
INTEGER :: i, j, k
INTEGER :: nofg, norg
REAL :: phiz(1), philz(1), phiuz(1)
REAL :: vol, volsum
REAL :: xsmactr(ng)

LOGICAL :: lresogrp, lres, lress
LOGICAL :: lCusping

REAL :: macstr
TYPE(XsMac_Type), SAVE :: XsMac(nThreadMax) 
TYPE(XsMac_Type), SAVE :: XsMac0(nThreadMax), XsMac1(nThreadMax)

Pin => Core%Pin
CellInfo => Core%CellInfo
nCoreFsr = Core%nCoreFsr
nCoreFxr = Core%nCoreFxr
nxy = Core%nxy
tid = 1
IF(lxslib) THEN
  nofg = GroupInfo%nofg; norg = GroupInfo%norg
  lresogrp = FALSE
  if(ig.gt.nofg.and.ig.le.(nofg+norg)) lresogrp =TRUE
ENDIF
SPHfac=1._8
!$  call omp_set_dynamic(.FALSE.)
!$  call omp_set_num_threads(PE%nThread) 
!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(i, j, k, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr, icel, ipin, ifxr, ifsr, ifsrlocal, itype, tid, xsmactr, lres, lress, SPHfac, phiz, philz, phiuz, vol, volsum)
!$  tid = omp_get_thread_num()+1
!$OMP DO
DO ipin = 1, nxy
  FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
  icel = Pin(ipin)%Cell(iz)
  nlocalFxr = CellInfo(icel)%nFxr 
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j -1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    !Get Transport XS
      !myFxr => Fxr(ifxr)
    IF(lxslib) THEN
      !Obtain Absoprtion Xs
      CALL MacXsBase(XSMac(tid), Fxr(ifxr), ig, ig, ng, 1._8, FALSE, TRUE, TRUE)  
      
      lres = lresogrp .and. Fxr(ifxr)%lres
      lress = lres .and. lRST
      
      !Self-Sheilding Effect of Absorption/Scattering XS      
      IF(lres) XsMac(tid)%XsMacA(ig) = XsMac(tid)%XsMacA(ig) * Fxr(ifxr)%fresoA(ig)
      !Obtain Outgoing total XS(trasport Corrected)
      IF(lTrCorrection) THEN
        CALL BaseMacSTr(XsMac(tid), Fxr(ifxr), ig, ng, TRUE)  
        IF(lress) XsMac(tid)%XsMacStr(ig) = XsMac(tid)%XsMacStr(ig) * Fxr(ifxr)%fresoStr(ig)
        xsmactr(ig) = XsMac(tid)%XsMacA(ig) + XsMac(tid)%XsMacstr(ig)
      ELSE
        IF(lress) XsMac(tid)%XsMacS(ig) = XsMac(tid)%XsMacS(ig) * Fxr(ifxr)%fresoS(ig)
        xsmactr(ig) = XsMac(tid)%XsMacA(ig) + XsMac(tid)%XsMacs(ig)
      ENDIF
#ifdef inflow      
      xsmactr(ig) = xsmactr(ig) + fxr(ifxr)%Delinflow(ig)
#endif      
!        CALL MacP1XsScatMatrix(XsMac, myFxr, 1, ng, ng, GroupInfo)
!        CONTINUE
      IF (lsSPH) THEN
        IF (ig.ge.igresb.and.ig.le.igrese) THEN
          IF (CellInfo(icel)%lsSPH) THEN
            IF (lsSPHreg) THEN
              SPHfac(j)=Fxr(ifxr)%SPHfactor(ig)
            ELSE
              SPHfac(j)=CellInfo(icel)%SPHfactor(j,ig)
            ENDIF
            xsmactr(ig)=xsmactr(ig)*SPHfac(j)
          ENDIF
        ENDIF
      ENDIF
    ELSE
      itype = Fxr(ifxr)%imix
      IF(TranCntl%lDynamicBen) THEN
        lCusping = DynMacXsBen(itype)%lCusping
      ELSE
        lCusping = MacXsBen(itype)%lCusping
      END IF
      IF(lCusping) THEN
        phiz = 0.
        philz = 0.
        phiuz = 0.
        volsum = 0. 
        DO i = 1, nFsrInFxr
          ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(i, j)
          ifsr = FsrIdxSt + ifsrlocal - 1
          vol = CellInfo(icel)%Vol(ifsrlocal)
          volsum = volsum + vol
          phiz = phiz + phis(ifsr, iz, ig) * vol
          IF(iz .EQ. PE%myzb) THEN
            philz = philz + FmInfo%neighphis(ifsr, ig, BOTTOM) * vol  
          ELSE
            philz = philz + phis(ifsr, iz-1, ig) * vol
          END IF
          IF(iz .EQ. PE%myze) THEN
            phiuz = phiuz + FmInfo%neighphis(ifsr, ig, TOP) * vol
          ELSE
            phiuz = phiuz + phis(ifsr, iz+1, ig) * vol
          END IF
        END DO 
        phiz = phiz / volsum
        philz = philz / volsum
        phiuz = phiuz / volsum
        IF(lTrCorrection) THEN
          IF(TranCntl%lDynamicBen) THEN
            xsmactr(ig:ig) = GetXsTrDynBen_Cusping(itype, TranInfo%fuelTemp(ipin, iz), ig, ig, phiz, philz, phiuz, &
              Core%hzfm(iz), Core%hzfm(iz-1), Core%hzfm(iz+1))
          ELSE
            xsmactr(ig:ig) = GetXsTrBen_Cusping(itype, ig, ig, phiz, philz, phiuz, &
              Core%hzfm(iz), Core%hzfm(iz-1), Core%hzfm(iz+1))
          END IF
        ELSE
          IF(TranCntl%lDynamicBen) THEN
            xsmactr(ig:ig) = GetXsTDynBen_Cusping(itype, TranInfo%fuelTemp(ipin, iz), ig, ig, phiz, philz, phiuz, &
              Core%hzfm(iz), Core%hzfm(iz-1), Core%hzfm(iz+1))
          ELSE
            xsmactr(ig:ig) = GetXsTBen_Cusping(itype, ig, ig, phiz, philz, phiuz, &
              Core%hzfm(iz), Core%hzfm(iz-1), Core%hzfm(iz+1))
          END IF
        ENDIF
      ELSE
        IF(lTrCorrection) THEN
          IF(TranCntl%lDynamicBen) THEN
            xsmactr(ig:ig) = GetXsTrDynBen(itype, TranInfo%fuelTemp(ipin, iz), ig, ig)
          ELSE
            xsmactr(ig:ig) = GetXsTrBen(itype, ig, ig)
          END IF
        ELSE
          IF(TranCntl%lDynamicBen) THEN
            xsmactr(ig:ig) = GetXsTDynBen(itype, TranInfo%fuelTemp(ipin, iz), ig, ig)
          ELSE
            xsmactr(ig:ig) = GetXsTBen(itype, ig, ig)
          END IF
        ENDIF
      END IF
    ENDIF
    !Assgin 
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      xstr(ifsr) = xsmactr(ig)
      IF (lsSPH) THEN
        IF (ig.ge.igresb.and.ig.le.igrese) THEN
          IF (CellInfo(icel)%lsSPH) ssphf(ifsr,iz,ig) = SPHfac(j)
        ENDIF
      ENDIF
    ENDDO
  ENDDO !End of Fxr Sweep
ENDDO !End of Pin Sweep
!$OMP END DO
!$OMP END PARALLEL
NULLIFY(Pin)
NULLIFY(CellInfo)

END SUBROUTINE SetRtMacXsGM_Cusping
! ------------------------------------------------------------------------------------------------------------