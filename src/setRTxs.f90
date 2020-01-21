#include <defines.h>
SUBROUTINE SetRtMacXs(core, Fxr, xstr, iz, ig, ng, lxslib, lTrCorrection, lRST, lsSPH, lsSPHreg, PE)
USE PARAM
USE TYPEDEF,  ONLY : coreinfo_type,    Fxrinfo_type,      Cell_Type, pin_Type,   &
                     PE_TYPE,                                                    &
                     XsMac_Type
USE Core_mod, ONLY : GroupInfo

USE BenchXs,  ONLY : GetXstrBen, GetXsTBen
USE MacXsLib_mod,   ONLY : MacXsBase,  BaseMacSTr
USE OMP_LIB
USE SPH_mod,  ONLY : ssphf
USE XSLIB_MOD,ONLY : igresb,igrese
IMPLICIT NONE
TYPE(coreinfo_type) :: Core
TYPE(Fxrinfo_type) :: Fxr(:)
TYPE(PE_TYPE) :: PE
REAL, POINTER :: xstr(:)
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
REAL :: xsmactr(ng)

LOGICAL :: lresogrp, lres, lress

REAL :: macstr
TYPE(XsMac_Type), SAVE :: XsMac(nThreadMax) 

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
!$OMP PRIVATE(i, j, k, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr, icel, ipin, ifxr, ifsr, ifsrlocal, itype, tid, xsmactr, lres, lress, SPHfac)
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
      ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(1,j)
      !itype = CellInfo(icel)%iReg(ifsrlocal)
      itype = Fxr(ifxr)%imix
      IF(lTrCorrection) THEN
        xsmactr(ig:ig) = GetXsTrBen(itype, ig, ig)
      ELSE
        xsmactr(ig:ig) = GetXsTBen(itype, ig, ig)
      ENDIF
          
      !CALL xstrben(itype, ig, ig, xsmactr)
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
!NULLIFY(myFxr)
END SUBROUTINE

!--- CNJ Edit : Node Majors    
SUBROUTINE SetRtMacXsNM(core, Fxr, xstnm, iz, ng, lxslib, lTrCorrection, lRST, lssph, lssphreg, PE)
USE PARAM
USE TYPEDEF,  ONLY : coreinfo_type,    Fxrinfo_type,      Cell_Type, pin_Type,   &
                     PE_TYPE,                                                    &
                     XsMac_Type
USE Core_mod, ONLY : GroupInfo

USE BenchXs,  ONLY : GetXstrBen
USE MacXsLib_mod,   ONLY : MacXsBase,  BaseMacSTr
USE OMP_LIB
USE SPH_mod,  ONLY : ssphfnm
USE XSLIB_MOD,ONLY : igresb,igrese
IMPLICIT NONE
TYPE(coreinfo_type) :: Core
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
REAL :: xsmactr(ng)
REAL :: SPHfac(ng,40)

LOGICAL :: lresogrp(ng), lres, lress

REAL :: macstr
TYPE(XsMac_Type), POINTER :: XsMac(:)

Pin => Core%Pin
CellInfo => Core%CellInfo
nCoreFsr = Core%nCoreFsr
nCoreFxr = Core%nCoreFxr
nxy = Core%nxy
xyb = PE%myPinBeg; xye = PE%myPinEnd   !--- CNJ Edit : Domain Decomposition + MPI

ALLOCATE(XsMac(PE%nThread))

IF(lxslib) THEN
  nofg = GroupInfo%nofg; norg = GroupInfo%norg
  DO ig = 1, ng
    lresogrp(ig) = FALSE
    IF(ig .gt. nofg .and. ig .le. (nofg + norg)) lresogrp(ig) =TRUE
  ENDDO
ENDIF
SPHfac=1._8

!$OMP PARALLEL DEFAULT(SHARED)                                                                                      &
!$OMP PRIVATE(i, j, k, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr,                                                    &
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
      ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(1,j)
      itype = Fxr(ifxr)%imix
      xsmactr = GetXsTrBen(itype, 1, ng)
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

NULLIFY(Pin)
NULLIFY(CellInfo)
END SUBROUTINE

SUBROUTINE AddBuckling(Core, Fxr, xstr1g, bsq, iz, ig, ng, lxslib, lRST)
USE PARAM
USE TYPEDEF,      ONLY : coreinfo_type,       Cell_Type,      pin_Type,      &
                         Fxrinfo_type,        XsMac_Type
USE Core_mod,     ONLY : GroupInfo
USE BenchXs,  ONLY : GetXstrBen
USE MacXsLib_mod,   ONLY : MacXsBase,  BaseMacSTr
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
      xsmactr(ig:ig) = GetXsTrBen(itype, ig, ig)
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
END SUBROUTINE
    
!--- CNJ Edit : Node Majors    
SUBROUTINE AddBucklingNM(Core, Fxr, xstnm, bsq, iz, ng, lxslib, lRST)
USE PARAM
USE TYPEDEF,      ONLY : coreinfo_type,       Cell_Type,      pin_Type,      &
                         Fxrinfo_type,        XsMac_Type
USE Core_mod,     ONLY : GroupInfo
USE BenchXs,      ONLY : GetXstrBen
USE MacXsLib_mod, ONLY : MacXsBase,           BaseMacSTr
USE PE_MOD,       ONLY : PE
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
      xsmactr = GetXsTrBen(itype, 1, ng)
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
END SUBROUTINE    