#include <defines.h>
MODULE Inflow_mod
IMPLICIT NONE
INTERFACE

SUBROUTINE InflowTrXS(Core, FmInfo, CmInfo, eigv, GroupInfo, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,             ONLY : CoreInfo_Type,         FmInfo_Type,       PE_Type,         &
                                CmInfo_TYpe,           GroupInfo_Type
USE CNTL,                ONLY : nTracerCntl_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(GroupInfo_Type) ::  GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
REAL :: eigv

END SUBROUTINE 

SUBROUTINE InflowDifCoeff(D, xst, p1sm, phi, ng)
IMPLICIT NONE
INTEGER :: ng
REAL :: D(ng)
REAL :: xst(ng)
REAL :: p1sm(ng, ng)
REAL :: phi(ng)

END SUBROUTINE

SUBROUTINE InflowDifCoeff2(D, xst, p1sm, bsq, phi, ng)
IMPLICIT NONE
INTEGER :: ng
REAL :: D(ng)
REAL :: xst(ng)
REAL :: p1sm(ng, ng)
REAL :: bsq(ng)
REAL :: phi(ng)

END SUBROUTINE

END INTERFACE

END MODULE

SUBROUTINE InflowTrXS(Core, FmInfo, CmInfo, eigv, GroupInfo, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,             ONLY : CoreInfo_Type,         FmInfo_Type,       PE_Type,         &
                                CmInfo_TYpe,           GroupInfo_Type,                     &
                                FxrInfo_Type,          Pin_Type,          Cell_Type,       &
                                XsMac_Type
USE MacXsLib_mod,       ONLY : MacXsBase,              BaseMacStr,      MacP1XsScatMatrix, &
                               MacXsScatMatrix    
USE MOC_Mod,            ONLY : FxrAvgPhi
USE CNTL,               ONLY : nTracerCntl_Type
USE Inflow_mod,         ONLY : InflowDifCoeff
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(GroupInfo_Type) ::  GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
REAL :: eigv

REAL, POINTER :: Phis(:, :, :)
TYPE(Fxrinfo_type), POINTER :: Fxr(:, :)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)

REAL,PARAMETER :: voidcrit = 1.E-4_8

INTEGER:: myzb, myze
INTEGER :: nFsr, nFxr,  nxy, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr, niso
INTEGER :: icel, ig, iz, ixy, ifxr, ifsr, ifsrlocal, itype, tid
INTEGER :: i, j, k
INTEGER :: nofg, norg, ng
LOGICAL :: lresogrp(1000), lres

REAL :: pnsum

REAL, POINTER :: P1sm(:,:)
REAL :: bsq(ngmax)
REAL :: xsmactr(ngmax), xsmact(ngmax), D(ngmax), PhiAvg(ngmax), PsiAvg, areasum
TYPE(XsMac_Type), SAVE :: XsMac(nThreadMax) 


Fxr => FmInfo%Fxr; Phis => FmInfo%Phis

Pin => Core%Pin
CellInfo => Core%CellInfo
nFsr = Core%nCoreFsr; nFxr = Core%nCoreFxr; nxy = Core%nxy
myzb = PE%myzb; myze = PE%myze
nofg = GroupInfo%nofg; norg = GroupInfo%norg; ng = GroupInfo%ng
tid = 1
DO ig = 1, ng
  lresogrp(ig) = FALSE
  if(ig .gt. nofg .and. ig .le. (nofg + norg)) lresogrp(ig) =TRUE
ENDDO
DO iz = myzb, myze
  DO ixy = 1, nxy
    FsrIdxSt = Pin(ixy)%FsrIdxSt; FxrIdxSt = Pin(ixy)%FxrIdxSt
    icel = Pin(ixy)%Cell(iz)
    nlocalFxr = CellInfo(icel)%nFxr
    CALL CalBsq(Bsq(1:ng), CmInfo%PinXS(ixy, iz), CmInfo%PhiC(ixy, iz, 1:ng), CmInfo%PsiC(ixy, iz), eigv, Core%PinVol(ixy, iz), ng)
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j -1
      
      niso = Fxr(ifxr, iz)%niso
      
      pnsum = SUM(Fxr(ifxr, iz)%pnum(1:niso))
      Fxr(ifxr, iz)%lvoid = .FALSE.
      IF(pnsum .LT. VoidCrit) Fxr(ifxr, iz)%lvoid = .TRUE.
      IF(Fxr(ifxr, iz)%lvoid) CYCLE
      nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)      
      !Calculating Total and Transport Cross Section
      
      CALL MacXsBase(XSMac(tid), Fxr(ifxr, iz), 1, ng, ng, 1._8, FALSE, TRUE, TRUE)  
      !Obtain Outgoing total XS(trasport Corrected)
      
      IF(Fxr(ifxr, iz)%lres) THEN
        DO ig = nofg+1, norg
          XsMac(tid)%XsMacA(ig) = XsMac(tid)%XsMacA(ig) * Fxr(ifxr, iz)%fresoa(ig)
          XsMac(tid)%XsMacS(ig) = XsMac(tid)%XsMacS(ig) * Fxr(ifxr, iz)%fresos(ig)
        ENDDO
      ENDIF
      DO ig = 1, ng
        lres = Fxr(ifxr, iz)%lres .and. lresogrp(ig)
        CALL BaseMacStr(XsMac(tid), Fxr(ifxr, iz), ig, ng, TRUE)
        if (lres) XsMac(tid)%XsMacstr(ig) = XsMac(tid)%XsMacstr(ig) * Fxr(ifxr, iz)%fresostr(ig)
        xsmactr(ig) = XsMac(tid)%XsMacA(ig) + XsMac(tid)%XsMacstr(ig)
        xsmact(ig) = XsMac(tid)%XsMacA(ig) + XsMac(tid)%XsMacs(ig)
        D(ig) = 1._8 / (3._8 * xsmactr(ig))
      ENDDO
      !P1 Scattering Matrix
      CALL MacP1XsScatMatrix(XsMac(tid), Fxr(ifxr, iz), 1, ng, ng, GroupInfo)
      P1sm => XsMac(tid)%XsMacP1Sm
      
      PhiAvg(1:ng) = 0
      PsiAvg = 0
      
      DO ig = 1, ng
        DO k= 1, nFsrInFxr
          ifsr = FsrIdxSt + CellInfo(icel)%MapFxr2FsrIdx(k, j) - 1
          areasum = areasum + CellInfo(icel)%vol(CellInfo(icel)%MapFxr2FsrIdx(k, j))
          PhiAvg(ig) = PhiAvg(ig) + CellInfo(icel)%vol(CellInfo(icel)%MapFxr2FsrIdx(k, j)) * phis(ifxr, iz, ig)
          PhiAvg(ig) = PhiAvg(ig) + CellInfo(icel)%vol(CellInfo(icel)%MapFxr2FsrIdx(k, j)) * FmInfo%psi(ifxr, iz)
        ENDDO
        PhiAvg(ig) = PhiAvg(ig) / areasum
      ENDDO
      !Spectrum
      PhiAvg(1:ng) = CmInfo%PhiC(ixy, iz, 1:ng)
      PhiAvg(1:ng) = FxrAvgPhi(Core, Fxr, Phis, ixy, j, iz, ng, PE)
      D(1:ng) = CmInfo%PinXS(ixy, iz)%XSD(1:ng)
      CALL InflowDifCoeff(D(1:ng), xsmact(1:ng), p1sm(1:ng, 1:ng), phiavg(1:ng), ng)
      !CALL InflowDifCoeff2(D(1:ng), xsmact(1:ng), p1sm(1:ng, 1:ng), Bsq(1:ng), phiavg(1:ng), ng)
      DO ig = 1, ng
        D(ig) = 1._8 / (3._8 * D(ig))
        Fxr(ifxr, iz)%delinflow(ig) = D(ig) - xsmactr(ig)
        !Fxr(ifxr, iz)%delinflow(ig) = 0 
      ENDDO
      CONTINUE
    ENDDO
    
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE InflowDifCoeff(D, xst, p1sm, phi, ng)
IMPLICIT NONE
INTEGER :: ng
REAL :: D(ng)
REAL :: xst(ng)
REAL :: p1sm(ng, ng)
REAL :: phi(ng)

REAL :: deno, numer, scat

INTEGER :: ig, ig2
INTEGER :: i
DO i = 1, 10
  DO ig = 1, ng
    deno = 3._8 *(xst(ig) - p1sm(ig, ig)) * phi(ig)
    numer = phi(ig)
    scat = 0
    DO ig2 = 1, ng
      IF(ig2 .EQ. ig) CYCLE
      scat = scat + D(ig2) * phi(ig2) * P1Sm(ig2, ig)
    ENDDO
    numer = numer + 3._8 * Scat
    D(ig) = numer / deno
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE InflowDifCoeff2(D, xst, p1sm, bsq, phi, ng)
IMPLICIT NONE
INTEGER :: ng
REAL :: D(ng)
REAL :: xst(ng), xstr(ng)
REAL :: p1sm(ng, ng)
REAL :: bsq(ng)
REAL :: phi(ng)

REAL :: D0(ng)
REAL :: deno, numer, scat

INTEGER :: ig, ig2
INTEGER :: i
D0 = D
DO i = 1, 1
  DO ig = 1, ng
    !deno = 3._8 *(xst(ig) - p1sm(ig, ig)) * bsq(ig)* phi(ig)
    !numer = phi(ig)* bsq(ig)
    !scat = 0
    !DO ig2 = 1, ng
    !  IF(ig2 .EQ. ig) CYCLE
    !  scat = scat + D(ig2) * phi(ig2) * P1Sm(ig2, ig) * bsq(ig2)
    !ENDDO
    !numer = numer + 3._8 * Scat
    !D(ig) = numer / deno
    deno = xst(ig) 
    DO ig2 = 1, ng
      !IF(ig2 .EQ. ig) CYCLE
      scat = scat + D0(ig2) * phi(ig2) * P1Sm(ig2, ig) * bsq(ig2)
    ENDDO
    scat = scat / D0(ig) * phi(ig) * bsq(ig)
    D(ig) = deno - scat
    D(ig) = 1._8 / (3._8 * D(ig))
  ENDDO
ENDDO
CONTINUE
END SUBROUTINE

SUBROUTINE CalBsq(Bsq, PinXS, PhiC, PsiC, eigv, Vol, ng)
USE PARAM
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE
REAL :: Bsq(ng)
TYPE(PinXS_TYPE) :: PinXS
REAL :: PhiC(ng), PsiC
REAL :: eigv, vol
INTEGER :: ng

REAL :: Src, LKG
INTEGER :: ig, ig2

DO ig = 1, ng
 ! LKG = PinXS%xsr(ig) * PhiC(ig)
  Src = 1._8/eigv * PsiC * PinXS%Chi(ig) / Vol
  DO ig2 = PinXS%XSS(ig)%ib, PinXS%XSS(ig)%ie
    Src = Src + PhiC(ig2) * PinXS%XSS(ig)%From(ig2)
  ENDDO
  Lkg = Src - PinXS%xsr(ig) * PhiC(ig)
  Bsq(ig) = Lkg / (PinXS%XSD(ig) * PhiC(ig))
   Bsq(ig)  = abs(Bsq(ig))
  Lkg = Lkg - Bsq(ig) *(PinXS%XSD(ig) * PhiC(ig)) 
ENDDO


END SUBROUTINE
