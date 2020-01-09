#include <defines.h>
SUBROUTINE UpdtBoronPPM(target_eigv, eigv, ppm, lreset, Master)
USE PARAM
USE Boron_mod,           ONLY : ppmd, eigvd, iter
USE FILES,               ONLY : io8
USE IOUTIL,              ONLY : message
IMPLICIT NONE
REAL :: target_eigv, eigv, ppm
REAL :: delppm
LOGICAL :: lreset, master

REAL :: reigvd, reigv, worthppm
IF(lreset) THEN
  eigvd = 1; iter = 0
ENDIF
iter = iter + 1
reigv = 1._8 / eigv; reigvd = 1._8 / eigvd
worthppm = 5.0_8
IF(iter .GT. 1 .AND. ppm .GT. 0) THEN
  worthppm = 1.0_8 * (reigvd - reigv) * 1.e+5 / (ppmd - ppm)
  worthppm = max(worthppm, 10.0_8)
ENDIF
ppmd =ppm
delppm = (1._8 / target_eigv - reigv) * 1.e+5_8 / worthppm
delppm = min(200._8, delppm)
delppm = max(-200._8, delppm)
ppm = ppm + delppm

IF(ppm .LT. 0) ppm = 0
IF(MASTER) THEN
  write(mesg,'(a,f12.2,a)') ' =======  Applied Boron PPM : ',ppm, ' PPM'
  CALL MESSAGE(io8, TRUE, TRUE, mesg)
ENDIF
eigvd = eigv
END SUBROUTINE

SUBROUTINE SetBoronCoolant(Core, Fxr, boronppm, myzb, myze)
USE PARAM
USE TypeDef,       ONLY : CoreInfo_type, Fxrinfo_type
USE NuclidMap_mod, ONLY : AtomicWeight
USE Boron_mod,     ONLY : b10frac
USE PE_MOD,        ONLY : PE
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
REAL :: boronppm
INTEGER :: myzb, myze

TYPE(FxrInfo_Type), POINTER :: myFxr

INTEGER :: nFxr
INTEGER :: fxrb, fxre   !--- CNJ Edit : Domain Decomposition + MPI
INTEGER :: iz, ifxr, id
INTEGER :: i, j, k
nFxr = Core%nCoreFxr

!--- CNJ Edit : Domain Decomposition + MPI
fxrb = PE%myFxrBeg; fxre = PE%myFxrEnd
IF (PE%RTMASTER) THEN
  fxrb = 1; fxre = nFxr
ENDIF

DO iz = myzb, myze
  DO ifxr = fxrb, fxre
    myFxr => Fxr(ifxr, iz)
    IF(.NOT. myFxr%lH2o) CYCLE
    CALL  FxrBoronUpdate(myFxr, BoronPpm)
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE FxrBoronUpdate(myFxr, BoronPpm)
USE PARAM
USE TypeDef,       ONLY : Fxrinfo_type
USE NuclidMap_mod, ONLY : AtomicWeight
USE Boron_mod,     ONLY : b10frac
USE CrCsp_Mod,     ONLY : CspFxrBoronUpdate
IMPLICIT NONE

TYPE(FxrInfo_Type) :: myFxr
REAL :: Boronppm

INTEGER :: id, i
REAL :: rho, aw
REAL :: tp1, tp2
LOGICAL :: lBoronIstope

rho = 0
DO i = 1, myFxr%niso
  id = myFxr%idiso(i)/1000
  IF(id .EQ. 1 .OR. id .EQ. 8) THEN
    aw = AtomicWeight(myFxr%idiso(i))
    rho = rho + myFxr%pnum(i) * aw / AVOGADRO
  ENDIF
ENDDO
lBoronIstope = FALSE
DO i = 1, myFxr%niso
  !id = myFxr%idiso(i)/1000
  id = myFxr%idiso(i)
  IF(id .EQ. 5000) THEN
    myFxr%pnum(i) = rho * boronppm * epsm6 * avogadro / awboron
    lBoronIstope = TRUE
#define natb
#ifndef natb        
  ELSEIF(id .EQ. 5010) THEN
    myFxr%pnum(i) = rho * boronppm * epsm6 * avogadro / awboron * b10frac
    lBoronIstope = TRUE      
  ELSEIF(id .EQ. 5011) THEN
    myFxr%pnum(i) = rho * boronppm * epsm6 * avogadro / awboron * (1._8-b10frac)
    lBoronIstope = TRUE      
#endif
  ENDIF
ENDDO
IF(.NOT. lBoronIstope) THEN
#ifdef natb
  myFxr%niso = myFxr%niso + 1; i = myFxr%niso
  myFxr%idiso(i) = 5000
  myFxr%pnum(i) = rho * boronppm * epsm6 * avogadro / awboron
#else
  myFxr%niso = myFxr%niso + 2; i = myFxr%niso-1
  myFxr%idiso(i) = 5010;  aw = AtomicWeight(myFxr%idiso(i))
  myFxr%pnum(i) = rho * boronppm * epsm6 * avogadro / awboron * b10frac
      
  i = myFxr%niso
  myFxr%idiso(i) = 5011;  aw = AtomicWeight(myFxr%idiso(i))
  myFxr%pnum(i) = rho * boronppm * epsm6 * avogadro / awboron * (1._8 - b10frac)
  myFxr%niso_depl=myFxr%niso
#endif
ENDIF

IF(myFxr%lCrCspFtn) THEN
  CALL CspFxrBoronUpdate(myFxr, BoronPPM)
ENDIF

END SUBROUTINE

SUBROUTINE SetBoronCoolantCTF(Core, Fxr,ThInfo, myzb, myze)
USE PARAM
USE TypeDef,       ONLY : CoreInfo_type, Fxrinfo_type, Pin_Type,&
                          THInfo_TYPE, Cell_Type
USE NuclidMap_mod, ONLY : AtomicWeight
USE Boron_mod,     ONLY : b10frac
IMPLICIT NONE
TYPE(THInfo_TYPE) :: THInfo
TYPE(CoreInfo_Type) :: Core
TYPE(Pin_Type), POINTER :: Pin (:)
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(Cell_Type), POINTER :: Cell(:)
INTEGER :: myzb, myze
TYPE(FxrInfo_Type), POINTER :: myFxr
INTEGER :: iz, ifxr, id, nlocalfxr, ixy, nxy, FxrIdxst, icel
INTEGER :: i, j, k
REAL :: boronppm
Pin=>Core%Pin; Cell=>Core%CellInfo
nxy=Core%nxy
DO iz = myzb, myze
	DO ixy = 1,nxy
		FxrIdxst=pin(ixy)%FxrIdxst

		icel=pin(ixy)%cell(iz)
		nlocalfxr=cell(icel)%nfxr
		DO j=1, nLocalFxr
			ifxr=FxrIdxSt +j -1
			myFxr => Fxr(ifxr, iz)

			IF(.NOT. myFxr%lH2o) CYCLE

			Boronppm = THInfo%CBMCool(iz,ixy)
			
			CALL  FxrBoronUpdate(myFxr, Boronppm)
		ENDDO
  ENDDO
ENDDO


END SUBROUTINE


SUBROUTINE UpdtBoronCmfdXS(Core, Fxr, PhiS, PinXS, boronppm, myzb, myze, ng)
USE PARAM
USE TypeDef,       ONLY : CoreInfo_type, Fxrinfo_type, PinXS_Type,   &
                          Pin_Type,      Cell_Type
USE NuclidMap_mod, ONLY : AtomicWeight
USE Boron_mod,     ONLY : b10frac
USE MacXsLib_mod,  ONLY : XsMac0,        MacXsBase
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(PinXS_Type), POINTER :: PinXS(:, :)
REAL, POINTER :: phis(:, :, :)

REAL :: boronppm
INTEGER :: myzb, myze, ng 

TYPE(FxrInfo_Type), POINTER :: myFxr
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)

REAL, POINTER :: PinVol(:, :)

INTEGER :: nFxr, nxy
INTEGER :: ig, ig2, igb, ige, iz, ixy, ifxr, ifsr, id, icel, ifsrlocal
INTEGER :: FxrIdxSt, nlocalFxr, nFsrinFxr, FsrIdxSt
INTEGER :: i, j, k
REAL :: rho, aw, tp1,tp2, temp, pnumB, volsum, delpnum
REAL :: DelXs(ng), RR(ng), FxrPhi(ng)
INTEGER :: idiso(2)
REAL :: pnum(2)
LOGICAL :: lBoronIstope, lh2o

nFxr = Core%nCoreFxr
nxy = core%nxy
Pin => Core%Pin; Cell => Core%CellInfo

DO iz = myzb, myze
  DO ixy = 1, nxy
    FxrIdxSt = Pin(ixy)%FxrIdxSt; FsrIdxSt = Pin(ixy)%FsrIdxSt
    icel = Pin(ixy)%Cell(iz)
    nlocalfxr = Cell(icel)%nFxr 
    lh2o = .FALSE.
    volsum = 0; DelXs = 0; RR = 0; temp = 0
    DO j = 1, nlocalFxr  
      ifxr = FxrIdxSt + j - 1
      myFxr => Fxr(ifxr, iz)
      lh2o = lh2o .or. myFxr%lh2o
      volsum = volsum + myFxr%area
      nFsrInFxr = Cell(icel)%nFsrInFxr(j)
      FxrPhi = 0
      DO i = 1, nFsrInFxr
        ifsrlocal = Cell(icel)%MapFxr2FsrIdx(i, j)
        ifsr = FsrIdxSt + ifsrlocal - 1
        FxrPhi(:) = FxrPhi + phis(ifsr, iz, :) * Cell(icel)%vol(ifsrlocal)
      ENDDO
      RR(1:ng) = RR(1:ng) +  FxrPhi(1:ng)
      IF(.not. myFxr%lh2o) cycle
      volsum = volsum + myFxr%area
      temp = temp + myFxr%area * myFxr%temp
      rho = 0
      DO i = 1, myFxr%niso
        id = myFxr%idiso(i)/1000
        IF(id .EQ. 1 .OR. id .EQ. 8) THEN
          aw = AtomicWeight(myFxr%idiso(i))
          rho = rho + myFxr%pnum(i) * aw / AVOGADRO
        ENDIF
      ENDDO
      lBoronIstope = FALSE
      pnumB = rho * boronppm * epsm6 * avogadro / awboron
      DO i = 1, myFxr%niso
        id = myFxr%idiso(i)
        IF(id .EQ. 5000) THEN
          delpnum = (pnumB-myFxr%pnum(i))
          myFxr%pnum(i) = pnumB
          lBoronIstope = TRUE   
        ENDIF
      ENDDO
      IF(.NOT. lBoronIstope) THEN
        myFxr%niso = myFxr%niso + 1; i = myFxr%niso
        myFxr%idiso(i) = 5000
        myFxr%pnum(i) = pnumB
        delpnum = pnumB
      ENDIF
      !Calculate the effective number density chagne
      DO ig = 1, ng
        DelXs(ig) = DelXs(ig) + delpnum * FxrPhi(ig)
      ENDDO
    ENDDO

    DO ig = 1, ng
      DelXs(ig) = DelXs(ig) / RR(ig)
    ENDDO
    IF(.NOT. lh2o) CYCLE
    temp = temp / volsum
    pnum(1) = 1; idiso(1) = 5000
    Call MacXsBase(XsMac0, temp, 1, idiso(1:1), pnum(1:1), 1, ng, ng, 1._8, .FALSE., .TRUE.)
    DO ig = 1, ng
      PinXS(ixy, iz)%XST(ig) = PinXS(ixy, iz)%XST(ig) + DelXS(ig) * XsMac0%xsmact(ig)
      PinXS(ixy, iz)%XSTR(ig) = PinXS(ixy, iz)%XSTR(ig) + DelXS(ig) * XsMac0%xsmactr(ig)
      PinXS(ixy, iz)%XSR(ig) = PinXS(ixy, iz)%XSR(ig) + DelXS(ig) * (XsMac0%xsmactr(ig)-XsMac0%xsmacsm(ig, ig))
      igb = PinXs(ixy,iz)%xss(ig)%ib; ige = PinXs(ixy,iz)%xss(ig)%ie
      DO ig2 = igb, ige
        PinXS(ixy, iz)%xss(ig)%from(ig2) = PinXS(ixy, iz)%xss(ig)%from(ig2) + DelXS(ig) * XsMac0%XsMacSm(ig2, ig)
      ENDDO
      PinXS(ixy, iz)%xss(ig)%withinGroupScat = PinXS(ixy, iz)%xss(ig)%withinGroupScat + PinXS(ixy, iz)%xss(ig)%from(ig)
      PinXS(ixy, iz)%xss(ig)%from(ig) = 0
      PinXs(ixy, iz)%XSD(ig) = 1._8/(3._8*PinXs(ixy, iz)%XSTR(ig))
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE