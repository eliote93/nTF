#include <defines.h>
SUBROUTINE CalB10XS(Core, FmInfo, istep, ng, GroupInfo, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,       ONLY : CoreInfo_Type,    FmInfo_Type,     GroupInfo_Type,     &
                          PE_Type,                                               &
                          FxrInfo_Type,     Pin_Type,        Cell_Type
USE Depl_mod,      ONLY : FluxNormalizeFactor, DeplCntl
USE MacXsLib_Mod,  ONLY : MacXS1gBase, MacXsAF
USE CNTL,          ONLY : nTracerCntl_Type
USE MOC_Mod,       ONLY : FxrAvgPhi
USE Boron_Mod,     ONLY : sigc1g_b10,        phi1g_mod,                          &
                          b10frac0
#ifdef MPI_ENV
USE MPIComm_Mod,   ONLY : REDUCE
#endif
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE
INTEGER :: ng
INTEGER :: istep

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(Pin_Type), POINTER :: Pin(:)
REAL, POINTER :: Phis(:, :, :)
REAL, POINTER :: hz(:)
REAL :: pnum(1)
INTEGER :: idiso(1)

TYPE(FxrInfo_Type), POINTER :: myFxr

REAL :: NormFactor, Phi1g, AvgPhi(ng)
REAL :: xsa(1), xsf(1), xsn2n(1), xsn3n(1), temp
REAL :: isoA(1,ng), isoF(1,ng)
REAL :: xsc1g_b10, phi1g_b10, sumvol, pnumsum
REAL :: buf0(4), buf(4)
INTEGER :: nfxr, nfsr, nxy
INTEGER :: myzb, myze
INTEGER :: FxrIdxst, FsrIdxst, nLocalFxr, nLocalFsr
INTEGER :: iz, ixy, ifxr, icel
INTEGER :: i, j, k


nfxr = Core%nCoreFxr; nFsr = Core%nCoreFsr
myzb = PE%myzb; myze = PE%myze
nxy = Core%nxy
Pin => Core%Pin; CellInfo => Core%CellInfo
Phis => FmInfo%Phis; Fxr => FmInfo%Fxr
hz => Core%hz;

NormFactor = FluxNormalizeFactor(Core, FmInfo, GroupInfo, nTracerCntl%PowerCore, FALSE, TRUE, PE)
NormFactor = NormFactor * nTracerCntl%PowerLevel
xsc1g_b10 = 0; phi1g_b10 = 0; sumvol = 0; pnumsum =0
DO iz = myzb, myze
  DO ixy = 1, nxy
    FsrIdxSt = Pin(ixy)%FsrIdxSt; FxrIdxSt = Pin(ixy)%FxrIdxSt
    icel = Pin(ixy)%Cell(iz);
    nlocalFxr = CellInfo(icel)%nFxr; nlocalFsr = CellInfo(icel)%nFsr
    DO j =1, nLocalFxr
      ifxr = FxrIdxSt + j - 1
      myFxr => Fxr(ifxr, iz)
      IF(.NOT. myFxr%lH2O) CYCLE
      AvgPhi = FxrAvgPhi(Core, Fxr, Phis, ixy, j, iz, ng, PE)
      DO i = 1, myFxr%niso
        IF (myFxr%idiso(i).EQ.5010) THEN
          CALL MacXsAF(isoA,isoF,myFxr%temp,1,myFxr%idiso(i:i),myFxr%pnum(i:i),1,ng,ng)
          CALL MacXS1gBase(xsa, xsf, xsn2n, xsn3n, myFxr%temp, 1, myFxr%idiso(i:i), myFxr%pnum(i:i), AvgPhi(1:ng), ng)
          xsa(1) = SUM(isoA(1,:)*AvgPhi)/SUM(AvgPhi)
          Phi1g = SUM(AvgPhi(1:ng))*hz(iz)*myFxr%area
          xsc1g_b10 = xsc1g_b10 + phi1g*xsa(1)
          !print*, xsc1g_b10, xsa(1)
          phi1g_b10 = phi1g_b10 + phi1g
          pnumsum = pnumsum + myFxr%pnum(i)*hz(iz)*myFxr%area
          sumvol = sumvol + hz(iz)*myFxr%area
        END IF

        !IF(myFxr%idiso(i) .NE. 5000) CYCLE
        !idiso(1) = 5010
        !pnum(1) = myFxr%pnum(i) * b10frac0
        !CALL MacXS1gBase(xsa, xsf, xsn2n, xsn3n, myFxr%Temp, 1, idiso(1:1), pnum(1:1), AvgPhi(1:ng), ng)
        !Phi1g = SUM(AvgPhi(1:ng)) * hz(iz) * myFxr%area
        !xsc1g_b10 = xsc1g_b10 + phi1g * (xsa(1) - xsf(1))
        !phi1g_b10 = phi1g_b10 + phi1g
        !pnumsum = pnumsum +  myFxr%pnum(i) * hz(iz) * myFxr%area
        !sumvol = sumvol + hz(iz) * myFxr%area

        !IF(myFxr%idiso(i) .NE. 5010) CYCLE
        !!For B-10
        !CALL MacXS1gBase(xsa, xsf, xsn2n, xsn3n, myFxr%Temp, 1, myFxr%idiso(i:i), myFxr%pnum(i:i), AvgPhi(1:ng), ng)
        !Phi1g = SUM(AvgPhi(1:ng)) * hz(iz) * myFxr%area
        !xsc1g_b10 = xsc1g_b10 + phi1g * (xsa(1) - xsf(1))
        !phi1g_b10 = phi1g_b10 + phi1g
        !pnumsum = pnumsum +  myFxr%pnum(i) * hz(iz) * myFxr%area
        !sumvol = sumvol + hz(iz) * myFxr%area
      ENDDO
    ENDDO
  ENDDO
ENDDO

#ifdef  MPI_ENV
buf0 = (/xsc1g_b10, phi1g_b10, sumvol, pnumsum/)
buf = 0
CALL REDUCE(buf0, buf, 4, PE%MPI_RTMASTER_COMM, .TRUE.)
xsc1g_b10 = buf(1); phi1g_b10 = buf(2); sumvol = buf(3); pnumsum = buf(4)
#endif

pnumsum = pnumsum / sumvol
xsc1g_b10 = xsc1g_b10 / phi1g_b10/pnumsum
phi1g_b10 = phi1g_b10 / sumvol


phi1g_mod(istep) = phi1g_b10 * 1.0E-24_8 * NormFactor
sigc1g_b10(istep) = xsc1g_b10
!print*, phi1g_mod(istep), sigc1g_b10(istep), istep
NULLIFY(Pin, CellInfo, Phis, Fxr, hz)

END SUBROUTINE

SUBROUTINE PostB10Depl0(nowstep, DeplCntl, PE)
USE PARAM
USE TYPEDEF,    ONLY : PE_TYPE
USE DeplType,   ONLY : DeplCntl_Type
USE FILES,               ONLY : io8
USE Boron_mod,  ONLY : sigc1g_b10,          phi1g_mod,                        &
                       BoronPPM,            DeplBoronPPM,    DeplB10Frac,     &
                       b10frac0!,            vratio
USE IOUTIL,     ONLY : message

IMPLICIT NONE

TYPE(DeplCntl_Type) :: DeplCntl
TYPE(PE_TYPE) :: PE
INTEGER :: nowstep

INTEGER :: nstep
INTEGER :: istep
REAL :: RR1, RR2, RR, vratio
REAL :: time(2),delt, theta, b10, B11
theta = 0.5
nStep = DeplCntl%nBurnupStep
!vratio = 0.053957
vratio = DeplCntl%vratio
time(1) = 0.00001_8 / DeplCntl%PowerCore * DeplCntl%Hm_Mass0_kg
b10 = b10frac0; b11=1._8 - b10
DeplB10Frac(0) = b10; DeplBoronPPM(0) = BoronPPM(0)
DO istep = 1, nowstep
 time(2) = DeplCntl%T_efpd(istep)
 DelT = (time(2) - time(1)) * 86400._8
 ! Only work for PC_OPT = 1 (semi PC)
 ! For full PC, use the online subroutine
 RR1 = sigc1g_b10(istep) * phi1g_mod(istep)
 RR2 = sigc1g_b10(istep-1) * phi1g_mod(istep-1)
 RR = (theta * RR1 + (1-theta) * RR2) * vratio
 b10 = b10 * exp(-DelT * RR)
 b10 = b10/(b10+b11)
 b11 = 1.-b10
 time(1) = time(2)
 DeplB10Frac(istep) = b10
 DeplBoronPPM(istep) = BoronPPM(istep) * b10frac0 / b10
 IF(PE%MASTER .and. istep .eq. nowstep) THEN
   write(88, *) istep,  BoronPPM(istep),  DeplBoronPPM(istep), b10
   write(mesg,'(a,f12.2,a, f12.5)') 'B-10 Depletion Correction :', DeplBoronPPM(istep), ' PPM', b10
   CALL MESSAGE(io8, TRUE, TRUE, mesg)
 ENDIF
ENDDO

END SUBROUTINE


SUBROUTINE PostB10Depl(DeplCntl, PE)
USE PARAM
USE TYPEDEF,    ONLY : PE_TYPE
USE DeplType,   ONLY : DeplCntl_Type
USE Boron_mod,  ONLY : sigc1g_b10,          phi1g_mod,                        &
                       BoronPPM,            DeplBoronPPM,    DeplB10Frac,     &
                       b10frac0!,            vratio

IMPLICIT NONE

TYPE(DeplCntl_Type) :: DeplCntl
TYPE(PE_TYPE) :: PE

INTEGER :: nstep
INTEGER :: istep
REAL :: RR1, RR2, RR, vratio
REAL :: time(2),delt, theta, b10, B11
theta = 0.5
nStep = DeplCntl%nBurnupStep
!vratio = 0.053957
vratio = DeplCntl%vratio
time(1) = 0.00001_8 / DeplCntl%PowerCore * DeplCntl%Hm_Mass0_kg
b10 = b10frac0; b11=1._8 - b10
DeplB10Frac(0) = b10; DeplBoronPPM(0) = BoronPPM(0)
DO istep = 1, nStep
 time(2) = DeplCntl%T_efpd(istep)
 DelT = (time(2) - time(1)) * 86400._8
 RR1 = sigc1g_b10(istep) * phi1g_mod(istep)
 RR2 = sigc1g_b10(istep-1) * phi1g_mod(istep-1)
 RR = (theta * RR1 + (1-theta) * RR2) * vratio
 b10 = b10 * exp(-DelT * RR)
 b10 = b10/(b10+b11)
 b11 = 1.-b10
 time(1) = time(2)
 DeplB10Frac(istep) = b10
 DeplBoronPPM(istep) = BoronPPM(istep) * b10frac0 / b10
 IF(PE%MASTER) write(88, *) istep,  BoronPPM(istep),  DeplBoronPPM(istep), b10
ENDDO

END SUBROUTINE

SUBROUTINE OnlineB10Depl(istep, lpredict, DeplCntl, PE)
USE PARAM
USE TYPEDEF,    ONLY : PE_TYPE
USE DeplType,   ONLY : DeplCntl_Type
USE Boron_mod,  ONLY : sigc1g_b10,          phi1g_mod,                        &
                       BoronPPM,            DeplBoronPPM,    DeplB10Frac,     &
                       b10frac0,            b10frac!,         vratio
USE FILES,               ONLY : io8
USE IOUTIL,              ONLY : message
IMPLICIT NONE

TYPE(DeplCntl_Type) :: DeplCntl
TYPE(PE_TYPE) :: PE
LOGICAL :: lPredict
INTEGER :: istep

INTEGER :: nstep
REAL :: RR1, RR2, RR, vratio
REAL :: time(2),delt, theta, b10, B11

theta = 0.5
IF(lPredict) theta = 0
nStep = DeplCntl%nBurnupStep
!vratio = 0.25_8
vratio = DeplCntl%vratio

IF(istep .EQ. 1) THEN
  time(1) = 0.00001_8 / DeplCntl%PowerCore * DeplCntl%Hm_Mass0_kg
  DeplB10Frac(0) = b10frac0; !DeplBoronPPM(0) = BoronPPM(0)
ELSE
  time(1) = DeplCntl%T_efpd(istep-1)
ENDIF
b10 =  DeplB10Frac(istep-1); b11=1-b10
time(2) = DeplCntl%T_efpd(istep)
DelT = (time(2) - time(1)) * 86400._8
RR1 = sigc1g_b10(istep) * phi1g_mod(istep)
RR2 = sigc1g_b10(istep-1) * phi1g_mod(istep-1)
IF(lPredict) RR1 = 0
RR = (theta * RR1 + (1-theta) * RR2) * vratio
b10 = b10 * exp(-DelT * RR)
b10 = b10/(b10+b11)
b11 = 1.-b10
b10frac = b10
DeplB10Frac(istep) = b10
IF(PE%MASTER) THEN
  write(mesg,'(a,f8.3)') "Boron-10 : ", B10
  CALL MESSAGE(io8, TRUE, TRUE, mesg)
ENDIF
!IF(PE%MASTER) print *, "Boron-10 : ", B10
!time(1) = 0.00001_8 / DeplCntl%CorePower * DeplCntl%Hm_Mass0_kg

END SUBROUTINE
