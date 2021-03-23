#include <defines.h>
#ifdef __GAMMA_TRANSPORT
!--- CNJ Edit : Driver Routine for Gamma MOC Calculation
SUBROUTINE GammaMOCSweep(Core, RayInfo, FmInfo, PE, nTracerCntl, ItrCntl)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,      RayInfo_Type,       FmInfo_Type,        &
                           PE_TYPE,            FxrInfo_Type
USE CNTL,           ONLY : nTracerCntl_Type
USE itrcntl_mod,    ONLY : ItrCntl_TYPE
USE CORE_MOD,       ONLY : GroupInfo
USE GammaCore_mod
USE GamMOC_MOD
USE IOUTIL,         ONLY : message
USE Timer,          ONLY : nTracer_dclock,     TimeChk
USE FILES,          ONLY : io8
#ifdef MPI_ENV
USE MPICOMM_MOD,    ONLY : BCAST,              MPI_SYNC
#endif
USE HexData,        ONLY : nInnMOCItr
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(RayInfo_Type) :: RayInfo
TYPE(FmInfo_Type) :: FmInfo
TYPE(PE_TYPE) :: PE
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(ItrCntl_TYPE) :: ItrCntl

REAL, POINTER :: phis(:, :, :), phisnm(:, :)
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
INTEGER :: ig, iz, i, l, jsweep, InIter, nInIter
INTEGER :: ng, ngg, nGroupInfo, nFsr
INTEGER :: myzb, myze
INTEGER :: GrpBeg, GrpEnd
INTEGER, SAVE :: iter = 0
REAL :: resconv, reserr
REAL :: TimeRt1gBeg, TimeRt1gEnd, TimeRtngBeg, TimeRtngEnd, TimeRtElapsed
REAL :: MocTimeBeg, MocTimeEnd
LOGICAL :: lJout, lxslib, l3dim, lscat1, lTrCorrection
LOGICAL :: MASTER, RTMASTER

MOCTimeBeg = nTracer_dclock(FALSE, FALSE)

Fxr => FmInfo%Fxr
phis => FmInfo%phis
Master = PE%Master
RTMaster = PE%RTMaster
ng = GroupInfo%ng
ngg = GroupInfo%ngg
nFsr = Core%nCoreFsr
myzb = PE%myzb; myze = PE%myze
nGroupInfo = 2; nInIter = 1
l3dim = nTracerCntl%l3dim
lscat1 = nTracerCntl%lGammaScat1
lTrCorrection = .NOT. lscat1
resconv = itrcntl%resconv
IF (.NOT. GroupInfo%lUpScat_PH) nGroupInfo = 1

IF (.NOT. lscat1) THEN
  ALLOCATE(phisnm(ng, nFsr))
ENDIF

IF (nTracerCntl%lHex) nInIter = nInnMOCItr
iter = iter + 1

IF (lscat1) THEN

  TimeRtElapsed = 0.
  DO jsweep = 1, nGroupInfo
    GrpBeg = 1; GrpEnd = ngg
    IF(jsweep .GT. 1) THEN
      GrpBeg = GroupInfo%UpScatRange_PH(1); GrpEnd = GroupInfo%UpScatRange_PH(2)
    ENDIF
    DO ig = GrpBeg, GrpEnd
      TimeRt1gBeg = nTracer_dclock(FALSE, FALSE)
      DO iz = myzb, myze
        ljout = FALSE
        gphis1g => gphis(:, ig, iz)
        gphim1g => gphim(:, :, ig, iz)
        gPhiAngIn1g => gPhiAngIn(:, :, ig, iz)
        gJout1g => gJout(:, :, :, ig, iz)
        gAxSrc1g => gAxSrc(:, ig, iz)
        gAxPXS1g => gAxPXS(:, ig, iz)
        DO InIter = 1, nInIter
          IF (InIter .EQ. nInIter) ljout = TRUE
          IF (RTMASTER) THEN
            CALL SetGamMacXs(Core, Fxr(:, iz), gxst1g, iz, ig, ngg, lTrCorrection, PE)
#ifdef LkgSplit
            CALL GamPseudoAbsorption(Core, Fxr(:, iz), gAxPXS1g, gxst1g, iz, l3dim)
#endif            
            CALL SetGamSrc(Core, Fxr(:, iz), gsrc1g, phis, gphis, gAxSrc1g, gxst1g, iz, ig, ng, ngg,                   &
                           l3dim, lscat1, PE, GroupInfo)
            CALL SetGamP1Src(Core, Fxr(:, iz), gsrcm1g, gphim, gxst1g, iz, ig, ngg, lscat1, nTracerCntl%GammaScatOd, PE)
          ENDIF
          CALL RayTraceGamma_Pn(RayInfo, Core, gphis1g, gphim1g, gPhiAngIn1g, gxst1g, gsrc1g, gsrcm1g, gJout1g,       &
                                iz, nTracerCntl%GammaScatOd, lJout)
        ENDDO
      ENDDO
#ifdef MPI_ENV
      CALL MPI_SYNC(PE%MPI_RTMASTER_COMM)
#endif
      TimeRt1gEnd = nTracer_dclock(FALSE, FALSE)
      TimeRtElapsed = TimeRtElapsed + TimeRt1gEnd - TimeRt1gBeg
    ENDDO
  ENDDO

ELSE

  TimeRtElapsed = 0.
  DO iz = myzb, myze
    IF (RTMASTER) THEN
      DO ig = 1, ng
        phisnm(ig, :) = phis(:, iz, ig)
      ENDDO
      gphisnm => gphis(:, :, iz)
      gPhiAngInnm => gPhiAngIn(:, :, :, iz)
      gJoutnm => gJout(:, :, :, :, iz)
      CALL SetGamMacXsNM(Core, Fxr(:, iz), gxstnm, iz, ngg, lTrCorrection, PE)
#ifdef LkgSplit
      CALL GamPseudoAbsorptionNM(Core, Fxr(:, iz), gAxPXS, gxstnm, iz, ngg, l3dim)
#endif       
    ENDIF
    DO jsweep = 1, nGroupInfo
      GrpBeg = 1; GrpEnd = ngg
      IF(jsweep .GT. 1) THEN
        GrpBeg = GroupInfo%UpScatRange_PH(1); GrpEnd = GroupInfo%UpScatRange_PH(2)
      ENDIF
      TimeRtngBeg = nTracer_dclock(FALSE, FALSE)
      DO InIter = 1, nInIter
        ljout = FALSE
        IF (InIter .EQ. nInIter) ljout = TRUE
        IF (RTMASTER) THEN
          CALL SetGamSrcNM(Core, Fxr(:, iz), gsrcnm, phisnm, gphisnm, gAxSrc, gxstnm, iz,                              &
                           GrpBeg, GrpEnd, ng, ngg, l3dim, lscat1, PE)
        ENDIF
        CALL RayTraceGamma(RayInfo, Core, gphisnm, gPhiAngInnm, gxstnm, gsrcnm, gJoutnm,                              &
                           iz, GrpBeg, GrpEnd, ljout)
      ENDDO
      TimeRtngEnd = nTracer_dclock(FALSE, FALSE)
      TimeRtElapsed = TimeRtElapsed + (TimeRtngEnd - TimeRtngBeg)
    ENDDO
#ifdef MPI_ENV
    CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
#endif
  ENDDO

ENDIF

IF (RTMASTER) reserr = GammaMocResidual(Core, FmInfo, GroupInfo, PE, nTracerCntl)

WRITE(mesg, '(2x, a10, i2, a14, 1p, e10.4)') 'Iteration ', iter, ' / Residual = ', reserr
IF(Master) CALL message(io8, TRUE, TRUE, mesg)

ItrCntl%lGammaConv = reserr .LT. resconv

#ifdef MPI_ENV
CALL BCAST(reserr, PE%MPI_RT_COMM, 0)
#endif

IF (RTMaster) THEN
 CALL NeutronLocalQUpdate(Core, Fxr, phis, LocalNPower,  GPowerGen, myzb, myze, ng, PE) !-- JSU EDIT 2017.09.15.
 CALL GammaPowerUpdate(Core, Fxr, gphis, gPower, myzb, myze, ngg, PE)
 CALL CompensateGPower(Core,Fxr,GPowerGen,GPower, myzb, myze,PE)
 print*, 'compensategpow'
END IF

MocTimeEnd = nTracer_dclock(FALSE, FALSE)
TimeRtElapsed = MocTimeEnd - MocTimeBeg
TimeChk%MocTime = TimeChk%MocTime + TimeRtElapsed

IF (.NOT. lscat1) THEN
  DEALLOCATE(phisnm)
ENDIF

!#define __PMH_MONITOR
!#ifdef __PMH_MONITOR
 print*, 'compensategpow1'
IF(ItrCntl%lGammaConv) THEN 
 print*, 'compensategpow2'
  CALL PMH_MONITOR(Core, PE, nTracerCntl)
END IF
!#endif

END SUBROUTINE

!#ifdef __PMH_MONITOR
SUBROUTINE PMH_MONITOR(Core, PE, nTracerCntl)
USE TYPEDEF,    ONLY : coreinfo_type, PE_TYPE, &
                       pin_type,  cell_type,  THCell_Type,  Asy_Type, &
                       ASYINFO_TYPE
USE cntl, ONLY : nTracerCntl_Type
USE GammaCore_mod
USE th_mod, ONLY : ThOpt
IMPLICIT NONE 
TYPE(coreinfo_type) :: Core
TYPE(PE_TYPE) :: PE
TYPE(nTracerCntl_Type) :: nTracerCntl

REAL, POINTER :: pinpw(:,:), pingpw(:,:)
REAL, POINTER :: gratio_pin(:,:), gratio_asy(:,:), gratio_core(:,:)
TYPE(pin_type), POINTER :: pin(:)
TYPE(cell_type), POINTER :: cell(:)
TYPE(Asy_Type), POINTER :: asy(:)
TYPE(ASYINFO_TYPE), POINTER :: asyinfo(:)
TYPE(THCell_Type), POINTER :: Thcell
REAL :: locasypw, locplnpw, volsum, locvol
INTEGER :: myzb, myze, nxy
INTEGER :: ixy, iz, icel, iasy, iasytyp, ixy0
INTEGER :: fxrIdxSt, fsrIdxst
INTEGER :: nlocalfxr, nFsrInFxr, locfsr
INTEGER :: ifxr, ifsr
INTEGER :: i, j

myzb = PE%myzb
myze = PE%myze
nxy = Core%nxy

print*, 'a0'
asy => core%Asy
asyinfo => core%AsyInfo
pin => Core%Pin
cell => Core%CellInfo

print*, 'a1'
ALLOCATE(pinpw(nxy, myzb:myze))
ALLOCATE(pingpw(nxy, myzb:myze))

pinpw = 0.
pingpw = 0.

DO iz = myzb, myze
  DO ixy = 1, nxy
CALL SetThCell(cell(Pin(ixy)%cell(iz)), thopt, nTracerCntl)
  END DO 
END DO 

DO iz = myzb, myze
  DO ixy = 1, nxy
    icel = Pin(ixy)%cell(iz)
    thcell => cell(icel)%THCell
    fxridxst = pin(ixy)%FxrIdxSt
    fsrIdxst = pin(ixy)%FsrIdxSt
    nlocalfxr = cell(icel)%nFXR
    volsum = 0.
    DO i = 1, nlocalfxr
      ifxr = fxrIdxSt + i - 1
      nFsrInFxr = Cell(icel)%nFsrInFxr(i)
      DO j = 1, nFsrInFxr
        locfsr = cell(icel)%MapFxr2FsrIdx(j,i)
        locvol = cell(icel)%vol(locfsr)
        ifsr = fsrIdxst + cell(icel)%MapFxr2FsrIdx(j,i) - 1
        pinpw(ixy,iz) = pinpw(ixy,iz) + (LocalNPower(ifsr, iz) + gPower(ifsr,iz))*locvol
        IF((i .LT. ThCell%FuelReg(1)) .OR. (i .GT. Thcell%FuelReg(2))) THEN 
          pingpw(ixy,iz) = pingpw(ixy,iz) + (localNPower(ifsr,iz) + gPower(ifsr, iz)) * locvol
        END IF
      END DO 
    END DO 
  END DO
END DO 

ALLOCATE(gratio_pin(nxy, myzb:myze))
ALLOCATE(gratio_asy(nxy, myzb:myze))
ALLOCATE(gratio_core(nxy, myzb:myze))
DO iz = myzb, myze
  DO ixy = 1, nxy
      IF(.NOT. cell(pin(ixy)%cell(iz))%lfuel) THEN 
        gratio_pin(ixy,iz) = 0.
      ELSE
        gratio_pin(ixy,iz) = pingpw(ixy,iz) / pinpw(ixy,iz)
      END IF
  END DO 
END DO

DO iz = myzb, myze
  DO iasy = 1, core%nxya
    iasytyp = asy(iasy)%AsyType
    locasypw = 0.
    volsum = 0. 
    DO ixy0 = 1, asyinfo(iasytyp)%nxy
      ixy = asy(iasy)%GlobalPinIdx(ixy0)
      IF(.NOT. cell(pin(ixy)%cell(iz))%lfuel) CYCLE 
      locasypw = locasypw + pinpw(ixy,iz) 
      volsum = volsum + core%pinvol(ixy,iz)
    END DO
    DO ixy0 = 1, asyinfo(iasytyp)%nxy
      ixy = asy(iasy)%GlobalPinIdx(ixy0)
      IF(.NOT. cell(pin(ixy)%cell(iz))%lfuel) THEN 
        gratio_asy(ixy,iz) = 0.
      ELSE
        gratio_asy(ixy, iz) = pingpw(ixy,iz) * volsum / (locasypw * core%pinvol(ixy,iz)) 
      END IF
    END DO
  END DO 
END DO 

DO iz = myzb, myze 
  locplnpw = 0.
  volsum = 0.
  DO ixy = 1, nxy
    IF(.NOT. cell(pin(ixy)%cell(iz))%lfuel) CYCLE 
    locplnpw = locplnpw + pinpw(ixy,iz) 
    volsum = volsum + core%pinvol(ixy,iz)
  END DO 
  DO ixy = 1, nxy
      IF(.NOT. cell(pin(ixy)%cell(iz))%lfuel) THEN 
        gratio_core(ixy,iz) = 0.
      ELSE
        gratio_core(ixy,iz) = pingpw(ixy,iz) * volsum / locplnpw / core%pinvol(ixy,iz)
      END IF
  ENDDO 
END DO 

open(118, FILE = 'geom.pin', status='replace')
DO iz = myzb, myze
  DO ixy = 1, nxy
    icel = pin(ixy)%cell(iz)
    write(118, '(3i5,4es14.6)') pin(ixy)%ix, pin(ixy)%iy, iz, cell(icel)%geom%lx, cell(icel)%geom%ly, cell(icel)%geom%cx, cell(icel)%geom%cy
  END DO 
END DO 
close(118)

open(118, FILE = 'gratio.pin', status='replace')
DO iz = myzb, myze
  DO ixy = 1, nxy
    write(118, '(3i5,es14.6)') pin(ixy)%ix, pin(ixy)%iy, iz, gratio_pin(ixy,iz) 
  END DO 
END DO 
close(118)
open(118, FILE = 'gratio.asy', status='replace')
DO iz = myzb, myze
  DO ixy = 1, nxy
    write(118, '(3i5,es14.6)') pin(ixy)%ix, pin(ixy)%iy, iz, gratio_asy(ixy,iz)
  END DO 
END DO 
close(118)
open(118, FILE = 'gratio.core', status='replace')
DO iz = myzb, myze
  DO ixy = 1, nxy
    write(118, '(3i5,es14.6)') pin(ixy)%ix, pin(ixy)%iy, iz, gratio_core(ixy,iz)
  END DO 
END DO 
close(118)

END SUBROUTINE 
!#endif
#endif
