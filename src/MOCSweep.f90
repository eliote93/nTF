#include <defines.h>
SUBROUTINE MOCSweep(Core, RayInfo, FMInfo, eigv, ng, PE, nTracerCntl, ItrCntl)

USE PARAM,       ONLY : TRUE, FALSE, ZERO, mesg
USE TYPEDEF,     ONLY : CoreInfo_Type, RayInfo_Type, FmInfo_Type, PE_TYPE, FxrInfo_Type
USE CNTL,        ONLY : nTracerCntl_Type
USE itrcntl_mod, ONLY : ItrCntl_TYPE
USE CORE_MOD,    ONLY : GroupInfo, srcSlope, phisSlope, psiSlope
USE MOC_MOD,     ONLY : SetRtMacXsGM, SetRtSrcGM, SetRtLinSrc, SetRtP1SrcGM, AddBucklingGM, RayTraceGM_AFSS, RayTraceGM_OMP, RayTraceP1GM_OMP, RayTraceP1GM_AFSS, &
                        PsiUpdate, CellPsiUpdate, UpdateEigv, MocResidual, PsiErr, PseudoAbsorptionGM, PowerUpdate, FluxUnderRelaxation, FluxInUnderRelaxation, CurrentUnderRelaxation, &
                        phis1g, phim1g, MocJout1g, xst1g, tSrc, AxSrc1g, PhiAngin1g, srcm, &
                        LinPsiUpdate, RayTraceLS_CASMO, RayTraceLS, SetRTLinSrc_CASMO, LinPsiUpdate_CASMO, LinSrc1g, LinPsi, &
                        SetRtMacXsNM, SetRtsrcNM, AddBucklingNM, SetRtP1srcNM, PseudoAbsorptionNM, RayTraceNM_OMP, RayTraceP1NM_OMP, phisNM, PhiAngInNM, MocJoutNM, xstNM, srcNM, phimNM, srcmNM, &
                        RayTrace_Dcmp, RayTraceP1_Dcmp, RayTraceLin_Dcmp, DcmpPhiAngIn, DcmpGatherCurrent
USE SUbGrp_Mod,  ONLY : FxrChiGen
USE IOUTIL,      ONLY : message
USE Timer,       ONLY : nTracer_dclock, TimeChk
USE FILES,       ONLY : io8
USE XSLIB_MOD,   ONLY : igresb, igrese
USE SPH_mod,     ONLY : ssphf, ssphfNM, calcPinSSPH
USE HexData,     ONLY : nInnMOCItr

use SubChCoupling_mod, ONLY : last_TH

#ifdef __INTEL_COMPILER
USE IFPORT,      ONLY : HOSTNM
#endif
#ifdef MPI_ENV
USE MPICOMM_MOD, ONLY : BCAST, MPI_SYNC, MPI_MAX_REAL
#endif

IMPLICIT NONE

TYPE (CoreInfo_Type)    :: Core
TYPE (RayInfo_Type)     :: RayInfo
TYPE (FmInfo_Type)      :: FmInfo
TYPE (PE_TYPE)          :: PE
TYPE (nTracerCntl_Type) :: nTracerCntl
TYPE (ItrCntl_TYPE)     :: ItrCntl

REAL :: eigv
INTEGER :: ng
! ----------------------------------------------------

INTEGER :: ig, iz, ist, ied, iout, jswp, iinn, ninn, nitermax, myzb, myze, nPhiAngSv, nPolarAngle, nginfo, GrpBeg, GrpEnd, nscttod, fmoclv
INTEGER :: grpbndy(2, 2)

REAL :: eigconv, psiconv, resconv, psipsi, psipsid, eigerr, fiserr, peigv, reserr, tmocst, tmoced, t1gst, t1ged, tngst, tnged, tdel
REAL :: errdat(3), tngdel(2)

LOGICAL :: lJout, lxslib, l3dim, lscat1, ltrc, lRST, lAFSS, lssph, lssphreg, MASTER, RTMaster, lDmesg, lLinSrc, lLSCASMO, lmocUR, lbsq, laxrefFDM, ldcmp

CHARACTER(80) :: hostname

#ifdef __PGI
#ifdef __linux
INTEGER :: hostnm
#endif
#endif

TYPE(FxrInfo_Type), POINTER, DIMENSION(:,:) :: Fxr

REAL, POINTER, DIMENSION(:,:)       :: psi, psid, psic
REAL, POINTER, DIMENSION(:,:,:)     :: phis, axsrc, axpxs
REAL, POINTER, DIMENSION(:,:,:,:)   :: linsrcslope, phim
REAL, POINTER, DIMENSION(:,:,:,:,:) :: radjout

! DEBUG
INTEGER :: iFSR
! ----------------------------------------------------

tmocst = nTracer_dclock(FALSE, FALSE)

! Pointing
FXR     => FmInfo%Fxr
phis    => FmInfo%PHIS
psi     => FmInfo%PSI
psid    => FmInfo%PSID
psic    => FmInfo%PSIC
RadJout => FmInfo%RadJout
AxSrc   => FmInfo%AxSrc
AxPXS   => FmInfo%AxPXS

Master   = PE%Master
RTMaster = PE%RTMaster

IF (nTracerCntl%lscat1)  phim        => FmInfo%phim
IF (nTracerCntl%lLinSrc) LinSrcSlope => FmInfo%LinSrcSlope

! sSPH
IF (nTracerCntl%lsSPHreg) THEN
  WRITE (mesg, '(a24)') 'Calculating Pin SSPH...'
  IF (Master) CALL message(io8, TRUE, TRUE, mesg)
  
  CALL calcPinSSPH(Core, Fxr, PE)
END IF

CALL omp_set_num_threads(PE%nThread)

! Basic
myzb    = PE%myzb
myze    = PE%myze
lJout   = TRUE
nscttod = nTracerCntl%scatod

nPolarAngle = RayInfo%nPolarAngle
nPhiAngSv   = RayInfo%nPhiAngSv

nginfo = 2
IF (.NOT. GroupInfo%lUpScat) nginfo = 1

ninn = 2
IF (nTracerCntl%lHex) ninn = nInnMOCItr

grpbndy(1, 1) = 1
grpbndy(2, 1) = ng
grpbndy(1, 2) = GroupInfo%UpScatRange(1)
grpbndy(2, 2) = GroupInfo%UpScatRange(2)

lxslib    = nTracerCntl%lXslib
l3dim     = nTracerCntl%l3dim
lscat1    = nTracerCntl%lscat1
ltrc      = nTracerCntl%lTrCorrection
lRST      = nTracerCntl%lRST
lssph     = nTracerCntl%lsSPH
lssphreg  = nTracerCntl%lsSPHreg
lLinSrc   = nTracerCntl%lLinSrc
lLSCASMO  = nTracerCntl%lLinSrcCASMO
lmocUR    = nTracerCntl%lmocUR
lbsq      = nTracerCntl%lbsq
laxrefFDM = nTracerCntl%laxrefFDM
ldcmp     = nTracerCntl%ldomaindcmp
lAFSS     = nTracerCntl%lAFSS
fmoclv    = nTracerCntl%FastMOCLv

nitermax = itrcntl%MOCItrCntl%nitermax
psiconv  = itrcntl%psiconv
eigconv  = itrcntl%eigconv
resconv  = itrcntl%resconv

lDmesg = TRUE
IF (ng .GT. 10) lDmesg = FALSE

! UPD : psi
IF (RTMASTER) THEN
  psid = psi
  
  CALL PsiUpdate(Core, Fxr, phis, psi, myzb, myze, ng, lxslib, GroupInfo)
  CALL CellPsiUpdate(Core, Psi, psic, myzb, myze)
END IF
! ----------------------------------------------------
IF (.NOT. nTracerCntl%lNodeMajor) THEN
  DO iout = 1, nitermax
    itrcntl%mocit = itrcntl%mocit + 1
    
    tdel = ZERO
    
    IF (Master) THEN
      WRITE (mesg, '(A22, I5, A3)') 'Performing Ray Tracing', itrcntl%mocit, '...'
      CALL message(io8, TRUE, TRUE, mesg)
    END IF
    
    IF (RTMaster) CALL FxrChiGen(Core, Fxr, FmInfo, GroupInfo, PE, nTracerCntl)
    
    IF (lLinSrc .AND. RTMaster) CALL LinPsiUpdate(Core, Fxr, LinPsi, LinSrcSlope, myzb, myze, ng, lxslib, GroupInfo)
    
    DO jswp = 1, nginfo
      GrpBeg = grpbndy(1, jswp)
      GrpEnd = grpbndy(2, jswp)
      
      DO ig = GrpBeg, GrpEnd
        t1gst = nTracer_dclock(FALSE, FALSE)
        
        DO iz = myzb, myze
          IF (.NOT. Core%lFuelPlane(iz) .AND. laxrefFDM) CYCLE
        
          IF (RTMASTER .AND. l3dim) AxSrc1g = AxSrc(:, iz, ig)
                    
          DO iinn = 1, ninn
            ljout = iinn.EQ.ninn .OR. lmocUR
            
            IF (RTMASTER) THEN
              ! SET : XS
              CALL SetRtMacXsGM(Core, Fxr(:, iz), xst1g, iz, ig, ng, lxslib, ltrc, lRST, lssph, lssphreg, PE)
#ifdef LkgSplit
              CALL PseudoAbsorptionGM(Core, Fxr(:, iz), tsrc, phis(:, iz, ig), AxPXS(:, iz, ig), xst1g, iz, ig, ng, GroupInfo, l3dim)
#endif
#ifdef Buckling
              IF (nTracerCntl%lBsq) CALL AddBucklingGM(Core, Fxr, xst1g, nTracerCntl%Bsq, iz, ig, ng, lxslib, lRST)
#endif
              ! SET : Src.
              CALL SetRtSrcGM(Core, Fxr(:, iz), tsrc, phis, psi, AxSrc1g, xst1g, eigv, iz, ig, ng, GroupInfo, l3dim, lXslib, lscat1, FALSE, PE)
              
              PhiAngin1g = FmInfo%PhiAngin(:, :, iz, ig)
              
              IF (lscat1) phim1g = phim(:, :, iz, ig)
              IF (lscat1) CALL SetRtP1SrcGM(Core, Fxr(:, iz), srcm, phim, xst1g, iz, ig, ng, GroupInfo, l3dim, lXsLib, lscat1, nscttod, PE)
            END IF
            
            ! Ray Trace
            IF (.NOT. lLinSrc) THEN
              IF (.NOT. lscat1) THEN
                IF (lAFSS) THEN
                  CALL RayTraceGM_AFSS(RayInfo, Core, phis1g, PhiAngIn1g, xst1g, tsrc, MocJout1g, iz, lJout, fmoclv)
                ELSE
                  CALL RayTraceGM_OMP (RayInfo, Core, phis1g, PhiAngIn1g, xst1g, tsrc, MocJout1g, iz, lJout, fmoclv)
                END IF
              ELSE
                IF (lAFSS) THEN
                  CALL RayTraceP1GM_AFSS(RayInfo, Core, phis1g, phim1g, PhiAngIn1g, xst1g, tsrc, Srcm, MocJout1g, iz, lJout, nscttod, fmoclv)
                ELSE
                  CALL RayTraceP1GM_OMP (RayInfo, Core, phis1g, phim1g, PhiAngIn1g, xst1g, tsrc, Srcm, MocJout1g, iz, lJout, nscttod, fmoclv)
                END IF
              END IF
            ELSE
              CALL LinPsiUpdate(Core, Fxr, LinPsi, LinSrcSlope, myzb, myze, ng, lxslib, GroupInfo)
              CALL SetRtLinSrc(Core, Fxr(:, iz), RayInfo, tsrc, LinSrc1g, LinPsi, LinSrcSlope, xst1g, eigv, iz, ig, ng, GroupInfo, l3dim, lXslib, TRUE)
              CALL RayTraceLS(RayInfo, Core, phis1g, PhiAngIn1g, xst1g, tsrc, LinSrc1g, LinSrcSlope(:, :, iz, ig), MocJout1g, iz, lJout)
            END IF
            
            ! Spectral SPH
            IF (lssph .AND. ig.GE.igresb .AND. ig.LE.igrese) phis1g = phis1g * ssphf(:, iz, ig)
            
            ! CnP
            IF (.NOT. RTMASTER) CYCLE
            
            IF (.NOT. lmocUR) THEN
              phis(:, iz, ig) = phis1g
              
              IF (lScat1) phim(:, :, iz, ig) = phim1g
              
              FMInfo%PhiAngIn(:, :, iz, ig) = PhiAngIn1g
            ELSE
              CALL FluxUnderRelaxation   (Core, Phis1g, Phis, FmInfo%w(ig), iz, ig, PE)
              CALL FluxInUnderRelaxation (Core, PhiANgIn1g, FmInfo%PhiAngIn, FmInfo%w(ig), nPolarAngle, nPhiAngSv, iz, ig, PE)
              CALL CurrentUnderRelaxation(Core, MocJout1g, RadJout, FmInfo%w(ig),iz, ig, PE)
            END IF
            
            IF (lJout .AND. RTMASTER .AND. .NOT. lmocUR) RadJout(:, :, :, iz, ig) = MocJout1g
          END DO
          
          t1ged = nTracer_dclock(FALSE, FALSE)
          tdel  = tdel + t1ged - t1gst
          
          IF (.NOT.lDmesg .AND. MOD(ig, 10).NE.0 .AND. ig.NE.ng) CYCLE
          
          CALL MPI_MAX_REAL(tdel, PE%MPI_RTMASTER_COMM, TRUE)
          
          IF (master) THEN
            WRITE (mesg, '(10X, A, I4, 2X, A, F10.3, 2X, A)') 'Group ', ig, ' finished in ', tdel, 'Sec'
            CALL message(io8, FALSE, TRUE, mesg)
          END IF
          
          tdel = ZERO
        END DO
      END DO
    END DO
  END DO
! ----------------------------------------------------
ELSE
  DO iout = 1, nitermax
    itrcntl%mocit = itrcntl%mocit + 1
    
    tngdel = ZERO
    
    IF (MASTER) THEN
      WRITE (mesg, '(A22, I5, A3)') 'Performing Ray Tracing', itrcntl%mocit, '...'
      CALL message(io8, TRUE, TRUE, mesg)
    END IF
    
    IF (RTMaster) CALL FxrChiGen(Core, Fxr, FmInfo, GroupInfo, PE, nTracerCntl)
    
    IF (lLSCASMO .AND. RTMaster) CALL LinPsiUpdate_CASMO(Core, Fxr, phisSlope, psiSlope, myzb, myze, ng, lxslib, GroupInfo)
    
    DO iz = myzb, myze
      IF (.NOT. Core%lFuelPlane(iz) .AND. laxrefFDM) CYCLE
      
      ! SET : XS
      IF (RTMASTER) THEN
        DO ig = 1, ng
          phisNM(ig, :) = phis(:, iz, ig)
          
          PhiAngInNM(:, ig, :) = FmInfo%PhiAngIn(:, :, iz, ig)
        END DO
        
        IF (lscat1) phimNM => phim(:, :, :, iz)
        
        IF (ldcmp) DcmpPhiAngIn => FMInfo%AsyPhiAngIn(:, :, :, :, :, iz)
        
        CALL SetRtMacXsNM(Core, Fxr(:, iz), xstNM, iz, ng, lxslib, ltrc, lRST, lssph, lssphreg, PE)
#ifdef LkgSplit
        CALL PseudoAbsorptionNM(Core, Fxr(:, iz), AxPXS, xstNM, iz, ng, GroupInfo, l3dim)
#endif
#ifdef Buckling
        IF (lbsq) CALL AddBucklingNM(Core, Fxr, xstNM, nTracerCntl%bsq, iz, ng, lxslib, lRST)
#endif
      END IF
      
      DO jswp = 1, nginfo
        GrpBeg = grpbndy(1, jswp)
        GrpEnd = grpbndy(2, jswp)
        
        tngst = nTracer_dclock(FALSE, FALSE)
        
        DO iinn = 1, ninn
          ljout = iinn .EQ. ninn
          
          ! SET : Src.
          IF (RTMASTER) THEN
            IF (.NOT. lLSCASMO) THEN
              CALL SetRtsrcNM(Core, Fxr(:, iz), srcNM, phisNM, psi, AxSrc, xstNM, eigv, iz, GrpBeg, GrpEnd, ng, GroupInfo, l3dim, lXslib, lscat1, FALSE, PE)
              
              IF (lScat1) CALL SetRtP1srcNM(Core, Fxr(:, iz), srcmNM, phimNM, xstNM, iz, GrpBeg, GrpEnd, ng, GroupInfo, lXsLib, nscttod, PE)
            ELSE
              CALL SetRtLinSrc_CASMO(Core, Fxr, RayInfo, phisNM, phisSlope, srcNM, srcSlope, psi, psiSlope, AxSrc, xstNM, eigv, iz, GrpBeg, GrpEnd, ng, GroupInfo, l3dim, lxslib, lscat1, FALSE)
            END IF
          END IF
          
          ! Ray Trace
          IF (.NOT. ldcmp) THEN
            IF (.NOT. lLSCASMO) THEN
              IF (lscat1) THEN
                CALL RayTraceP1NM_OMP(RayInfo, Core, phisNM, phimNM, PhiAngInNM, xstNM, srcNM, srcmNM, MocJoutNM, iz, GrpBeg, GrpEnd, ljout)
              ELSE
                CALL RayTraceNM_OMP  (RayInfo, Core, phisNM,         PhiAngInNM, xstNM, srcNM,         MocJoutNM, iz, GrpBeg, GrpEnd, ljout)
              END IF
            ELSE
              CALL RayTraceLS_CASMO(RayInfo, Core, phisNM, phisSlope, PhiAngInNM, srcNM, srcSlope, xstNM, MocJoutNM, iz, GrpBeg, GrpEnd, lJout)
            END IF
          ELSE
            IF (lScat1) THEN
              CALL RayTraceP1_Dcmp(RayInfo, Core, phisNM, phimNM, PhiAngInNM, xstNM, srcNM, srcmNM, MocJoutNM, iz, GrpBeg, GrpEnd, lJout)
            ELSE
              IF (.NOT.lLSCASMO) THEN
                CALL RayTrace_Dcmp(RayInfo, Core, phisNM,         PhiAngInNM, xstNM, srcNM,         MocJoutNM, iz, GrpBeg, GrpEnd, lJout)
              ELSE
                CALL RayTraceLin_Dcmp(RayInfo, Core, iz, GrpBeg, GrpEnd, lJout, nTracerCntl%lHybrid)
              END IF
            END IF
          END IF
          
          ! Spectral SPH
          IF (.NOT.lsSPH .OR. GrpBeg.GT.igrese) CYCLE
          
          ist = max(igresb, GrpBeg)
          ied = min(igrese, GrpEnd)
          
          phisNM(ist:ied, :) = phisNM(ist:ied, :) * ssphfNM(ist:ied, :, iz)
        END DO
        
        tnged        = nTracer_dclock(FALSE, FALSE)
        tngdel(jswp) = tngdel(jswp) + (tnged - tngst)
      END DO
      
#ifdef MPI_ENV
      IF (PE%nRTProc .GT. 1) CALL DcmpGatherCurrent(Core, MocJoutNM)
#endif

      ! CnP
      IF (.NOT. RTMASTER) CYCLE
      
      DO ig = 1, ng
        DO iFSR = 1, Core%nCoreFSR
          !IF (phisNM(ig, iFSR) .NE. phisNM(ig, iFSR)) THEN
          !  STOP
          !END IF
        END DO
        
        phis              (:, iz, ig) = phisNM       (ig, :)
        FmInfo%PhiAngIn(:, :, iz, ig) = PhiAngInNM(:, ig, :)
        RadJout     (:, :, :, iz, ig) = MocJoutNM (:, ig, :, :)
      END DO
    END DO
    
#ifdef MPI_ENV
    CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
#endif
    
    DO jswp = 1, nginfo
      tdel = tngdel(jswp)
      
      CALL MPI_MAX_REAL(tdel, PE%MPI_NTRACER_COMM, TRUE)
      
      GrpBeg = grpbndy(1, jswp)
      GrpEnd = grpbndy(2, jswp)
      
      IF (master) THEN
        WRITE (mesg, '(10X, A, I4, 2X, A, I4, 2X, A, F10.3, 2X, A)') 'Group ', GrpBeg, ' to ', GrpEnd, ' finished in ', tdel, 'Sec'
        CALL message(io8, FALSE, TRUE, mesg)
      END IF
    END DO
  END DO
END IF
! ----------------------------------------------------

! UPD : psi
IF (RTMASTER) THEN
  psid = psi
  
  CALL PsiUpdate(Core, Fxr, phis, psi, myzb, myze, ng, lxslib, GroupInfo)
  CALL CellPsiUpdate(Core, Psi, psic, myzb, myze)
END IF

! UDD : eigv
CALL UpdateEigv(Core, psi, psid, eigv, peigv, myzb, myze, PE)

! CAL : Res. Err.
IF (RTmaster) THEN
  fiserr = PsiErr(Core, psi, psid, myzb, myze, PE)
  eigerr = abs((eigv - peigv)) / eigv
  reserr = MocResidual(Core, FmInfo, eigv, GroupInfo, ng, PE, nTracerCntl)
END IF

#ifdef MPI_ENV
errdat = (/fiserr, eigerr, reserr/)

CALL BCAST(errdat, 3, PE%MPI_RT_COMM)

fiserr = errdat(1)
eigerr = errdat(2)
reserr = errdat(3)
#endif

IF (MASTER) THEN
  WRITE (mesg ,'(A5,I7,F15.6, 3(1pE12.3))') 'RT', itrcntl%mocit, eigv, eigerr, fiserr, reserr
  CALL message(io8, TRUE, TRUE, mesg)
END IF

ItrCntl%MocItrCntl%ResErr = ResErr

! UPD : power
IF (RTMASTER) CALL PowerUpdate(Core, Fxr, phis, FmInfo%Power, myzb, myze, ng, lxslib, GroupInfo, PE)

! CHK : Cnv.
IF (fiserr.LT.psiconv .AND. eigerr.LT.eigconv .AND. reserr.LT.resconv) THEN
  itrcntl%lconv = TRUE
  
  IF (nTracerCntl%lFeedBack) last_TH = TRUE
END IF
! ----------------------------------------------------
IF (nTracerCntl%lCusping_MPI .AND. fiserr.LT.ItrCntl%decuspconv) nTracerCntl%lCusping_MPI = FALSE

tmoced = nTracer_dclock(FALSE, FALSE)

TimeChk%MocTime = TimeChk%MocTime + tmoced - tmocst

NULLIFY (FXR)
NULLIFY (psi)
NULLIFY (psid)
NULLIFY (psic)
NULLIFY (phis)
NULLIFY (phim)
NULLIFY (RadJout)
NULLIFY (axsrc)
! ----------------------------------------------------

END SUBROUTINE MOCSweep