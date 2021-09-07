#include <defines.h> 
SUBROUTINE MOCSweep(Core, RayInfo, FMInfo, eigv, ng, PE, nTracerCntl, ItrCntl)

USE PARAM,       ONLY : TRUE, FALSE, ZERO, mesg
USE TYPEDEF,     ONLY : CoreInfo_Type, RayInfo_Type, FmInfo_Type, PE_TYPE, FxrInfo_Type
USE CNTL,        ONLY : nTracerCntl_Type
USE itrcntl_mod, ONLY : ItrCntl_TYPE
USE CORE_MOD,    ONLY : GroupInfo, srcSlope, phisSlope, psiSlope
USE MOC_MOD,     ONLY : SetRtMacXsGM, SetRtSrcGM, SetRtP1SrcGM, AddBucklingGM, PseudoAbsorptionGM, RayTrace_GM, RayTraceP1_GM, phis1g, phim1g, MocJout1g, xst1g, src1g, PhiAngin1g, srcm1g, AxSrc1g, &
                        SetRtMacXsNM, SetRtsrcNM, SetRtP1srcNM, AddBucklingNM, PseudoAbsorptionNM, RayTrace_NM, RayTraceP1_NM, phisNg, phimNg, MocJoutNg, xstNg, srcNg, PhiAngInNg, srcmNg, &
                        PsiUpdate, CellPsiUpdate, UpdateEigv, MocResidual, PsiErr, PowerUpdate, FluxUnderRelaxation, FluxInUnderRelaxation, CurrentUnderRelaxation, &
                        SetRtLinSrc, LinPsiUpdate, RayTraceLS_CASMO, RayTraceLS, SetRTLinSrc_CASMO, LinPsiUpdate_CASMO, LinSrc1g, LinPsi, &
                        RayTraceDcmp_NM, RayTraceDcmp_GM, RayTraceDcmpP1_GM, RayTraceDcmpP1_NM, RayTraceLin_Dcmp, DcmpPhiAngInNg, DcmpPhiAngIn1g, DcmpGatherCurrentNg, DcmpGatherCurrent1g
USE SUbGrp_Mod,  ONLY : FxrChiGen
USE IOUTIL,      ONLY : message
USE Timer,       ONLY : nTracer_dclock, TimeChk
USE FILES,       ONLY : io8
USE XSLIB_MOD,   ONLY : igresb, igrese
USE SPH_mod,     ONLY : ssphf, ssphfNM, calcPinSSPH
USE HexData,     ONLY : nInnMOCItr

USE SubChCoupling_mod, ONLY : last_TH

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
INTEGER :: ig, iz, ist, ied, iout, jswp, iinn, ninn, nitermax, myzb, myze, nPhiAngSv, nPolarAngle, nginfo, GrpBeg, GrpEnd, ScatOd, fmoclv
INTEGER :: grpbndy(2, 2)

REAL :: eigconv, psiconv, resconv, psipsi, psipsid, eigerr, fiserr, peigv, reserr, tmocst, tmoced, tgmst, tgmed, tgmdel, tnmst, tnmed
REAL :: errdat(3), tnmdel(2)

LOGICAL :: lJout, lxslib, l3dim, lscat1, ltrc, lRST, lssph, lssphreg, MASTER, RTMaster, lDmesg, lLinSrc, lLSCASMO, lmocUR, lbsq, laxrefFDM, ldcmp

TYPE(FxrInfo_Type), POINTER, DIMENSION(:,:) :: Fxr

REAL, POINTER, DIMENSION(:)           :: wmoc
REAL, POINTER, DIMENSION(:,:)         :: psi, psid, psic
REAL, POINTER, DIMENSION(:,:,:)       :: phis, axsrc, axpxs
REAL, POINTER, DIMENSION(:,:,:,:)     :: linsrcslope, phim, phiangin
REAL, POINTER, DIMENSION(:,:,:,:,:)   :: radjout
REAL, POINTER, DIMENSION(:,:,:,:,:,:) :: AsyPhiAngIn
! ----------------------------------------------------

tmocst = nTracer_dclock(FALSE, FALSE)

CALL omp_set_num_threads(PE%nThread)

Master   = PE%Master
RTMaster = PE%RTMaster
myzb     = PE%myzb
myze     = PE%myze

! Pointing
FXR      => FmInfo%Fxr
phis     => FmInfo%PHIS
psi      => FmInfo%PSI
psid     => FmInfo%PSID
psic     => FmInfo%PSIC
RadJout  => FmInfo%RadJout
AxSrc    => FmInfo%AxSrc
AxPXS    => FmInfo%AxPXS
phiangin => FmInfo%phiangin

IF (nTracerCntl%lmocUR)      wmoc        => FmInfo%w
IF (nTracerCntl%lscat1)      phim        => FmInfo%phim
IF (nTracerCntl%lLinSrc)     LinSrcSlope => FmInfo%LinSrcSlope
IF (nTracerCntl%lDomainDcmp) AsyPhiAngIn => FmInfo%AsyPhiAngIn

! sSPH
IF (nTracerCntl%lsSPHreg) THEN
  WRITE (mesg, '(a24)') 'Calculating Pin SSPH...'
  IF (Master) CALL message(io8, TRUE, TRUE, mesg)
  
  CALL calcPinSSPH(Core, Fxr, PE)
END IF

! Basic
lJout = TRUE

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
fmoclv    = nTracerCntl%FastMOCLv
ScatOd    = nTracerCntl%scatod

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
    
    tgmdel = ZERO
    
    WRITE (mesg, '(A22, I5, A3)') 'Performing Ray Tracing', itrcntl%mocit, '...'
    IF (Master) CALL message(io8, TRUE, TRUE, mesg)
    
    IF (RTMaster) THEN
      CALL FxrChiGen(Core, Fxr, FmInfo, GroupInfo, PE, nTracerCntl)
      
      IF (lLinSrc) CALL LinPsiUpdate(Core, Fxr, LinPsi, LinSrcSlope, myzb, myze, ng, lxslib, GroupInfo)
    END IF
    
    DO jswp = 1, nginfo
      GrpBeg = grpbndy(1, jswp)
      GrpEnd = grpbndy(2, jswp)
      
      DO ig = GrpBeg, GrpEnd
        tgmst = nTracer_dclock(FALSE, FALSE)
        
        DO iz = myzb, myze
          IF (.NOT. Core%lFuelPlane(iz) .AND. laxrefFDM) CYCLE
          
          ! Pointing
          IF (RTMASTER) THEN
            phis1g     => phis(:, iz, ig)
            PhiAngin1g => PhiAngin(:, :, ig, iz)
            MocJout1g  => RadJout(:, :, :, iz, ig)
            
            IF (l3dim)  AxSrc1g         = AxSrc(:, iz, ig)
            IF (lscat1) phim1g         => phim(:, :, ig, iz)
            IF (lscat1) phimNg         => phim(:, :, :, iz)
            IF (ldcmp)  DcmpPhiAngIn1g => AsyPhiAngIn(:, :, :, :, ig, iz)
          END IF
          
          DO iinn = 1, ninn
            ljout = iinn.EQ.ninn .OR. lmocUR
            
            ! SET : XS, Src.
            IF (RTMASTER) THEN
              CALL SetRtMacXsGM(Core, Fxr(:, iz), xst1g, iz, ig, ng, lxslib, ltrc, lRST, lssph, lssphreg, PE)
#ifdef LkgSplit
              CALL PseudoAbsorptionGM(Core, Fxr(:, iz), phis(:, iz, ig), AxPXS(:, iz, ig), xst1g, iz, ig, ng, GroupInfo, l3dim)
#endif
#ifdef Buckling
              IF (nTracerCntl%lBsq) CALL AddBucklingGM(Core, Fxr, xst1g, nTracerCntl%Bsq, iz, ig, ng, lxslib, lRST)
#endif
              CALL SetRtSrcGM(Core, Fxr(:, iz), src1g, phis, psi, AxSrc1g, xst1g, eigv, iz, ig, ng, GroupInfo, l3dim, lXslib, lscat1, FALSE, PE)
              
              IF (lscat1) CALL SetRtP1SrcGM(Core, Fxr(:, iz), srcm1g, phimNg, xst1g, iz, ig, ng, GroupInfo, l3dim, lXsLib, lscat1, ScatOd, PE)
            END IF
            
            ! Ray Trace
            IF (.NOT. lLinSrc) THEN
              IF (.NOT. ldcmp) THEN
                IF (.NOT. lscat1) THEN
                  CALL RayTrace_GM      (RayInfo, Core, phis1g,         PhiAngIn1g, xst1g, src1g,         MocJout1g, iz, lJout, fmoclv)
                ELSE
                  CALL RayTraceP1_GM    (RayInfo, Core, phis1g, phim1g, PhiAngIn1g, xst1g, src1g, srcm1g, MocJout1g, iz, lJout, fmoclv)
                END IF
              ELSE
                IF (.NOT. lscat1) THEN
                  CALL RayTraceDcmp_GM  (RayInfo, Core, phis1g,         PhiAngIn1g, xst1g, src1g,         MocJout1g, iz, lJout)
                ELSE
                  CALL RayTraceDcmpP1_GM(RayInfo, Core, phis1g, phim1g, PhiAngIn1g, xst1g, src1g, srcm1g, MocJout1g, iz, lJout)
                END IF
              END IF
            ELSE
              CALL LinPsiUpdate(Core, Fxr, LinPsi, LinSrcSlope, myzb, myze, ng, lxslib, GroupInfo)
              CALL SetRtLinSrc(Core, Fxr(:, iz), RayInfo, src1g, LinSrc1g, LinPsi, LinSrcSlope, xst1g, eigv, iz, ig, ng, GroupInfo, l3dim, lXslib, TRUE)
              CALL RayTraceLS(RayInfo, Core, phis1g, PhiAngIn1g, xst1g, src1g, LinSrc1g, LinSrcSlope(:, :, iz, ig), MocJout1g, iz, lJout)
            END IF
            
            ! Spectral SPH
            IF (lssph .AND. ig.GE.igresb .AND. ig.LE.igrese) phis1g = phis1g * ssphf(:, iz, ig)
            
            ! CnP
#ifdef MPI_ENV
            IF (PE%nRTProc .GT. 1) CALL DcmpGatherCurrent1g(Core, MocJout1g)
#endif
            IF (.NOT. RTMASTER) CYCLE
            IF (.NOT. lmocUR)   CYCLE
            
            CALL FluxUnderRelaxation   (Core, Phis1g, Phis, wmoc(ig), iz, ig, PE)
            CALL FluxInUnderRelaxation (Core, PhiANgIn1g, PhiAngIn, wmoc(ig), nPolarAngle, nPhiAngSv, iz, ig, PE)
            CALL CurrentUnderRelaxation(Core, MocJout1g, RadJout, wmoc(ig),iz, ig, PE)
          END DO
        END DO
        
        ! Time
        tgmed  = nTracer_dclock(FALSE, FALSE)
        tgmdel = tgmdel + tgmed - tgmst
        
        IF (.NOT.lDmesg .AND. MOD(ig, 10).NE.0 .AND. ig.NE.ng) CYCLE
        
        CALL MPI_MAX_REAL(tgmdel, PE%MPI_RTMASTER_COMM, TRUE)
        
        WRITE (mesg, '(10X, A, I4, 2X, A, F10.3, 2X, A)') 'Group ', ig, ' finished in ', tgmdel, 'Sec'
        IF (master) CALL message(io8, FALSE, TRUE, mesg)
        
        tgmdel = ZERO
      END DO
    END DO
  END DO
! ----------------------------------------------------
ELSE
  DO iout = 1, nitermax
    itrcntl%mocit = itrcntl%mocit + 1
    
    WRITE (mesg, '(A22, I5, A3)') 'Performing Ray Tracing', itrcntl%mocit, '...'
    IF (MASTER) CALL message(io8, TRUE, TRUE, mesg)
    
    IF (RTMaster) CALL FxrChiGen(Core, Fxr, FmInfo, GroupInfo, PE, nTracerCntl)
    IF (RTMaster .AND. lLSCASMO) CALL LinPsiUpdate_CASMO(Core, Fxr, phisSlope, psiSlope, myzb, myze, ng, lxslib, GroupInfo)
        
    DO iz = myzb, myze
      IF (.NOT. Core%lFuelPlane(iz) .AND. laxrefFDM) CYCLE
      
      ! Pointing & SET : XS
      IF (RTMASTER) THEN
        DO ig = 1, ng
          phisNg(ig, :) = phis(:, iz, ig)
        END DO
        
        PhiAngInNg => PhiAngIn(:, :, :, iz)
        
        IF (lscat1) phimNg         => phim(:, :, :, iz)
        IF (ldcmp)  DcmpPhiAngInNg => AsyPhiAngIn(:, :, :, :, :, iz)
        
        CALL SetRtMacXsNM(Core, Fxr(:, iz), xstNg, iz, ng, lxslib, ltrc, lRST, lssph, lssphreg, PE)
#ifdef LkgSplit
        CALL PseudoAbsorptionNM(Core, Fxr(:, iz), AxPXS, xstNg, iz, ng, GroupInfo, l3dim)
#endif
#ifdef Buckling
        IF (lbsq) CALL AddBucklingNM(Core, Fxr, xstNg, nTracerCntl%bsq, iz, ng, lxslib, lRST)
#endif
      END IF
      
      DO jswp = 1, nginfo
        GrpBeg = grpbndy(1, jswp)
        GrpEnd = grpbndy(2, jswp)
        
        tnmst = nTracer_dclock(FALSE, FALSE)
        
        DO iinn = 1, ninn
          ljout = iinn .EQ. ninn
          
          ! SET : Src.
          IF (RTMASTER) THEN
            IF (.NOT. lLSCASMO) THEN
              CALL SetRtsrcNM(Core, Fxr(:, iz), srcNg, phisNg, psi, AxSrc, xstNg, eigv, iz, GrpBeg, GrpEnd, ng, GroupInfo, l3dim, lXslib, lscat1, FALSE, PE)
              
              IF (lScat1) CALL SetRtP1srcNM(Core, Fxr(:, iz), srcmNg, phimNg, xstNg, iz, GrpBeg, GrpEnd, ng, GroupInfo, lXsLib, ScatOd, PE)
            ELSE
              CALL SetRtLinSrc_CASMO(Core, Fxr, RayInfo, phisNg, phisSlope, srcNg, srcSlope, psi, psiSlope, AxSrc, xstNg, eigv, iz, GrpBeg, GrpEnd, ng, GroupInfo, l3dim, lxslib, lscat1, FALSE)
            END IF
          END IF
          
          ! Ray Trace
          IF (.NOT. ldcmp) THEN
            IF (.NOT. lLSCASMO) THEN
              IF (.NOT. lscat1) THEN
                CALL RayTrace_NM    (RayInfo, Core, phisNg,         PhiAngInNg, xstNg, srcNg,         MocJoutNg, iz, GrpBeg, GrpEnd, ljout)
              ELSE
                CALL RayTraceP1_NM  (RayInfo, Core, phisNg, phimNg, PhiAngInNg, xstNg, srcNg, srcmNg, MocJoutNg, iz, GrpBeg, GrpEnd, ljout)
              END IF
            ELSE
              CALL RayTraceLS_CASMO (RayInfo, Core, phisNg, phisSlope, PhiAngInNg, srcNg, srcSlope, xstNg, MocJoutNg, iz, GrpBeg, GrpEnd, lJout)
            END IF
          ELSE
            IF (.NOT. lscat1) THEN
              IF (.NOT.lLSCASMO) THEN
                CALL RayTraceDcmp_NM(RayInfo, Core, phisNg,         PhiAngInNg, xstNg, srcNg,         MocJoutNg, iz, GrpBeg, GrpEnd, lJout)
              ELSE
                CALL RayTraceLin_Dcmp(RayInfo, Core, iz, GrpBeg, GrpEnd, lJout, nTracerCntl%lHybrid)
              END IF
            ELSE
              CALL RayTraceDcmpP1_NM(RayInfo, Core, phisNg, phimNg, PhiAngInNg, xstNg, srcNg, srcmNg, MocJoutNg, iz, GrpBeg, GrpEnd, lJout)
            END IF
          END IF
          
          ! Spectral SPH
          IF (.NOT.lsSPH .OR. GrpBeg.GT.igrese) CYCLE
          
          ist = max(igresb, GrpBeg)
          ied = min(igrese, GrpEnd)
          
          phisNg(ist:ied, :) = phisNg(ist:ied, :) * ssphfNM(ist:ied, :, iz)
        END DO
        
        tnmed        = nTracer_dclock(FALSE, FALSE)
        tnmdel(jswp) = (tnmed - tnmst)
      END DO
      
      ! CnP
#ifdef MPI_ENV
      IF (PE%nRTProc .GT. 1) CALL DcmpGatherCurrentNg(Core, MocJoutNg)
#endif
      IF (.NOT. RTMASTER) CYCLE
      
      DO ig = 1, ng
        phis(:, iz, ig) = phisNg(ig, :)
        
        RadJout(:, :, :, iz, ig) = MocJoutNg(:, ig, :, :)
      END DO
      
      ! Time
      IF (.NOT. l3dim) THEN
        DO jswp = 1, nginfo
          CALL MPI_MAX_REAL(tnmdel(jswp), PE%MPI_NTRACER_COMM, TRUE)
          
          GrpBeg = grpbndy(1, jswp)
          GrpEnd = grpbndy(2, jswp)
          
          WRITE (mesg, '(10X, A, I4, 2X, A, I4, 2X, A, F10.3, 2X, A)') 'Group ', GrpBeg, ' to ', GrpEnd, ' finished in ', tnmdel(jswp), 'Sec'
          IF (master) CALL message(io8, FALSE, TRUE, mesg)
        END DO
      ELSE
        tgmdel = sum(tnmdel(1:nginfo))
        
        CALL MPI_MAX_REAL(tgmdel, PE%MPI_NTRACER_COMM, TRUE)
        
        WRITE (mesg, '(10X, A5, I4, 2X, A, F10.3, 2X, A)') 'Plane ', iz, ' finished in ', tgmdel, 'Sec'
        IF (master) CALL message(io8, FALSE, TRUE, mesg)
      END IF
    END DO
    
#ifdef MPI_ENV
    CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
#endif
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

WRITE (mesg ,'(A5,I7,F15.6, 3(1pE12.3))') 'RT', itrcntl%mocit, eigv, eigerr, fiserr, reserr
IF (MASTER) CALL message(io8, TRUE, TRUE, mesg)

ItrCntl%MocItrCntl%ResErr = ResErr

! UPD : power
IF (RTMASTER) CALL PowerUpdate(Core, Fxr, phis, FmInfo%Power, myzb, myze, ng, lxslib, GroupInfo, PE)

! CHK : Cnv.
IF (fiserr.LT.psiconv .AND. eigerr.LT.eigconv .AND. reserr.LT.resconv) THEN
  itrcntl%lconv = TRUE
  
  IF (nTracerCntl%lFeedBack) last_TH = TRUE
END IF

! FIN
IF (nTracerCntl%lCusping_MPI .AND. fiserr.LT.ItrCntl%decuspconv) nTracerCntl%lCusping_MPI = FALSE

tmoced = nTracer_dclock(FALSE, FALSE)

TimeChk%MocTime = TimeChk%MocTime + tmoced - tmocst
! ----------------------------------------------------
NULLIFY (FXR)
NULLIFY (phis)
NULLIFY (psi)
NULLIFY (psid)
NULLIFY (psic)
NULLIFY (RadJout)
NULLIFY (AxSrc)
NULLIFY (AxPxs)
NULLIFY (phiangin)

IF (lmocur)  NULLIFY (wmoc)
IF (lscat1)  NULLIFY (phim)
IF (lLinSrc) NULLIFY (LinSrcSlope)
IF (ldcmp)   NULLIFY (AsyPhiAngIn)

IF (RTMASTER) THEN
  IF (.NOT. nTracerCntl%lNodeMajor) THEN
    NULLIFY (phis1g)
    NULLIFY (PhiAngin1g)
    NULLIFY (MocJout1g)
    
    IF (lscat1) NULLIFY (phim1g)
    IF (ldcmp)  NULLIFY (DcmpPhiAngIn1g)
  ELSE
    NULLIFY (PhiAngInNg)
    
    IF (lscat1) NULLIFY (phimNg)
    IF (ldcmp)  NULLIFY (DcmpPhiAngInNg)
  END IF
END IF
! ----------------------------------------------------

END SUBROUTINE MOCSweep