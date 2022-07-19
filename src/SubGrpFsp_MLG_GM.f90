MODULE SubGrpFspMLG_GM

IMPLICIT NONE

CONTAINS

! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SubGrpFsp_MLG_GM(Core, Fxr, THInfo, RayInfo, GroupInfo, nTracerCntl, PE)

USE allocs
USE PARAM,      ONLY : FALSE, ONE, ZERO, EPSM3, TRUE, VoidCell, mesg
USE TYPEDEF,    ONLY : coreinfo_type, Fxrinfo_type, GroupInfo_Type, RayInfo_Type, PE_TYPE, THInfo_Type
USE CNTL,       ONLY : nTracerCntl_Type
USE FILES,      ONLY : IO8
USE SUBGRP_MOD, ONLY : SetPlnLsigP_MLG, SetPlnLsigP_1gMLG, SubGrpFspErr, EquipXSGen_MLG, EquipXsGen_1gMLG, SetSubGrpSrc1g, UpdtFnAdj, UpdtFtAdj
USE MOC_MOD,    ONLY : RayTrace_GM, RayTraceDcmp_GM, DcmpPhiAngIn1g
USE IOUTIL,     ONLY : message
USE TIMER,      ONLY : nTracer_dclock, TimeChk
USE XSLib_mod,  ONLY : mlgdata, mlgdata0

#ifdef MPI_ENV
USE MPICOMM_MOD, ONLY : MPI_SYNC, MPI_MAX_REAL, MPI_MAX_INT, BCAST, REDUCE
#endif

IMPLICIT NONE

TYPE(CoreInfo_Type)    :: Core
TYPE(THInfo_Type)      :: THInfo
TYPE(RayInfo_Type)     :: RayInfo
TYPE(GroupInfo_Type)   :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE)          :: PE

TYPE(FxrInfo_Type), POINTER, DIMENSION(:,:) :: Fxr

INTEGER :: nfsr, nfxr, myzb, myze, iz, ilv, nlv, myitersum
INTEGER :: ig, ig1, ig2, nofg, norg, niter, iter, itermax, itersum, nPhiAngSv, nPolar

REAL :: lv, errmax, errmaxlv, Tbeg, Tend, rt1, rt2, rtnet

REAL, POINTER, DIMENSION(:)     :: phis1g, phis1gd, siglamPot, xstr1g, src1g
REAL, POINTER, DIMENSION(:,:)   :: PhiAngIn1g
REAL, POINTER, DIMENSION(:,:,:) :: jout1g

LOGICAL :: lCLD, lAIC, master, RTmaster, lDcpl, lSilent, ldcmp
! ----------------------------------------------------

nofg = GroupInfo%nofg
norg = GroupInfo%norg

ig1 = nofg + 1
ig2 = nofg + norg

nFxr = Core%nCoreFxr
nFsr = Core%nCoreFsr

myzb     = PE%myzb
myze     = PE%myze
master   = PE%master
RTmaster = PE%RTmaster

nPolar    = RayInfo%nPolarAngle
nPhiAngSv = RayInfo%nPhiAngSv

ldcmp   = nTracerCntl%lDomainDcmp
lSilent = FALSE

#ifdef MPI_ENV
lDcpl   = nTracerCntl%lDcpl
lSilent = lDcpl
#endif

Tbeg = nTracer_Dclock(FALSE, FALSE)

itermax = 100
IF (.NOT.ldcmp .AND. any(Core%RadBC(1:4) .EQ. VoidCell)) itermax = 1

CALL dmalloc(phis1g,  nFsr)
CALL dmalloc(phis1gd, nFsr)
CALL dmalloc(xstr1g,  nFsr)
CALL dmalloc(src1g,   nFsr)
CALL dmalloc(PhiAngIn1g, nPolar, nPhiAngSv)
CALL dmalloc(SigLamPot, nFxr)

IF (ldcmp) CALL dmalloc(DcmpPhiAngIn1g, nPolar, 2, RayInfo%nModRay, Core%nxya)

WRITE (mesg,'(A, F10.2, A)') "Reference Fuel Temperature", THInfo%RefFuelTemp(0), " C"
IF (master) CALL message(io8,TRUE,TRUE,mesg)

WRITE (mesg,'(A)') 'Solving Subgroup FSP (MLG)'
IF (master) CALL message(io8, TRUE, TRUE, mesg)
! ----------------------------------------------------
itersum   = 0
myitersum = 0

DO iz = myzb, myze
  rtnet = 0
  nlv   = mlgdata(iz)%f_nmaclv
  
  phis1g     = ONE
  PhiAngIn1g = ONE
  PhiAngIn1g(:, 1) = ZERO
  
  IF (Core%lFuelPlane(iz)) THEN
    DO ig = ig1, ig2
      errmax = 0._8
      niter  = 0
      
      CALL UpdtFnAdj(Core, Fxr, ig, iz, PE) ! NUMBER DENSITY CONSIDERATION FACTOR
      
      DO ilv = 1, nlv
        lv = mlgdata(iz)%f_maclv(ilv)
        
        ! Transport xs and source setup for Subgroup Fixed Source Problem
        IF (RTMASTER) THEN
          CALL SetPlnLsigP_MLG(SigLamPot, xstr1g, lv, Core, Fxr, ilv, iz, ig, PE)
          CALL SetSubGrpSrc1g(src1g, SigLamPot, xstr1g, Core, Fxr, iz, PE) ! SRC1G(ifsr) = SigLampot(ifxr) / xstr1g(ifsr)
        END IF
        
        DO iter = 1, itermax
          ! Solving FSP
          IF (RtMaster) CALL CopyFlux(phis1g, phis1gd, nFsr) ! Save Previous data
          
          rt1   = nTracer_dclock(FALSE, FALSE)
          IF (ldcmp) THEN
            CALL RayTraceDcmp_GM(RayInfo, Core, phis1g, PhiAngIn1g, xstr1g, src1g, jout1g, iz, FALSE)
          ELSE
            CALL RayTrace_GM    (RayInfo, Core, phis1g, PhiAngIn1g, xstr1g, src1g, jout1g, iz, FALSE, nTracerCntl%FastMocLv)
          END IF
          rt2   = nTracer_dclock(FALSE, FALSE)
          rtnet = rtnet + (rt2 - rt1)
          
          ! Update Equivalence XS
          errmaxlv = SubGrpFspErr(phis1g, phis1gd, nfsr, PE)
          niter    = niter + 1
          
          IF (RTMASTER) THEN
            CALL EquipXSGen_MLG(Phis1g, SigLamPot, xstr1g, ilv, Core, Fxr, iz, ig, PE) ! equivalence XS calculation with flux from FSP
            CALL UpdtFtAdj(Core, Fxr, ilv, ig, iz, PE) ! TCF update
          END IF
          
          IF (errmaxlv .LT. EPSM3) EXIT
          
          IF (RTMASTER) THEN ! with updated TCF
            CALL SetPlnLsigP_MLG(SigLamPot, xstr1g, lv, Core, Fxr, ilv, iz, ig, PE)
            CALL SetSubGrpSrc1g(src1g, SigLamPot, xstr1g, Core, Fxr, iz, PE)
          END IF
        END DO
        
        errmax = max(errmax,errmaxlv)
        TimeChk%SubGrpFSPNum = TimeChk%SubGrpFSPNum + 1
      END DO
      
      WRITE (mesg,'(2X, 2(A, I3), A, I5, 1P, E13.3)') '[Fuel] Pln', iz, ' Group ', ig,'  In_itr', niter, errmax
      IF (master .AND. .NOT.lSilent) CALL MESSAGE(io8, TRUE, TRUE, mesg)
      
      myitersum = myitersum + niter
    END DO
  END IF
  
  IF (Core%lcladPlane(iz)) THEN
    lCLD = TRUE
    lAIC = FALSE
    nlv  = mlgdata0%c_nmaclv1G

    errmax = 0._8
    niter = 0
    
    DO ilv = 1, nlv
      lv = mlgdata0%c_maclv1G(ilv)
      
      ! transport xs and source setup for Subgroup Fixed Source Problem
      IF (RTMASTER) THEN
        CALL SetPlnLsigP_1gMLG(SigLamPot, xstr1g, lv, Core, Fxr, iz, lCLD, lAIC, PE)
        CALL SetSubGrpSrc1g(src1g, SigLamPot, xstr1g, Core, Fxr, iz, PE)
      END IF
      
      DO iter = 1, itermax
        ! Solving FSP
        IF (RtMaster) CALL CopyFlux(phis1g, phis1gd, nFsr) ! Save Previous data
                
        rt1   = nTracer_dclock(FALSE, FALSE)
        IF (ldcmp) THEN
          CALL RayTraceDcmp_GM(RayInfo, Core, phis1g, PhiAngIn1g, xstr1g, src1g, jout1g, iz, FALSE)
        ELSE
          CALL RayTrace_GM    (RayInfo, Core, phis1g, PhiAngIn1g, xstr1g, src1g, jout1g, iz, FALSE, nTracerCntl%FastMocLv)
        END IF
        rt2   = nTracer_dclock(FALSE, FALSE)
        rtnet = rtnet + (rt2 - rt1)
        
        ! Update Equivalence XS
        errmaxlv = SubGrpFspErr(phis1g, phis1gd, nfsr, PE)
        niter    = niter + 1
        
        IF (errmaxlv .LT. EPSM3) EXIT
      END DO
      
      IF (RTMASTER) CALL EquipXSGen_1gMLG(Phis1g, SigLamPot, xstr1g, ilv, Core, Fxr, iz, lCLD, lAIC, PE)
      
      errmax = max(errmax, errmaxlv)
      TimeChk%SubGrpFSPNum = TimeChk%SubGrpFSPNum + 1
    END DO
    
    WRITE (mesg,'(2X, A, I3, A, I5, 1P, E13.3)') '[Clad] Pln', iz, '  In_itr', niter, errmax
    IF (master .AND. .NOT.lSilent) CALL MESSAGE(io8, TRUE, TRUE, mesg)
    
    myitersum = myitersum + niter
  END IF
  
  IF (Core%lAICPlane(iz)) THEN
    nlv    = mlgdata(iz)%f_nmaclv1G
    lCLD   = FALSE
    lAIC   = TRUE
    errmax = 0._8
    niter  = 0
    
    DO ilv = 1, nlv
      lv = mlgdata(iz)%f_maclv1G(ilv)
      
      ! transport xs and source setup for Subgroup Fixed Source Problem
      IF (RTMASTER) THEN
        CALL SetPlnLsigP_1gMLG(SigLamPot, xstr1g, lv, Core, Fxr, iz, lCLD, lAIC, PE)
        CALL SetSubGrpSrc1g(src1g, SigLamPot, xstr1g, Core, Fxr, iz, PE)
      END IF
      
      DO iter = 1, itermax
        ! Solving FSP
        IF (RtMaster) CALL CopyFlux(phis1g, phis1gd, nFsr) ! Save Previous data
                
        rt1   = nTracer_dclock(FALSE, FALSE)
        IF (ldcmp) THEN
          CALL RayTraceDcmp_GM(RayInfo, Core, phis1g, PhiAngIn1g, xstr1g, src1g, jout1g, iz, FALSE)
        ELSE
          CALL RayTrace_GM    (RayInfo, Core, phis1g, PhiAngIn1g, xstr1g, src1g, jout1g, iz, FALSE, nTracerCntl%FastMocLv)
        END IF
        rt2   = nTracer_dclock(FALSE, FALSE)
        rtnet = rtnet + (rt2 - rt1)
        
        ! Update Equivalence XS
        errmaxlv = SubGrpFspErr(phis1g, phis1gd, nfsr, PE)
        niter    = niter + 1
        
        IF (errmaxlv .LT. EPSM3) EXIT
      END DO
      
      IF (RTMASTER) CALL EquipXSGen_1gMLG(Phis1g, SigLamPot, xstr1g, ilv, Core, Fxr, iz, lCLD, lAIC, PE)
      
      errmax = max(errmax,errmaxlv)
      TimeChk%SubGrpFSPNum = TimeChk%SubGrpFSPNum + 1
    END DO
    
    WRITE (mesg,'(2X, A, I3, A, I5, 1P, E13.3)') '[AIC]  Pln', iz, '  In_itr', niter, errmax
    IF (master .AND. .NOT.lSilent) CALL MESSAGE(io8, TRUE, TRUE, mesg)
    
    myitersum = myitersum + niter
  END IF
END DO
! ----------------------------------------------------
DEALLOCATE (xstr1g, SigLamPot, phis1g, phis1gd, PhiAngIn1g)
IF (ldcmp) DEALLOCATE (DcmpPhiAngIn1g)

nTracerCntl%lSubGrpSweep = TRUE
itersum = myitersum

#ifdef MPI_ENV
CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
CALL REDUCE(myitersum, itersum, PE%MPI_NTRACER_COMM, FALSE)
#endif

Tend = nTracer_Dclock(FALSE, FALSE)

TimeChk%SubGrpRTniter   = TimeChk%SubGrpRTniter   + itersum
TimeChk%SubGrpTime      = TimeChk%SubGrpTime      + (Tend - Tbeg)
TimeChk%NetRTSubGrpTime = TimeChk%NetRTSubGrpTime + rtnet
! ----------------------------------------------------

END SUBROUTINE SubGrpFsp_MLG_GM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE CopyFlux(phis, phisd, nFsr)

IMPLICIT NONE

REAL, POINTER, DIMENSION(:) :: phis, phisd
INTEGER :: nFsr
INTEGER :: ifsr
! ----------------------------------------------------

#ifdef __INTEL_MKL
CALL dcopy(nFsr, phis, 1, phisd, 1)
#else
!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED)
DO ifsr = 1, nFsr
  phisd(ifsr) = phis(ifsr)
ENDDO
!$OMP END DO
!$OMP END PARALLEL
#endif
! ----------------------------------------------------

END SUBROUTINE CopyFlux
! ------------------------------------------------------------------------------------------------------------

END MODULE SubGrpFspMLG_GM