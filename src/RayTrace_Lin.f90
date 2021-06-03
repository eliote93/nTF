!#define newslope
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceLS(RayInfo, CoreInfo, phis, PhiAngIn, xst, src, LinSrc, Slope1g, jout, iz, ljout)

USE PARAM
USE TYPEDEF,  ONLY :  RayInfo_Type,      coreinfo_type,                                       &
                      Pin_Type,          Asy_Type,        AsyInfo_Type,     PinInfo_Type,     &
                      Cell_Type,                                                              &
                      AziAngleInfo_Type, PolarAngle_Type, ModRayInfo_type,  AsyRayInfo_type,  &
                      CoreRayInfo_Type,  RotRayInfo_Type, CellRayInfo_type
USE Moc_Mod, ONLY :   nMaxRaySeg,        nMaxCellRay,     nMaxAsyRay,       nMaxCoreRay,      &
                      ExpA,              ExpB,                                                &
                      ApproxExp
USE BasicOperation, ONLY : CP_CA, CP_VA                    
USE ALLOCS
IMPLICIT NONE
!Input Arguments
TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: phis(:), PhiAngIn(:, :), xst(:), src(:), LinSrc(:, :), jout(:, :, :)
REAL :: Slope1g(:, :)
!REAL :: phis(:), PhiAngIn(:, :), xst(:), src(:), jout(:, :, :)
INTEGER :: iz
LOGICAL :: ljout
!Pointing variablez
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Asy_Type), POINTER :: Asy(:)
!TYPE(AsyInfo_Type), POINTER :: AsyInfo
TYPE(PinInfo_Type), POINTER :: PinInfo(:)
TYPE(Cell_Type), POINTER :: Cell(:)
TYPE(AziAngleInfo_Type), POINTER :: AziAng(:)
TYPE(PolarAngle_Type), POINTER :: PolarAng(:)
!TYPE(ModRayInfo_type), POINTER :: ModRay
TYPE(AsyRayInfo_type), POINTER :: AsyRay(:)
TYPE(CoreRayInfo_Type), POINTER :: CoreRay(:)
TYPE(RotRayInfo_Type), POINTER :: RotRay(:)
TYPE(CellRayInfo_Type), POINTER :: CellRay

REAL, POINTER :: LenSeg(:)
INTEGER, POINTER :: LocalFsrIdx(:)

!Local Variable
LOGICAL, SAVE :: lfirst
INTEGER :: mp(2)
INTEGER :: iazi, ipol, irotray, iCoreRay, iasyray, iceray, irayseg, irot, itype, idir     !Ray related index
INTEGER :: nRotRay, nCoreRay, nAsyRay, nPinRay, nRaySeg, FsrIdxSt
INTEGER :: nPolarAng, nAziAng
INTEGER :: ipin, icel, iasy, ireg, isurf                                                  !Geometries index
INTEGER :: irsegidx, icellrayidx, PhiAnginSvIdx, PhiAngOutSvIdx
INTEGER :: nFsr, nxy
INTEGER :: i, j, k, l, m, jbeg, jend, jinc, ir, ir1
!Ray Data Variables
!#define stackvar

!INTEGER, POINTER, SAVE :: CellRayIdxSt(:, :, :),  PinIdx(:, :)             !Starting Adress of CellRay, Fsr Index, Pin Index, Exp function approx. index
!INTEGER, POINTER, SAVE :: SurfIdx(:, :, :), nTotRaySeg(:), nTotCellRay(:)   !Surface Index list, # tot ray segments in rot ray, # total Cell Ray index
INTEGER :: CellRayIdxSt(nMaxCellRay, nMaxCoreRay, 2)
INTEGER :: PinIdx(nMaxCellRay, nMaxCoreRay)
INTEGER :: SurfIdx(nMaxCellRay, nMaxCoreRay, 2)
INTEGER :: nTotRaySeg(nMaxCoreRay)
INTEGER :: nTotCellRay(nMaxCoreRay)

!REAL :: PhiAngOut(0:nMaxRaySeg + 1)
!
#ifndef  stackvar
REAL, POINTER :: PhiAngOut(:)
INTEGER, POINTER :: FsrIdx(:, :),  ExpAppIdx(:, :)
REAL, POINTER :: OptLenList(:, :), ExpApp(:, :)  !Optical Length, Approximated Exponential function list
REAL, POINTER :: SegLen(:, :), ExpApp1(:, :), ExpApp2(:, :)     !
#else
INTEGER :: FsrIdx(nMaxRaySeg, nMaxCoreRay)
INTEGER :: ExpAppIdx(nMaxRaySeg, nMaxCoreRay)
REAL :: OptLenList(nMaxRaySeg, nMaxCoreRay), ExpApp(nMaxRaySeg, nMaxCoreRay)
#endif
REAL, SAVE, POINTER :: asrc(:)


REAL :: wtang(100, 100), wttemp
REAL :: tau, phiobd, phid, wt
REAL :: rsinv, sinv2                   !Inverse of sine, sine square

REAL :: wtx, wty
REAL :: SlopeWt(2, 100, 100) 
!REAL, POINTER, SAVE :: Slope1g(:, :)
!REAL, POINTER, SAVE :: !SlopeX1g
!REAL, POINTER, SAVE :: SegLen(:, :)


DATA mp /2, 1/
DATA lfirst /TRUE/



!Ray Info Pointing
AziAng => RayInfo%AziAngle;
PolarAng => RayInfo%PolarAngle;
AsyRay => RayInfo%AsyRay
CoreRay => RayInfo%CoreRay
RotRay => RayInfo%RotRay

!Geometry Info Pointing
Asy => CoreInfo%Asy
Pin => CoreInfo%Pin
PinInfo => CoreInfo%Pininfo
Cell => CoreInfo%CellInfo

nAziAng = RayInfo%nAziAngle; nPolarAng = RayInfo%nPolarAngle
nFsr = CoreInfo%nCoreFsr; nxy = CoreInfo%nxy

ALLOCATE(PhiAngOut(0:nMaxRaySeg + 1))


IF(lFirst) THEN
  !nRotRay = RayInfo%nRotRay
  CALL Dmalloc(asrc, nFsr)

#ifndef  stackvar
  CALL Dmalloc(FsrIdx, nMaxRaySeg, nMaxCoreRay)
  CALL Dmalloc(ExpAppIdx, nMaxRaySeg, nMaxCoreRay)

  CALL Dmalloc(OptLenList, nMaxRaySeg, nMaxCoreRay)
  CALL Dmalloc(ExpApp, nMaxRaySeg, nMaxCoreRay)
  CALL Dmalloc(ExpApp1, nMaxRaySeg, nMaxCoreRay)

  CALL Dmalloc0(PhiAngOut, 0, nMaxRaySeg + 1)
!  CALL Dmalloc(PinIdx, nMaxCellRay, nMaxCoreRay);  CALL Dmalloc(SurfIdx, nMaxCellRay, nMaxCoreRay, 2)
!  CALL Dmalloc(CellRayIdxSt, nMaxCellRay, nMaxCoreRay, 2)
!  CALL Dmalloc(nTotRaySeg, nMaxCoreRay);  CALL Dmalloc(nTotCellRay, nMaxCoreRay)
#endif  
  !CALL ApproxExp(PolarAng, nPolarAng, ExpA, ExpB)
  CALL ApproxExp(PolarAng, nPolarAng)
#ifdef newslope
  CALL Dmalloc(ExpApp2, nMaxRaySeg, nMaxCoreRay)
#endif
  lFirst = FALSE
ENDIF

!Calculation Weighting 
DO ipol = 1, nPolarAng
  wttemp = PolarAng(ipol)%weight  * PolarAng(ipol)%sinv
  DO iazi = 1, nAziAng
    wtang(ipol, iazi) = wttemp  * AziAng(iazi)%weight * AziAng(iazi)%del
  ENDDO
ENDDO

DO ipol = 1, nPolarAng
#ifndef newslope
  wttemp = PolarAng(ipol)%weight * PolarAng(ipol)%sinv !* PolarAng(ipol)%sinv
#else
  wttemp = PolarAng(ipol)%weight * PolarAng(ipol)%sinv * PolarAng(ipol)%sinv
#endif  
  DO iazi = 1, nAziAng
    SlopeWt(1, ipol, iazi) = wttemp  * AziAng(iazi)%weight* AziAng(iazi)%cosv * AziAng(iazi)%del
    SlopeWt(2, ipol, iazi) = wttemp  * AziAng(iazi)%weight* AziAng(iazi)%sinv * AziAng(iazi)%del 
  ENDDO
ENDDO

!Initialize of scalar flux
CALL CP_CA(phis, ZERO, nFsr)
IF(ljout) CALL CP_CA(Jout, ZERO, 2, 4, CoreInfo%nxy)

CALL CP_CA(Slope1g, ZERO, 2, nFsr)
!CALL CP_CA(src, 0.5_8, nfsr)
wt =0
wttemp =1000
DO i = 1, nFSR
  wt = max(wt, src(i))
  wttemp = min(wttemp, src(i))
ENDDO
nRotRay = RayInfo%nRotRay
DO i = 1, nRotRay    !Rotational Ray Sweep 
  !Make CoreRay 1-D Array
  iRotRay = i
  nCoreRay = RotRay(irotRay)%nRay
  DO j = 1, nCoreRay    !Core Ray Sweep
    irsegidx = 0;   icellrayidx = 0
    iCoreRay = RotRay(iRotRay)%RayIdx(j)
    nAsyRay = CoreRay(iCoreRay)%nRay
    DO k = 1, nAsyRay    !Assembly Ray Sweep
      iasyray = CoreRay(iCoreRay)%AsyRayIdx(k)
      iasy = CoreRay(iCoreRay)%AsyIdx(k)
      IF(iasy .EQ. 0)  CYCLE   !Skip Dummy Assembly
      nPinRay = AsyRay(iAsyRay)%nCellRay
      itype = Asy(iasy)%PartialAsyFlag
      DO l = 1, nPinRay   !Pin Ray Sweep
        ipin = AsyRay(iAsyRay)%PinIdx(l)      !Local Pin Idx(within Assembly)
        iceray = AsyRay(iAsyRay)%PinRayIdx(l) !Cell Ray Index
        ipin = Asy(iAsy)%GlobalPinIdx(ipin)   !Global Pin Index
        icel = Pin(ipin)%Cell(iz)             !Cell Type
        FsrIdxSt = Pin(ipin)%FsrIdxSt
        CellRay => Cell(icel)%CellRay(iceray) !Pointing Cell Ray
        icellrayidx = icellrayidx + 1
        PinIdx(icellrayidx, j) = ipin
        CellRayIdxSt(icellrayidx, j, 2) = irsegidx + 1
        nRaySeg = CellRay%nSeg
        LocalFsrIdx => CellRay%LocalFsrIdx
        LenSeg => CellRay%LenSeg
        DO iRaySeg = 1, nRaySeg
          ireg = FsrIdxSt + CellRay%LocalFsrIdx(iRaySeg) - 1 
          tau = - LenSeg(iRaySeg) * xst(ireg)   !
          !Global Fsr Index
          !tau = -1000._8 * CellRay%LenSeg(iRaySeg) * xst(ireg)   !
          tau = - CellRay%LenSeg(iRaySeg) * xst(ireg)   !
          irsegidx = irsegidx + 1
          FsrIdx(irsegidx, j) = ireg
          OptLenList(irsegidx, j) = tau
          ExpAppIdx(irsegidx, j) = max(INT(tau), -40000)
          !OptLenList(irsegidx, j) = OptLenList(irsegidx, j)/1000._8
          !SegLen(irsegidx, j) = LenSeg(iRaySeg) * 0.001_8
          continue
        ENDDO   !End of Ray Segments Sweep, irayseg
        CellRayIdxSt(icellrayidx, j, 1) = irsegidx
        SurfIdx(icellRayIdx, j, 1) = AsyRay(iAsyRay)%PinRaySurf(2, l) !OutSurface
        SurfIdx(icellRayIdx, j, 2) = AsyRay(iAsyRay)%PinRaySurf(1, l) !Insurface
      ENDDO    !End of Pin Ray Seep, l
    ENDDO    !End of Asy Ray Swwep, k
    nTotRaySeg(j) = irsegidx
    nTotCellRay(j) = icellRayIdx
  ENDDO    !End of Core Ray Sweep, j
  
  DO ipol = 1, nPolarAng
    rsinv = 1._8 / PolarAng(ipol)%sinv
    sinv2 = PolarAng(ipol)%sinv 
    !Approximate 1-exp
    DO j = 1, nCoreRay 
      DO l = 1, nTotRaySeg(j)
        tau = -0.001_8 * OptLenList(l, j) * rsinv
        ExpApp(l, j) = ExpA(ExpAppIdx(l, j), ipol)*optlenlist(l, j) + ExpB(ExpAppIdx(l, j), ipol)
        ExpApp1(l, j) = 2._8 * (tau - ExpApp(l, j)) - tau * ExpApp(l, j)
        ExpApp1(l, j) = ExpApp1(l, j) * sinv2
       ! ExpApp2(l, j) = ExpApp(l, j) -tau + 0.5_8 * tau * tau
       ! ExpApp2(l, j) = ExpApp2(l, j) * sinv2 ! / rsinv
      ENDDO 
    ENDDO
    
    DO irot = 1, 2
      PhiAnginSvIdx = RayInfo%PhiAngInSvIdx(iRotRay ,irot)
      PhiAngOutSvIdx = RayInfo%PhiangOutSvIdx(iRotRay ,irot)
      phiobd = PhiAngIn(ipol,PhiAnginSvIdx)
      
      jinc = 1; jbeg = 1; jend = nCoreRay
      IF(irot .eq. 2) THEN !Backward Sweep
        jinc = -1; jbeg = nCoreRay; jend = 1
      ENDIF
      DO j = jbeg, jend, jinc
        idir = RotRay(i)%DIR(j); iazi = CoreRay(RotRay(irotray)%RayIdx(j))%iang
        IF(irot .eq. 2) idir = mp(idir)  !Reverse the sweep direction
        wt = wtang(ipol, iazi)
        nRaySeg = nTotRaySeg(j)
        IF(idir .eq. 1) THEN  !Forward Sweep
          PhiAngOut(0) = phiobd
          DO ir = 1, nRaySeg
            ireg = FsrIdx(ir, j)            
            phid = (PhiAngOut(ir - 1) - src(ireg)) * ExpApp(ir, j)
            phid = phid - LinSrc(ireg, iazi) * ExpApp1(ir, j)
            PhiAngOut(ir) = PhiAngOut(ir - 1) - phid
            phis(ireg) = phis(ireg) + wt*phid  ! w_pol * Sin_Pol * w_azi * del * (phi_in - phi_out)
          ENDDO
          wtx = SlopeWt(1, ipol, iazi)
          wty = SlopeWt(2, ipol, iazi)
          DO ir = 1, nRaySeg
            ireg = FsrIdx(ir, j)
#ifndef newslope            
            Slope1g(1, ireg) = Slope1g(1, ireg) + (PhiAngOut(ir) - PhiAngOut(ir - 1)) * wtx   
            Slope1g(2, ireg) = Slope1g(2, ireg) + (PhiAngOut(ir) - PhiAngOut(ir - 1)) * wty
#else
            phid =(PhiAngOut(ir - 1) - PhiAngOut(ir)) + LinSrc(ireg, iazi) * ExpApp2(ir, j)
            Slope1g(1, ireg) = Slope1g(1, ireg) + phid * wtx
            Slope1g(2, ireg) = Slope1g(2, ireg) + phid * wty
#endif            
          ENDDO          
          
          phiobd = PhiAngOut(nRaySeg)
          !Surface 
          IF(ljout) THEN
            DO ir = 1, nTotCellRay(j)
              icel = PinIdx(ir, j); isurf = SurfIdx(ir, j, 1)
              Jout(2, isurf, icel) = Jout(2, isurf, icel) + wt * PhiAngOut(CellRayIdxSt(ir, j, 1))
              isurf = SurfIdx(ir, j, 2)
              Jout(1, isurf, icel) = Jout(1, isurf, icel) + wt * PhiAngOut(CellRayIdxSt(ir, j, 2)-1)
            ENDDO
          ENDIF
        ELSE
          PhiAngOut(nRaySeg+1) = phiobd
          ir = nRaySeg + 1
          DO ir1 = 1, nRaySeg
            ir = ir - 1
            ireg = FsrIdx(ir, j)
            phid = (PhiAngOut(ir + 1) - src(ireg)) * ExpApp(ir, j)
            phid = phid + LinSrc(ireg, iazi) * ExpApp1(ir, j)
            PhiAngOut(ir) = PhiAngOut(ir + 1) - phid
            phis(ireg) = phis(ireg) + wt * phid
          ENDDO
          phiobd = PhiAngOut(1)
          
          wtx = SlopeWt(1, ipol, iazi)
          wty = SlopeWt(2, ipol, iazi)
          ir = nRaySeg + 1
          DO ir1 = 1, nRaySeg
            ir = ir - 1          
            ireg = FsrIdx(ir, j)
#ifndef newslope 
            Slope1g(1, ireg) = Slope1g(1, ireg) - (PhiAngOut(ir) - PhiAngOut(ir + 1)) * wtx   
            Slope1g(2, ireg) = Slope1g(2, ireg) - (PhiAngOut(ir) - PhiAngOut(ir + 1)) * wty
#else
            phid = (PhiAngOut(ir + 1) - PhiAngOut(ir)) - LinSrc(ireg, iazi) * ExpApp2(ir, j)
            Slope1g(1, ireg) = Slope1g(1, ireg) - phid * wtx
            Slope1g(2, ireg) = Slope1g(2, ireg) - phid * wty
#endif
          ENDDO            
          IF(lJout) THEN
            DO ir = 1, nTotCellRay(j)
              icel = PinIdx(ir, j); isurf = SurfIdx(ir, j, 2)
              Jout(2, isurf, icel) = Jout(2, isurf, icel) + wt * PhiAngOut(CellRayIdxSt(ir, j, 2))
              isurf = SurfIdx(ir, j, 1)
              Jout(1, isurf, icel) = Jout(1, isurf, icel) + wt * PhiAngOut(CellRayIdxSt(ir, j, 1) + 1)
            ENDDO
          ENDIF
        ENDIF
      ENDDO !End of CoreRay Sweep
      PhiAngIn(ipol,PhiAngOutSvIdx) = phiobd
    ENDDO !Backward and forward Sweep
  ENDDO !Polar Angle Sweep
ENDDO    !End of Rotation Ray Sweep, i

!Add Q contribution 
DO l = 1, nxy
  FsrIdxSt = Pin(l)%FsrIdxSt; icel = Pin(l)%Cell(iz);
  DO j = 1, Cell(icel)%nFsr
    ireg = FsrIdxSt + j - 1
    phis(ireg) = phis(ireg)/xst(ireg)/Cell(icel)%vol(j) + src(ireg)
    !phis(ireg) =  src(ireg)
  ENDDO  
ENDDO

!
DO l = 1, nxy
  FsrIdxSt = Pin(l)%FsrIdxSt; icel = Pin(l)%Cell(iz);
  DO j = 1, Cell(icel)%nFsr
    ireg = FsrIdxSt + j - 1
#ifndef newslope           
    Slope1g(1, ireg) = Slope1g(1, ireg) / (Cell(icel)%vol(j))
    Slope1g(2, ireg) = Slope1g(2, ireg) / (Cell(icel)%vol(j))
#else
    Slope1g(1, ireg) = -Slope1g(1, ireg) / (Cell(icel)%vol(j) * xst(ireg))
    Slope1g(2, ireg) = -Slope1g(2, ireg) / (Cell(icel)%vol(j) * xst(ireg))
#endif    
  ENDDO  
ENDDO


!Free the pointing variables
NULLIFY(AziAng); NULLIFY(PolarAng)
NULLIFY(AsyRay); NULLIFY(CoreRay)
NULLIFY(RotRay); NULLIFY(Asy)
NULLIFY(Pin); NULLIFY(PinInfo)
NULLIFY(Cell); NULLIFY(CellRay)
NULLIFY(LenSeg); NULLIFY(LocalFsrIdx)
DEALLOCATE(PhiAngOut)

END SUBROUTINE RayTraceLS
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceLS_CASMO(RayInfo, CoreInfo, phisnm, phisSlope, PhiAngInnm, srcnm, srcSlope, xstnm, joutnm, iz, gb, ge, ljout)

USE PARAM
USE TYPEDEF,    ONLY :  RayInfo_Type,       Coreinfo_type,      Pin_Type,           Asy_Type,               &
                        AsyInfo_Type,       PinInfo_Type,       Cell_Type,          AziAngleInfo_Type,      &
                        PolarAngle_Type,    ModRayInfo_type,    AsyRayInfo_type,    CoreRayInfo_Type,       &
                        RotRayInfo_Type,    CellRayInfo_type,   DcmpAsyRayInfo_Type
USE Moc_Mod,    ONLY :  nMaxRaySeg,         nMaxCellRay,        nMaxAsyRay,         nMaxCoreRay,            &
                        ExpA,               ExpB,               ApproxExp,          TrackingDat,            &
                        DcmpPhiAngIn,       DcmpPhiAngOut,      TrackRotRayLSDcmp_CASMO
USE geom,       ONLY :  ng
USE PE_Mod,     ONLY :  PE
USE ALLOCS
USE OMP_LIB
IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: phisnm(:, :), PhiAngInnm(:, :, :)
REAL, POINTER :: srcnm(:, :), xstnm(:, :), joutnm(:, :, :, :)
REAL, POINTER :: phisSlope(:, :, :, :), srcSlope(:, :, :, :)
INTEGER :: iz, gb, ge
LOGICAL :: ljout

LOGICAL, SAVE :: lfirst
DATA lfirst /.TRUE./

TYPE(Cell_Type), POINTER :: Cell(:)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(DcmpAsyRayInfo_Type), POINTER :: DcmpAsyRay(:, :)
INTEGER, POINTER :: DcmpAsyRayCount(:)
REAL, POINTER :: FsrMxx(:), FsrMyy(:), FsrMxy(:)

INTEGER :: nAziAng, nPolarAng, nPhiAngSv
INTEGER :: nRotray
INTEGER :: nFsr, nAsy, nxy
INTEGER :: nThread

INTEGER :: iRotRay, iAsyRay, iAsy
INTEGER :: tid
INTEGER :: FsrIdxSt, icel, ireg
INTEGER :: i, j, l, g
REAL :: detinv

REAL, POINTER :: phimx(:, :, :), phimy(:, :, :)

Cell => CoreInfo%CellInfo; Pin => CoreInfo%Pin
FsrMxx => CoreInfo%CoreFsrMxx(:, iz)
FsrMyy => CoreInfo%CoreFsrMyy(:, iz)
FsrMxy => CoreInfo%CoreFsrMxy(:, iz)

nAsy = CoreInfo%nxya
nAziAng = RayInfo%nAziAngle; nPolarAng = RayInfo%nPolarAngle
nPhiAngSv = RayInfo%nPhiAngSv
nRotRay = RayInfo%nRotRay
nFsr = CoreInfo%nCoreFsr; nxy = CoreInfo%nxy
nThread = PE%nThread

DcmpAsyRay => RayInfo%DcmpAsyRay
DcmpAsyRayCount => RayInfo%DcmpAsyRayCount

ALLOCATE(phimx(2, gb : ge, nFsr), phimy(2, gb : ge, nFsr))
phimx(:, :, :) = zero; phimy(:, :, :) = zero

IF(lfirst) THEN
  lFirst = FALSE
  CALL ApproxExp(RayInfo%PolarAngle, nPolarAng)
  DO tid = 1, nThread
    IF(TrackingDat(tid)%lAllocLinSrc) CYCLE
    CALL Dmalloc(TrackingDat(tid)%FsrIdx, nMaxRaySeg, nMaxCoreRay)
    CALL Dmalloc(TrackingDat(tid)%ExpAppIdxnm, ng, nMaxRaySeg, nMaxCoreRay)
    CALL Dmalloc(TrackingDat(tid)%OptLenListnm, ng, nMaxRaySeg, nMaxCoreRay)
    CALL Dmalloc(TrackingDat(tid)%cmOptLen, nPolarAng, ng, nMaxRaySeg, nMaxCoreRay)
    CALL Dmalloc(TrackingDat(tid)%cmOptLenInv, nPolarAng, ng, nMaxRaySeg, nMaxCoreRay)
    CALL Dmalloc(TrackingDat(tid)%PhiAngOutnm, nPolarAng, ng, nMaxRaySeg + 2)
    CALL Dmalloc(TrackingDat(tid)%Joutnm, 3, ng, 4, nxy)
    CALL Dmalloc(TrackingDat(tid)%q0, ng, nMaxRaySeg, nMaxCoreRay)
    CALL Dmalloc(TrackingDat(tid)%q1, nPolarAng, ng, nMaxRaySeg, nMaxCoreRay)
    CALL Dmalloc(TrackingDat(tid)%phisnm, ng, nFsr)
    CALL Dmalloc(TrackingDat(tid)%phimx, 2, ng, nFsr)
    CALL Dmalloc(TrackingDat(tid)%phimy, 2, ng, nFsr)
    CALL Dmalloc(TrackingDat(tid)%E1, nPolarAng, ng, nMaxRaySeg, nMaxCoreRay)
    CALL Dmalloc(TrackingDat(tid)%E3, nPolarAng, ng, nMaxRaySeg, nMaxCoreRay)
    CALL Dmalloc(TrackingDat(tid)%R1, nPolarAng, ng, nMaxRaySeg, nMaxCoreRay)
    CALL Dmalloc(TrackingDat(tid)%R3, nPolarAng, ng, nMaxRaySeg, nMaxCoreRay)
    CALL Dmalloc(TrackingDat(tid)%x0, 2, nMaxRaySeg, nMaxCoreRay)
    CALL Dmalloc(TrackingDat(tid)%y0, 2, nMaxRaySeg, nMaxCoreRay)
    TrackingDat(tid)%ExpA => ExpA; TrackingDat(tid)%ExpB => ExpB
    TrackingDat(tid)%PhiAngInnm => PhiAngInnm
    TrackingDat(tid)%srcnm => srcnm; TrackingDat(tid)%srcSlope => srcSlope(:, :, :, iz)
    TrackingDat(tid)%xstnm => xstnm
    TrackingDat(tid)%DcmpPhiAngIn => DcmpPhiAngIn
    TrackingDat(tid)%DcmpPhiAngOut => DcmpPhiAngOut
    TrackingDat(tid)%lAllocLinSrc = .TRUE.
  ENDDO
ENDIF

CALL omp_set_dynamic(.FALSE.)
CALL omp_set_num_threads(nThread)

phisnm(gb : ge, :) = zero
IF(ljout) joutnm(:, gb : ge, :, :) = zero !--- BYS edit / 150612 Surface flux

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(irotray, tid, i, j)
tid = omp_get_thread_num() + 1
!$OMP MASTER
DO j = 1, nThread
  TrackingDat(j)%phisnm(:, :) = zero
  TrackingDat(j)%phimx(:, :, :) = zero
  TrackingDat(j)%phimy(:, :, :) = zero
  TrackingDat(j)%joutnm(:, :, :, :) = zero
ENDDO
!$OMP END MASTER
!$OMP BARRIER

!$OMP DO SCHEDULE(GUIDED)
DO i = 1, nRotRay
  CALL TrackRotRayLS_CASMO(RayInfo, CoreInfo, TrackingDat(tid), ljout, i, iz, gb, ge)
ENDDO
!$OMP END DO
!$OMP END PARALLEL

DO j = 1, nThread
  DO i = 1, nFsr
    DO g = gb, ge
      phisnm(g, i) = phisnm(g, i) + TrackingDat(j)%phisnm(g, i)
      phimx(:, g, i) = phimx(:, g, i) + TrackingDat(j)%phimx(:, g, i)
      phimy(:, g, i) = phimy(:, g, i) + TrackingDat(j)%phimy(:, g, i)
    ENDDO
  ENDDO
ENDDO

DO i = 1, nFsr
  DO g = gb, ge
    phimx(1, g, i) = phimx(1, g, i) / xstnm(g, i)
    phimy(1, g, i) = phimy(1, g, i) / xstnm(g, i)
  ENDDO
ENDDO

IF(ljout) THEN
  DO j = 1, nThread
    joutnm(:, gb : ge, :, :) = joutnm(:, gb : ge, :, :) + TrackingDat(j)%joutnm(:, gb : ge, :, :)
  ENDDO
ENDIF

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l, j, FsrIdxSt, icel, ireg)
!$OMP DO
DO l = 1, nxy
  FsrIdxSt = Pin(l)%FsrIdxSt; icel = Pin(l)%Cell(iz);
  DO j = 1, Cell(icel)%nFsr
    ireg = FsrIdxSt + j - 1
    DO g = gb, ge
      phisnm(g, ireg) = phisnm(g, ireg) / (xstnm(g, ireg) * Cell(icel)%vol(j)) + srcnm(g, ireg)
      phimx(:, g, ireg) = phimx(:, g, ireg) / (xstnm(g, ireg) * Cell(icel)%vol(j))
      phimy(:, g, ireg) = phimy(:, g, ireg) / (xstnm(g, ireg) * Cell(icel)%vol(j))
    ENDDO
  ENDDO  
ENDDO
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(detinv, i, g)
DO i = 1, nFsr
  detinv = 1.0 / (FsrMxx(i) * FsrMyy(i) - FsrMxy(i) * FsrMxy(i))
  DO g = gb, ge
    phisSlope(1, g, i, iz) = detinv * (FsrMyy(i) * SUM(phimx(:, g, i)) - FsrMxy(i) * SUM(phimy(:, g, i)))
    phisSlope(2, g, i, iz) = detinv * (FsrMxx(i) * SUM(phimy(:, g, i)) - FsrMxy(i) * SUM(phimx(:, g, i)))
  ENDDO
ENDDO
!$OMP END PARALLEL DO

DEALLOCATE(phimx, phimy)

NULLIFY(Cell); NULLIFY(Pin)
NULLIFY(FsrMxx); NULLIFY(FsrMyy); NULLIFY(FsrMxy)

END SUBROUTINE RayTraceLS_CASMO
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE TrackRotRayLS_CASMO(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, gb, ge)
USE PARAM
USE TYPEDEF,    ONLY :  RayInfo_Type,       Coreinfo_type,      Pin_Type,           Asy_Type,               &
                        AsyInfo_Type,       PinInfo_Type,       Cell_Type,          AziAngleInfo_Type,      &
                        PolarAngle_Type,    AsyRayInfo_type,    CoreRayInfo_Type,                           &
                        RotRayInfo_Type,    CellRayInfo_type,   TrackingDat_Type
USE Moc_Mod,    ONLY :  nMaxRaySeg,         nMaxCellRay,        nMaxAsyRay,         nMaxCoreRay
USE BasicOperation, ONLY : CP_CA, CP_VA
IMPLICIT NONE

TYPE(RayInfo_Type), INTENT(INOUT) :: RayInfo
TYPE(CoreInfo_Type), INTENT(INOUT) :: CoreInfo
TYPE(TrackingDat_Type), INTENT(INOUT) :: TrackingDat
LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: irotray, iz, gb, ge

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Asy_Type), POINTER :: Asy(:)
TYPE(PinInfo_Type), POINTER :: PinInfo(:)
TYPE(Cell_Type), POINTER :: Cell(:)
TYPE(AziAngleInfo_Type), POINTER :: AziAng(:)
TYPE(PolarAngle_Type), POINTER :: PolarAng(:)
TYPE(AsyRayInfo_type), POINTER :: AsyRay(:)
TYPE(CoreRayInfo_Type), POINTER :: CoreRay(:)
TYPE(RotRayInfo_Type), POINTER :: RotRay(:)
TYPE(CellRayInfo_Type), POINTER :: CellRay

REAL, POINTER :: LenSeg(:)
INTEGER, POINTER :: LocalFsrIdx(:)
INTEGER, POINTER :: FsrIdx(:, :),  ExpAppIdx(:, :, :)
REAL, POINTER :: OptLenList(:, :, :)
REAL, POINTER :: phisC(:, :), phimx(:, :, :), phimy(:, :, :)
REAL, POINTER :: srcC(:, :), srcSlope(:, :, :), xst(:, :), jout(:, :, :, :)
REAL, POINTER :: PhiAngOut(:, :, :), PhiAngIn(:, :, :)
REAL, POINTER :: ExpA(:, :), ExpB(:, :)
REAL, POINTER :: E1(:, :, :, :), E3(:, :, :, :), R1(:, :, :, :), R3(:, :, :, :)
REAL, POINTER :: cmOptLen(:, :, :, :), cmOptLenInv(:, :, :, :)
REAL, POINTER :: q0(:, :, :), q1(:, :, :, :)
REAL, POINTER :: x0(:, :, :), y0(:, :, :)
REAL, POINTER :: FsrCentroid(:, :)

INTEGER :: mp(2)
INTEGER :: iazi, ipol, ig, iCoreRay, iasyray, iceray, irayseg, irot, itype, idir
INTEGER :: nRotRay, nCoreRay, nAsyRay, nPinRay, nRaySeg, FsrIdxSt
INTEGER :: nPolarAng, nAziAng, nPhiAngSv
INTEGER :: ipin, icel, iasy, ireg, isurf
INTEGER :: irsegidx, icellrayidx, PhiAnginSvIdx, PhiAngOutSvIdx
INTEGER :: nFsr, nxy
INTEGER :: i, j, k, l, m, jbeg, jend, jinc, ir, ir1
INTEGER :: ibcel

INTEGER :: CellRayIdxSt(nMaxCellRay, nMaxCoreRay, 2)
INTEGER :: PinIdx(nMaxCellRay, nMaxCoreRay)
INTEGER :: SurfIdx(nMaxCellRay, nMaxCoreRay, 2)
INTEGER :: nTotRaySeg(nMaxCoreRay)
INTEGER :: nTotCellRay(nMaxCoreRay)

REAL :: wtang(10, 100), wttemp, wtang2(10, 100, 4)
REAL :: tau, phiobd(10, gb : ge), phid, phim, wt(10), wt2(10, 4)
REAL :: ax(10), ay(10), polsininv(10)
REAL :: segCenter(2), segSlope(2, 10)

DATA mp /2, 1/

!Ray Info Pointing
AziAng => RayInfo%AziAngle;
PolarAng => RayInfo%PolarAngle;
AsyRay => RayInfo%AsyRay
CoreRay => RayInfo%CoreRay
RotRay => RayInfo%RotRay

!Geometry Info Pointing
Asy => CoreInfo%Asy
Pin => CoreInfo%Pin
PinInfo => CoreInfo%Pininfo
Cell => CoreInfo%CellInfo
FsrCentroid => CoreInfo%CoreFsrCentroid(:, :, iz)

!Tracking Dat Pointing
FsrIdx => TrackingDat%FsrIdx
ExpAppIdx => TrackingDat%ExpAppIdxnm
OptLenList => TrackingDat%OptLenListnm
cmOptLen => TrackingDat%cmOptLen; cmOptLenInv => TrackingDat%cmOptLenInv
phisC => TrackingDat%phisnm; phimx => TrackingDat%phimx; phimy => TrackingDat%phimy
srcC => TrackingDat%srcnm; srcSlope => TrackingDat%srcSlope
xst => TrackingDat%xstnm; jout => TrackingDat%joutnm
PhiAngOut => TrackingDat%PhiAngOutnm
PhiAngIn => TrackingDat%phiAngInnm
ExpA => TrackingDat%ExpA
ExpB => TrackingDat%ExpB
E1 => TrackingDat%E1; E3 => TrackingDat%E3
R1 => TrackingDat%R1; R3 => TrackingDat%R3
q0 => TrackingDat%q0; q1 => TrackingDat%q1
x0 => TrackingDat%x0; y0 => TrackingDat%y0

nAziAng = RayInfo%nAziAngle; nPolarAng = RayInfo%nPolarAngle
nPhiAngSv = RayInfo%nPhiAngSv
nFsr = CoreInfo%nCoreFsr; nxy = CoreInfo%nxy

DO ipol = 1, nPolarAng
  wttemp = PolarAng(ipol)%weight * PolarAng(ipol)%sinv
  polsininv(ipol) = 1._8 / PolarAng(ipol)%sinv
  DO iazi = 1, nAziAng
    wtang(ipol, iazi) = wttemp * AziAng(iazi)%weight * AziAng(iazi)%del
  ENDDO
ENDDO
IF(lJout)THEN
  DO ipol = 1, nPolarAng
    DO iazi = 1, nAziAng
      wtang2(ipol, iazi, 1) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%sinv)
      wtang2(ipol, iazi, 3) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%sinv)
      wtang2(ipol, iazi, 2) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%cosv)
      wtang2(ipol, iazi, 4) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%cosv)
    ENDDO
  ENDDO
ENDIF

i = iRotRay
nCoreRay = RotRay(irotRay)%nRay
DO j = 1, nCoreRay
  irsegidx = 0;   icellrayidx = 0
  iCoreRay = RotRay(iRotRay)%RayIdx(j)
  nAsyRay = CoreRay(iCoreRay)%nRay
  iazi = CoreRay(iCoreRay)%iang
  DO ipol = 1, nPolarAng
    segSlope(:, ipol) = (/ AziAng(iazi)%cosv * PolarAng(ipol)%sinv, AziAng(iazi)%sinv * PolarAng(ipol)%sinv /)
  ENDDO
  DO k = 1, nAsyRay
    iasyray = CoreRay(iCoreRay)%AsyRayIdx(k)
    iasy = CoreRay(iCoreRay)%AsyIdx(k)
    IF(iasy .EQ. 0)  CYCLE
    nPinRay = AsyRay(iAsyRay)%nCellRay
    itype = Asy(iasy)%PartialAsyFlag
    DO l = 1, nPinRay
      ipin = AsyRay(iAsyRay)%PinIdx(l)
      iceray = AsyRay(iAsyRay)%PinRayIdx(l)
      ipin = Asy(iAsy)%GlobalPinIdx(ipin)
      icel = Pin(ipin)%Cell(iz)
      FsrIdxSt = Pin(ipin)%FsrIdxSt
      ibcel = Cell(icel)%basecellstr
      CellRay => Cell(ibcel)%CellRay(iceray)
      icellrayidx = icellrayidx + 1
      PinIdx(icellrayidx, j) = ipin
      CellRayIdxSt(icellrayidx, j, 2) = irsegidx + 1
      nRaySeg = CellRay%nSeg
      LocalFsrIdx => CellRay%LocalFsrIdx
      LenSeg => CellRay%LenSeg
      DO iRaySeg = 1, nRaySeg
        ireg = FsrIdxSt + CellRay%LocalFsrIdx(iRaySeg) - 1
        irsegidx = irsegidx + 1
        FsrIdx(irsegidx, j) = ireg
        segCenter(:) = half * (CellRay%pts(:, iRaySeg) + CellRay%pts(:, iRaySeg + 1)) - FsrCentroid(:, ireg)
        x0(1, irsegidx, j) = segCenter(1) - LenSeg(iRaySeg) * 0.0005_8 * AziAng(iazi)%cosv
        x0(2, irsegidx, j) = segCenter(1) + LenSeg(iRaySeg) * 0.0005_8 * AziAng(iazi)%cosv
        y0(1, irsegidx, j) = segCenter(2) - LenSeg(iRaySeg) * 0.0005_8 * AziAng(iazi)%sinv
        y0(2, irsegidx, j) = segCenter(2) + LenSeg(iRaySeg) * 0.0005_8 * AziAng(iazi)%sinv
        DO ig = gb, ge
          tau = - LenSeg(iRaySeg) * xst(ig, ireg)
          OptLenList(ig, irsegidx, j) = tau
          q0(ig, irsegidx, j) = srcC(ig, ireg) + sum(srcSlope(1 : 2, ig, ireg) * segCenter(:))
          DO ipol = 1, nPolarAng
            cmOptLen(ipol, ig, irsegidx, j) = - tau * polsininv(ipol) * 0.001_8
            cmOptLenInv(ipol, ig, irsegidx, j) = 1._8 / cmOptLen(ipol, ig, irsegidx, j)
            q1(ipol, ig, irsegidx, j) = half * sum(srcSlope(3 : 4, ig, ireg) * segSlope(:, ipol))
          ENDDO
          ExpAppIdx(ig, irsegidx, j) = max(INT(tau), -40000)
          ExpAppIdx(ig, irsegidx, j) = min(0, ExpAppIdx(ig, irsegidx, j))
        ENDDO
      ENDDO
      CellRayIdxSt(icellrayidx, j, 1) = irsegidx
      SurfIdx(icellRayIdx, j, 1) = AsyRay(iAsyRay)%PinRaySurf(2, l) !OutSurface
      SurfIdx(icellRayIdx, j, 2) = AsyRay(iAsyRay)%PinRaySurf(1, l) !Insurface
    ENDDO
  ENDDO
  nTotRaySeg(j) = irsegidx
  nTotCellRay(j) = icellRayIdx
ENDDO

DO ipol = 1, nPolarAng
  DO j = 1, nCoreRay 
    DO l = 1, nTotRaySeg(j)
      DO ig = gb, ge
        E1(ipol, ig, l, j) = ExpA(ExpAppIdx(ig, l, j), ipol) * OptLenList(ig, l, j) + ExpB(ExpAppIdx(ig, l, j), ipol)
        E3(ipol, ig, l, j) = 2._8 * (cmOptLen(ipol, ig, l, j) - E1(ipol, ig, l, j)) - cmOptLen(ipol, ig, l, j) * E1(ipol, ig, l, j)
        R1(ipol, ig, l, j) = 1._8 + half * cmOptLen(ipol, ig, l, j) - (1._8 + cmOptLenInv(ipol, ig, l, j)) * E1(ipol, ig, l, j)
        R3(ipol, ig, l, j) = rsix * cmOptLen(ipol, ig, l, j) - 2._8 * cmOptLenInv(ipol, ig, l, j) - 2._8  &
                             + (1._8 + cmOptLenInv(ipol, ig, l, j)) * (1._8 + 2._8 * cmOptLenInv(ipol, ig, l, j)) * E1(ipol, ig, l, j)
      ENDDO
    ENDDO
  ENDDO
ENDDO

DO irot = 1, 2
  PhiAnginSvIdx = RayInfo%PhiAngInSvIdx(iRotRay, irot)
  PhiAngOutSvIdx = RayInfo%PhiangOutSvIdx(iRotRay, irot)
  phiobd(1 : nPolarAng, gb : ge) = PhiAngIn(1 : nPolarAng, gb : ge, PhiAnginSvIdx)
  jinc = 1; jbeg = 1; jend = nCoreRay
  IF(irot .eq. 2) THEN
    jinc = -1; jbeg = nCoreRay; jend = 1
  ENDIF
  DO j = jbeg, jend, jinc
    idir = RotRay(i)%DIR(j); iazi = CoreRay(RotRay(irotray)%RayIdx(j))%iang
    IF(irot .eq. 2) idir = mp(idir)
    DO ipol = 1, nPolarAng
      wt(ipol) = wtang(ipol, iazi)
      ax(ipol) = wt(ipol) * AziAng(iazi)%cosv * PolarAng(ipol)%sinv
      ay(ipol) = wt(ipol) * AziAng(iazi)%sinv * PolarAng(ipol)%sinv
    ENDDO
    IF(lJout) THEN
      wt2(1 : nPolarAng, :) = wtang2(1 : nPolarAng, iazi, :)
    ENDIF
    nRaySeg = nTotRaySeg(j)
    IF(idir .eq. 1) THEN
      PhiAngOut(1 : nPolarAng, gb : ge, 1) = phiobd(1 : nPolarAng, gb : ge)
      DO ir = 1, nRaySeg
        ireg = FsrIdx(ir, j)
        DO ig = gb, ge
          DO ipol = 1, nPolarAng          
            phid = (q0(ig, ir, j) - PhiAngOut(ipol, ig, ir)) * E1(ipol, ig, ir, j) + q1(ipol, ig, ir, j) * E3(ipol, ig, ir, j)
            PhiAngOut(ipol, ig, ir + 1) = PhiAngOut(ipol, ig, ir) + phid
            phisC(ig, ireg) = phisC(ig, ireg) - wt(ipol) * phid
            phim = PhiAngOut(ipol, ig, ir) * (cmOptLen(ipol, ig, ir, j) * half) + (q0(ig, ir, j) - PhiAngOut(ipol, ig, ir)) * R1(ipol, ig, ir, j)  &
                   + q1(ipol, ig, ir, j) * cmOptLen(ipol, ig, ir, j) * R3(ipol, ig, ir, j)
            phimx(1, ig, ireg) = phimx(1, ig, ireg) + ax(ipol) * phim * cmOptLen(ipol, ig, ir, j)
            phimy(1, ig, ireg) = phimy(1, ig, ireg) + ay(ipol) * phim * cmOptLen(ipol, ig, ir, j)
            phimx(2, ig, ireg) = phimx(2, ig, ireg) + wt(ipol) * x0(idir, ir, j) * (-phid + q0(ig, ir, j) * cmOptLen(ipol, ig, ir, j))
            phimy(2, ig, ireg) = phimy(2, ig, ireg) + wt(ipol) * y0(idir, ir, j) * (-phid + q0(ig, ir, j) * cmOptLen(ipol, ig, ir, j))
          ENDDO
        ENDDO
      ENDDO
      phiobd(1 : nPolarAng, gb : ge) = PhiAngOut(1 : nPolarAng, gb : ge, nRaySeg + 1)
      IF(ljout) THEN
        DO ir = 1, nTotCellRay(j)
          icel = PinIdx(ir, j); isurf = SurfIdx(ir, j, 1)
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(2, ig, isurf, icel) = Jout(2, ig, isurf, icel) + wt(ipol) * PhiAngOut(ipol, ig, CellRayIdxSt(ir, j, 1) + 1)
              Jout(3, ig, isurf, icel) = Jout(3, ig, isurf, icel) + wt2(ipol, isurf) * PhiAngOut(ipol, ig, CellRayIdxSt(ir, j, 1) + 1)
            ENDDO
          ENDDO
          isurf = SurfIdx(ir, j, 2)
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(1, ig, isurf, icel) = Jout(1, ig, isurf, icel) + wt(ipol) * PhiAngOut(ipol, ig, CellRayIdxSt(ir, j, 2))
              Jout(3, ig, isurf, icel) = Jout(3, ig, isurf, icel) + wt2(ipol, isurf) * PhiAngOut(ipol, ig, CellRayIdxSt(ir, j, 2))
            ENDDO
          ENDDO
        ENDDO
      ENDIF
    ELSE
      PhiAngOut(1 : nPolarAng, gb : ge, nRaySeg + 2) = phiobd(1 : nPolarAng, gb : ge)
      ir = nRaySeg + 1
      DO ir1 = 1, nRaySeg
        ir = ir - 1
        ireg = FsrIdx(ir, j)
        DO ig = gb, ge
          DO ipol = 1, nPolarAng
            phid = (q0(ig, ir, j) - PhiAngOut(ipol, ig, ir + 2)) * E1(ipol, ig, ir, j) - q1(ipol, ig, ir, j) * E3(ipol, ig, ir, j)
            PhiAngOut(ipol, ig, ir + 1) = PhiAngOut(ipol, ig, ir + 2) + phid
            phisC(ig, ireg) = phisC(ig, ireg) - wt(ipol) * phid
            phim = PhiAngOut(ipol, ig, ir + 2) * (cmOptLen(ipol, ig, ir, j) * half) + (q0(ig, ir, j) - PhiAngOut(ipol, ig, ir + 2)) * R1(ipol, ig, ir, j)  &
                   - q1(ipol, ig, ir, j) * cmOptLen(ipol, ig, ir, j) * R3(ipol, ig, ir, j)
            phimx(1, ig, ireg) = phimx(1, ig, ireg) - ax(ipol) * phim * cmOptLen(ipol, ig, ir, j)
            phimy(1, ig, ireg) = phimy(1, ig, ireg) - ay(ipol) * phim * cmOptLen(ipol, ig, ir, j)
            phimx(2, ig, ireg) = phimx(2, ig, ireg) + wt(ipol) * x0(idir, ir, j) * (-phid + q0(ig, ir, j) * cmOptLen(ipol, ig, ir, j))
            phimy(2, ig, ireg) = phimy(2, ig, ireg) + wt(ipol) * y0(idir, ir, j) * (-phid + q0(ig, ir, j) * cmOptLen(ipol, ig, ir, j))
          ENDDO
        ENDDO
      ENDDO
      phiobd(1 : nPolarAng, gb : ge) = PhiAngOut(1 : nPolarAng, gb : ge, 2)
      IF(lJout) THEN
        DO ir = 1, nTotCellRay(j)
          icel = PinIdx(ir, j); isurf = SurfIdx(ir, j, 2)
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(2, ig, isurf, icel) = Jout(2, ig, isurf, icel) + wt(ipol) * PhiAngOut(ipol, ig, CellRayIdxSt(ir, j, 2) + 1)
              Jout(3, ig, isurf, icel) = Jout(3, ig, isurf, icel) + wt2(ipol, isurf) * PhiAngOut(ipol, ig, CellRayIdxSt(ir, j, 2) + 1)
            ENDDO
          ENDDO
          isurf = SurfIdx(ir, j, 1)
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(1, ig, isurf, icel) = Jout(1, ig, isurf, icel) + wt(ipol) * PhiAngOut(ipol, ig, CellRayIdxSt(ir, j, 1) + 2)
              Jout(3, ig, isurf, icel) = Jout(3, ig, isurf, icel) + wt2(ipol, isurf) * PhiAngOut(ipol, ig, CellRayIdxSt(ir, j, 1) + 2)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
    ENDIF
  ENDDO
  PhiAngIn(1 : nPolarAng, gb : ge, PhiAngOutSvIdx) = phiobd(1 : nPolarAng, gb : ge)
ENDDO

NULLIFY(AziAng); NULLIFY(PolarAng)
NULLIFY(AsyRay); NULLIFY(CoreRay)
NULLIFY(RotRay); NULLIFY(CellRay)

!Geometry Info Pointing
NULLIFY(Asy); NULLIFY(Pin)
NULLIFY(PinInfo); NULLIFY(Cell)
NULLIFY(FsrCentroid)

!Tracking Dat Pointing
NULLIFY(FsrIdx); NULLIFY(ExpAppIdx)
NULLIFY(OptLenList); NULLIFY(cmOptLen); NULLIFY(cmOptLenInv)
NULLIFY(LenSeg); NULLIFY(LocalFsrIdx)
NULLIFY(phisC); NULLIFY(phimx); NULLIFY(phimy)
NULLIFY(srcC); NULLIFY(srcSlope)
NULLIFY(xst); NULLIFY(jout)
NULLIFY(PhiAngOut); NULLIFY(PhiAngIn)
NULLIFY(ExpA); NULLIFY(ExpB)
NULLIFY(E1); NULLIFY(E3); NULLIFY(R1); NULLIFY(R3)
NULLIFY(q0); NULLIFY(q1); NULLIFY(x0); NULLIFY(y0)

END SUBROUTINE TrackRotRayLS_CASMO
! ------------------------------------------------------------------------------------------------------------