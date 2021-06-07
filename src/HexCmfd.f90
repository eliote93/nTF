MODULE HexCmfd
  
IMPLICIT NONE

CONTAINS
! ------------------------------------------------------------------------------------------------------------
!                                     01. HEX : Super Pin Currnet
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSuperPinCurrent(Pin, Jout, superJout, ng, nxy, myzb, myze)

USE TYPEDEF_COMMON, ONLY : superPin_Type

IMPLICIT NONE

TYPE(superPin_Type), POINTER, DIMENSION(:) :: Pin

REAL, POINTER, DIMENSION(:,:,:,:,:) :: Jout, superJout

INTEGER :: ig, ix, iy, iz, ixy, jxy, ipin, iNgh, iBndy, jBndy
INTEGER :: ng, nxy, myzb, myze
REAL    :: ratio
! ----------------------------------------------------

superJout = 0.0

!$OMP PARALLEL PRIVATE(ig, iz, ixy, iBndy, jxy, iPin, jBndy)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO ig = 1, ng
  DO iz = myzb, myze
    DO ixy = 1, nxy
      DO iNgh = 1, Pin(ixy)%nNgh
        iBndy = Pin(ixy)%NghBd(iNgh)
        ratio = Pin(ixy)%NghLgh(iNgh) / Pin(ixy)%BdLength(iBndy)
        
        DO jxy = 1, Pin(ixy)%nBdmPin(iBndy)
          iPin  = Pin(ixy)%BdMPidx(jxy, iBndy) ! MOC Pin
          jBndy = Pin(ixy)%BdMPsuf(jxy, iBndy) ! MOC Suf
          
          superJout(:, iNgh, ixy, iz, ig) = superJout(:, iNgh, ixy, iz, ig) + Jout(:, jBndy, iPin, iz, ig) * ratio
        END DO
      END DO
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL
! ----------------------------------------------------

END SUBROUTINE HexSuperPinCurrent
! ------------------------------------------------------------------------------------------------------------
!                                     02. HEX : Ray Info 4 CMFD Gen
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexRayInfo4CmfdGen(RayInfo, RayInfo4Cmfd)

USE ALLOCS
USE PARAM,   ONLY : FORWARD
USE TYPEDEF, ONLY : RayInfo_type, RayInfo4Cmfd_Type, DcmpAsyRayInfo_Type
USE CNTL,    ONLY : nTracerCntl
USE HexType, ONLY : Type_HexAsyRay
USE HexData, ONLY : hRotRay, hcRay, hAsy, hAsyTypInfo, haRay, nhAsy
USE HexUtil, ONLY : SetSgn_INT

IMPLICIT NONE

TYPE(RayInfo_type) :: RayInfo
TYPE(RayInfo4Cmfd_type), POINTER :: RayInfo4Cmfd

TYPE(DcmpAsyRayInfo_Type), POINTER, DIMENSION(:,:) :: DcmpAsyRay

INTEGER :: iRotRay, tDir, kDir, isv, icNum, imNum, ipNum, nRotRay, ncRay, nmRay, npRay
INTEGER :: jcRay, jmRay, jAsy, jaTyp, jGeo, jcBss, jPin, jPin1, jPin2, jPin3, iCnt, iAsy, iAsyTyp, iGeoTyp, icBss, iAsyRayBeg, iAsyRayEnd

TYPE (Type_HexAsyRay), POINTER :: haRay_Loc
! ----------------------------------------------------

nRotRay = RayInfo%nRotRay

ALLOCATE (RayInfo4CMFD)

RayInfo4CMFD%nRotRay   = nRotRay
RayInfo4CMFD%nPolAngle = RayInfo%nPolarAngle

CALL dmalloc(RayInfo4CMFD%RotRayInOutCell, nRotRay, 2)
! ----------------------------------------------------
DO iRotRay = 1, nRotRay
  ncRay = hRotRay(iRotRay)%ncRay
  
  DO tDir = 1, 2
    isv = RayInfo%PhiangInSvIdx(iRotRay, tDir)
    
    IF (isv .EQ. 1) CYCLE
        
    icNum = (2 - tDir) * 1 &   ! Beginning
          + (tDir - 1) * ncRay ! End
    jcRay = abs(hRotRay(iRotRay)%cRayIdx(icNum))
    
    nmRay = hcRay(jcRay)%nmRay
    kDir  = (2*tDir - 3) * SetSgn_INT(hRotRay(iRotRay)%cRayIdx(icNum))
    imNum = ((1 - kDir)/2) * 1 &   ! Beginning
          + ((kDir + 1)/2) * nmRay ! End
    jmRay = hcRay(jcRay)%mRayIdx(imNum)
    
    jAsy  = hcRay(jcRay)%AsyIdx(imNum)
    jaTyp = hAsy(jAsy)%AsyTyp
    jcBss = hAsyTypInfo(jaTyp)%iBss
    jGeo  = hAsy(jAsy)%GeoTyp
    
    npRay = haRay(jGeo, jcBss, jmRay)%nhpRay
    ipNum = ((1 - kDir)/2) * 1 &   ! Beginning
          + ((kDir + 1)/2) * npRay ! End
    jPin1 = haRay(jGeo, jcBss, jmRay)%CelRay(ipNum)%hPinIdx ! Local Pin Idx in iGeo = 1
    jPin2 = hAsyTypInfo(jaTyp)%PinLocIdx(jGeo, jPin1)       ! Local Pin Idx in iGeo = jGeo
    jPin3 = hAsy(jAsy)%PinIdxSt + jPin2 - 1                 ! Global Pin Idx
    
    RayInfo4CMfd%RotRayInOutCell(iRotRay, tDir) = jPin3
  END DO
END DO
! ----------------------------------------------------
IF (nTracerCntl%lDomainDcmp) THEN
  DcmpAsyRay => RayInfo%DcmpAsyRay
  
  CALL dmalloc(RayInfo4CMFD%DcmpAsyRayInCell, 2, RayInfo%nModRay, nhAsy)
  CALL dmalloc(RayInfo4CMFD%DcmpAsyRayInSurf, 2, RayInfo%nModRay, nhAsy)
  
  RayInfo4CMFD%DcmpAsyRayCount => RayInfo%DcmpAsyRayCount
  
  DO iAsy = 1, nhAsy
    iAsyTyp = hAsy(iAsy)%AsyTyp
    iGeoTyp = hAsy(iAsy)%GeoTyp
    icBss   = hAsyTypInfo(iAsyTyp)%iBss
    
    DO iCnt = 1, RayInfo%DcmpAsyRayCount(iAsy)
      iAsyRayBeg = DcmpAsyRay(iCnt, iAsy)%AsyRayList(1)
      iAsyRayEnd = DcmpAsyRay(iCnt, iAsy)%AsyRayList(DcmpAsyRay(iCnt, iAsy)%nAsyRay)
      
      haRay_Loc => haRay(iGeoTyp, icBss, iAsyRayBeg)
      npRay      = haRay_Loc%nhpRay
      
      IF (DcmpAsyRay(iCnt, iAsy)%DirList(1) .EQ. FORWARD) THEN
        jPin = haRay_Loc%CelRay(1)%hPinIdx
        
        RayInfo4CMFD%DcmpAsyRayInSurf(1, iCnt, iAsy) = haRay_Loc%CelRay(1)%hSufIdx(1)
        RayInfo4CMFD%DcmpAsyRayInCell(1, iCnt, iAsy) = hAsy(iAsy)%PinIdxSt + hAsyTypInfo(iAsyTyp)%PinLocIdx(iGeoTyp, jPin) - 1
      ELSE
        jPin = haRay_Loc%CelRay(npRay)%hPinIdx
        
        RayInfo4CMFD%DcmpAsyRayInSurf(1, iCnt, iAsy) = haRay_Loc%CelRay(npRay)%hSufIdx(2)
        RayInfo4CMFD%DcmpAsyRayInCell(1, iCnt, iAsy) = hAsy(iAsy)%PinIdxSt + hAsyTypInfo(iAsyTyp)%PinLocIdx(iGeoTyp, jPin) - 1
      END IF
      
      haRay_Loc => haRay(iGeoTyp, icBss, iAsyRayEnd)
      npRay     = haRay_Loc%nhpRay
      
      IF (DcmpAsyRay(iCnt, iAsy)%DirList(DcmpAsyRay(iCnt, iAsy)%nAsyRay) .EQ. 1) THEN
        jPin = haRay_Loc%CelRay(npRay)%hPinIdx
        
        RayInfo4CMFD%DcmpAsyRayInSurf(1, iCnt, iAsy) = haRay_Loc%CelRay(npRay)%hsufIdx(2)
        RayInfo4CMFD%DcmpAsyRayInCell(1, iCnt, iAsy) = hAsy(iAsy)%PinIdxSt + hAsyTypInfo(iAsyTyp)%PinLocIdx(iGeoTyp, jPin) - 1
      ELSE
        jPin = haRay_Loc%CelRay(1)%hPinIdx
        
        RayInfo4CMFD%DcmpAsyRayInSurf(1, iCnt, iAsy) = haRay_Loc%CelRay(1)%hSufIdx(1)
        RayInfo4CMFD%DcmpAsyRayInCell(1, iCnt, iAsy) = hAsy(iAsy)%PinIdxSt + hAsyTypInfo(iAsyTyp)%PinLocIdx(iGeoTyp, jPin) - 1
      END IF
    END DO
  END DO
END IF

RayInfo%RayInfo4Cmfd       => RayInfo4CMfd
RayInfo4CMFD%PhiAngInSvIdx => RayInfo%PhiAngInSvIdx
! ----------------------------------------------------

END SUBROUTINE HexRayInfo4CmfdGen
! ------------------------------------------------------------------------------------------------------------
!                                     03. SET : MOC Phi In
! ------------------------------------------------------------------------------------------------------------
#include <defines.h>
#ifdef __INTEL_MKL
SUBROUTINE HexSetMocPhiIn(Core, PinXS)

USE MKL_3D,      ONLY : mklGeom, mklCMFD, superPin_Type
USE PARAM,       ONLY : ZERO
USE TYPEDEF,     ONLY : CoreInfo_Type, Pin_Type, PinXs_Type, RayInfo4CMFD_Type, CmInfo_Type
USE RAYS,        ONLY : RayInfo4CMFD
USE CMFD_MOD,    ONLY : PinNeighIdx, CMFDPinXS, SubPlaneMap
USE CNTL,        ONLY : nTracerCntl
USE itrcntl_mod, ONLY : ItrCntl
USE HexData,     ONLY : hPinInfo

IMPLICIT NONE

TYPE (CoreInfo_Type) :: Core
TYPE (CmInfo_Type)   :: CMInfo

INTEGER :: myzb, myze, ng

TYPE (superPin_Type), POINTER, DIMENSION(:)   :: superPin
TYPE (PinXs_Type),    POINTER, DIMENSION(:,:) :: PinXs

REAL, POINTER, DIMENSION(:,:,:) :: PHIC, PhiFm

INTEGER, POINTER, DIMENSION(:)     :: DcmpAsyRayCount, pinMap
INTEGER, POINTER, DIMENSION(:,:)   :: RotRayInOutCell, PhiangInSvIdx, fmRange
INTEGER, POINTER, DIMENSION(:,:,:) :: DcmpAsyRayInSurf, DcmpAsyRayInCell

INTEGER :: iz, izf, ig, ipol, ixy, isv, ipin, ingh, nRotRay, nCoreRay, nPolar, nAsy, icnt, idir, iRotRay, isxy, jsxy, iAsy, isurf

REAL, POINTER, DIMENSION(:)           :: hz, hzfm
REAL, POINTER, DIMENSION(:,:,:,:)     :: PhiAngIn
REAL, POINTER, DIMENSION(:,:,:,:,:)   :: RadJout
REAL, POINTER, DIMENSION(:,:,:,:,:,:) :: AsyPhiAngIn

REAL :: myphi, nghphi, surfphi, atil, ahat, smy, fmult, phisum
! ----------------------------------------------------

IF (.NOT.nTracerCntl%lHex .OR. .NOT.nTracerCntl%lDomainDcmp) RETURN

nAsy = Core%nxya

pinMap   => mklGeom%pinMap
fmRange  => mklGeom%fmRange
hzfm     => mklGeom%hzfm
hz       => mklGeom%hz
superPin => mklGeom%superPin

nRotRay           = RayInfo4CMFD%nRotRay
nCoreRay          = RayInfo4CMFD%nCoreRay
nPolar            = RayInfo4CMFD%nPolAngle
RotRayInOutCell  => RayInfo4CMFD%RotRayInOutCell
PhiAngInSvIdx    => RayInfo4CMFD%PhiAngInSvIdx
PhiAngIn         => RayInfo4CMFD%PhiAngIn
DcmpAsyRayCount  => RayInfo4CMFD%DcmpAsyRayCount
DcmpAsyRayInSurf => RayInfo4CMFD%DcmpAsyRayInSurf
DcmpAsyRayInCell => RayInfo4CMFD%DcmpAsyRayInCell
AsyPhiAngIn      => RayInfo4CMFD%AsyPhiAngIn
! ----------------------------------------------------
DO idir = 1, 2
  DO iRotRay = 1, nRotRay
    ixy = RotRayInOutCell(iRotRay, idir) ! MoC
    
    IF (ixy .EQ. 0) CYCLE
    
    isxy = hPinInfo(ixy)%ihcPin ! Suer-Pin
    isv  = PhiAngInSvIdx(iRotRay, idir)
    
    DO iz = myzb, myze
      DO ig = 1, ng
        phisum = ZERO
        
        DO izf = fmRange(iz, 1), fmRange(iz, 2)
          phisum = phisum + mklCMFD%phis(isxy, izf, ig) * (hzfm(izf) / hz(iz))
        END DO
        
        fmult = phisum / PinXS(ixy, iz)%phi(ig)
        
        DO ipol = 1, nPolar
          PhiAngIn(ipol, isv, iz, ig) = PhiAngIn(ipol, isv, iz, ig) * fmult
        END DO
      END DO
    END DO
  END DO
END DO
! ----------------------------------------------------
IF (nTracerCntl%lDomainDcmp) THEN
  DO iz = myzb, myze
    DO iAsy = 1, nAsy
      DO icnt = 1, DcmpAsyRayCount(iAsy)
        DO idir = 1, 2
          ixy   = DcmpAsyRayInCell(idir, icnt, iAsy)
          isurf = DcmpAsyRayInSurf(idir, icnt, iAsy)
          isxy  = hPinInfo(ixy)%ihcPin
          ingh  = hPinInfo(ixy)%DcmpMP2slfSPngh(isurf)
          jsxy  = hPinInfo(ixy)%DcmpMP2nghSPidx(isurf)
          smy   = superPin(isxy)%BdLength(ingh)
          
          DO ig = 1, ng
            myphi = ZERO
            DO izf = fmRange(iz, 1), fmRange(iz, 2)
              myphi = myphi + mklCMFD%phis(isxy, izf, ig) * (hzfm(izf) / hz(iz))
            END DO
            
            nghphi = ZERO
            IF (jsxy .GT. 0) THEN
              DO izf = fmRange(iz, 1), fmRange(iz, 2)
                nghphi = nghphi + mklCMFD%phis(ingh, izf, ig) * (hzfm(izf) / hz(iz))
              END DO
            END IF
            
            atil = PinXS(ixy, iz)%atil(ingh, ig)
            ahat = PinXS(ixy, iz)%ahat(ingh, ig)
            
            surfphi = atil * myphi + (smy - atil) * nghphi
            
            IF (ItrCntl%mocit .EQ. 0) THEN
              AsyPhiAngIn(:, ig, idir, icnt, iAsy, iz) = surfphi / smy
            ELSE
              surfphi = surfphi + ahat * (myphi + nghphi)
              fmult   = surfphi / RadJout(3, isurf, ixy, iz, ig)
              
              DO ipol = 1, nPolar
                AsyPhiAngIn(ipol, ig, idir, icnt, iAsy, iz) = AsyPhiAngIn(ipol, ig, idir, icnt, iAsy, iz) * fmult
              END DO
            END IF
          END DO
        END DO
      END DO
    END DO
  END DO
END IF
! ----------------------------------------------------
! Dcmp.
NULLIFY (RotRayInOutCell)
NULLIFY (PhiAngInSvIdx)
NULLIFY (PhiAngIn)
NULLIFY (DcmpAsyRayCount)
NULLIFY (DcmpAsyRayInSurf)
NULLIFY (DcmpAsyRayInCell)
NULLIFY (AsyPhiAngIn)

! Geo
NULLIFY (superPin)
NULLIFY (RayInfo4Cmfd)
NULLIFY (PhiC)
NULLIFY (PhiFm)
NULLIFY (PinXS)
NULLIFY (RadJout)
! ----------------------------------------------------

END SUBROUTINE HexSetMocPhiIn
#endif
! ------------------------------------------------------------------------------------------------------------

END MODULE HexCmfd