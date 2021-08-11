MODULE HexCmfd
  
IMPLICIT NONE

CONTAINS
! ------------------------------------------------------------------------------------------------------------
!                                     01. HEX : Super Pin Currnet
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSuperPinCurrent(Pin, RadJout, superJout, ng, nxy, myzb, myze)

USE TYPEDEF_COMMON, ONLY : superPin_Type

IMPLICIT NONE

TYPE(superPin_Type), POINTER, DIMENSION(:) :: Pin

REAL, POINTER, DIMENSION(:,:,:,:,:) :: RadJout, superJout

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
          
          superJout(:, iNgh, ixy, iz, ig) = superJout(:, iNgh, ixy, iz, ig) + RadJout(:, jBndy, iPin, iz, ig) * ratio
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
! ----------------------------------------------------
TYPE(DcmpAsyRayInfo_Type), POINTER, DIMENSION(:,:) :: DcmpAsyRay

INTEGER :: iRotRay, tDir, kDir, isv, icNum, imNum, ipNum, nRotRay, ncRay, nmRay, npRay
INTEGER :: jcRay, jmRay, jAsy, jaTyp, jGeo, jcBss, jPin, jPin1, jPin2, jPin3, iCnt, iAsy, iAsyTyp, iGeoTyp, icBss, imRaySt, imRayEd

TYPE (Type_HexAsyRay), POINTER :: haRay_Loc

INTEGER, POINTER, DIMENSION(:,:,:) :: DcmpAsyRayInCell, DcmpAsyRayInSurf
! ----------------------------------------------------

nRotRay = RayInfo%nRotRay

ALLOCATE (RayInfo4CMFD)

RayInfo4CMFD%nRotRay   = nRotRay
RayInfo4CMFD%nPolAngle = RayInfo%nPolarAngle

RayInfo4CMFD%PhiAngInSvIdx => RayInfo%PhiAngInSvIdx

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
  CALL dmalloc(RayInfo4CMFD%DcmpAsyRayInCell, 2, RayInfo%nModRay, nhAsy)
  CALL dmalloc(RayInfo4CMFD%DcmpAsyRayInSurf, 2, RayInfo%nModRay, nhAsy)
  
  RayInfo4CMFD%DcmpAsyRayCount => RayInfo%DcmpAsyRayCount
  
  DcmpAsyRay       => RayInfo%DcmpAsyRay
  DcmpAsyRayInCell => RayInfo4CMFD%DcmpAsyRayInCell
  DcmpAsyRayInSurf => RayInfo4CMFD%DcmpAsyRayInSurf
  
  DO iAsy = 1, nhAsy
    iAsyTyp = hAsy(iAsy)%AsyTyp
    iGeoTyp = hAsy(iAsy)%GeoTyp
    icBss   = hAsyTypInfo(iAsyTyp)%iBss
    
    DO iCnt = 1, RayInfo%DcmpAsyRayCount(iAsy)
      imRaySt = DcmpAsyRay(iCnt, iAsy)%AsyRayList(1)
      imRayEd = DcmpAsyRay(iCnt, iAsy)%AsyRayList(DcmpAsyRay(iCnt, iAsy)%nAsyRay)
      
      ! In-coming Surf.
      haRay_Loc => haRay(iGeoTyp, icBss, imRaySt)
      npRay      = haRay_Loc%nhpRay
      
      IF (DcmpAsyRay(iCnt, iAsy)%DirList(1) .EQ. FORWARD) THEN
        jPin = haRay_Loc%CelRay(1)%hPinIdx
        
        DcmpAsyRayInSurf(1, iCnt, iAsy) = haRay_Loc%CelRay(1)%hSufIdx(1)
        DcmpAsyRayInCell(1, iCnt, iAsy) = hAsy(iAsy)%PinIdxSt + hAsyTypInfo(iAsyTyp)%PinLocIdx(iGeoTyp, jPin) - 1
      ELSE
        jPin = haRay_Loc%CelRay(npRay)%hPinIdx
        
        DcmpAsyRayInSurf(1, iCnt, iAsy) = haRay_Loc%CelRay(npRay)%hSufIdx(2)
        DcmpAsyRayInCell(1, iCnt, iAsy) = hAsy(iAsy)%PinIdxSt + hAsyTypInfo(iAsyTyp)%PinLocIdx(iGeoTyp, jPin) - 1
      END IF
      
      ! Out-going Surf.
      haRay_Loc => haRay(iGeoTyp, icBss, imRayEd)
      npRay      = haRay_Loc%nhpRay
      
      IF (DcmpAsyRay(iCnt, iAsy)%DirList(DcmpAsyRay(iCnt, iAsy)%nAsyRay) .EQ. 1) THEN
        jPin = haRay_Loc%CelRay(npRay)%hPinIdx
        
        DcmpAsyRayInSurf(2, iCnt, iAsy) = haRay_Loc%CelRay(npRay)%hsufIdx(2)
        DcmpAsyRayInCell(2, iCnt, iAsy) = hAsy(iAsy)%PinIdxSt + hAsyTypInfo(iAsyTyp)%PinLocIdx(iGeoTyp, jPin) - 1
      ELSE
        jPin = haRay_Loc%CelRay(1)%hPinIdx
        
        DcmpAsyRayInSurf(2, iCnt, iAsy) = haRay_Loc%CelRay(1)%hSufIdx(1)
        DcmpAsyRayInCell(2, iCnt, iAsy) = hAsy(iAsy)%PinIdxSt + hAsyTypInfo(iAsyTyp)%PinLocIdx(iGeoTyp, jPin) - 1
      END IF
    END DO
  END DO
  
  NULLIFY (DcmpAsyRay)
  NULLIFY (DcmpAsyRayInCell)
  NULLIFY (DcmpAsyRayInSurf)
END IF
! ----------------------------------------------------

END SUBROUTINE HexRayInfo4CmfdGen
! ------------------------------------------------------------------------------------------------------------
!                                     03. SET : MOC Phi In
! ------------------------------------------------------------------------------------------------------------
#include <defines.h>
#ifdef __INTEL_MKL
SUBROUTINE HexSetMocPhiIn(Core, Pin, superJout, PinXS, ng, nxy, myzb, myze, ItrCntl, nTracerCntl)

USE MKL_3D,      ONLY : mklGeom, mklCMFD, superPin_Type
USE PARAM,       ONLY : ZERO, HALF, RFOUR
USE TYPEDEF,     ONLY : CoreInfo_Type, PinXs_Type
USE geom,        ONLY : ncbd
USE ioutil,      ONLY : terminate
USE RAYS,        ONLY : RayInfo4CMFD
USE CNTL,        ONLY : nTracerCntl_Type
USE itrcntl_mod, ONLY : ItrCntl_TYPE
USE HexData,     ONLY : hPinInfo, hLgc

IMPLICIT NONE

TYPE (CoreInfo_Type)     :: Core
TYPE (ItrCntl_TYPE)      :: ItrCntl
TYPE (nTracerCntl_TYPE)  :: nTracerCntl

REAL, POINTER, DIMENSION(:,:,:,:,:) :: superJout

INTEGER :: ng, nxy, myzb, myze

TYPE (superPin_Type), POINTER, DIMENSION(:)   :: Pin
TYPE (PinXs_Type),    POINTER, DIMENSION(:,:) :: PinXs
! ----------------------------------------------------
INTEGER :: iz, izf, ig, imxy, isv, ingh, icnt, idir, iRotRay, isxy, jsxy, iAsy, isurf
INTEGER :: nRotRay, nAsy, nPolarAng

INTEGER, POINTER, DIMENSION(:)     :: DcmpAsyRayCount
INTEGER, POINTER, DIMENSION(:,:)   :: RotRayInOutCell, PhiangInSvIdx, fmRange
INTEGER, POINTER, DIMENSION(:,:,:) :: DcmpAsyRayInSurf, DcmpAsyRayInCell

REAL :: myphi, nghphi, surfphi, atil, ahat, slgh, fmult, phisum

REAL, POINTER, DIMENSION(:)           :: hz, hzfm
REAL, POINTER, DIMENSION(:,:,:)       :: phis
REAL, POINTER, DIMENSION(:,:,:,:)     :: PhiAngIn
REAL, POINTER, DIMENSION(:,:,:,:,:,:) :: AsyPhiAngIn

TYPE (superPin_Type), POINTER, DIMENSION(:) :: superPin
! ----------------------------------------------------

IF (.NOT.nTracerCntl%lHex .OR. .NOT.nTracerCntl%lDomainDcmp) RETURN
IF (.NOT.hLgc%lRadRef .AND. .NOT.nTracerCntl%lDomainDcmp)    RETURN

nAsy = Core%nxya

fmRange  => mklGeom%fmRange
hzfm     => mklGeom%hzfm
hz       => mklGeom%hz
superPin => mklGeom%superPin

phis => mklCMFD%phis

nPolarAng         = RayInfo4CMFD%nPolAngle
nRotRay           = RayInfo4CMFD%nRotRay
RotRayInOutCell  => RayInfo4CMFD%RotRayInOutCell
PhiAngInSvIdx    => RayInfo4CMFD%PhiAngInSvIdx
PhiAngIn         => RayInfo4CMFD%PhiAngIn
DcmpAsyRayCount  => RayInfo4CMFD%DcmpAsyRayCount
DcmpAsyRayInSurf => RayInfo4CMFD%DcmpAsyRayInSurf
DcmpAsyRayInCell => RayInfo4CMFD%DcmpAsyRayInCell
AsyPhiAngIn      => RayInfo4CMFD%AsyPhiAngIn
! ----------------------------------------------------
IF (hLgc%lRadRef) THEN
  !$OMP PARALLEL PRIVATE(idir, iRotRay, imxy, isxy, isv, iz, ig, phisum, izf, fmult)
  !$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
  DO idir = 1, 2
    DO iRotRay = 1, nRotRay
      imxy = RotRayInOutCell(iRotRay, idir) ! Global Idx. of MoC Pin
      
      IF (imxy .EQ. 0) CYCLE
      
      isxy = hPinInfo(imxy)%ihcPin          ! Global Idx. of Suer-Pin
      isv  = PhiAngInSvIdx(iRotRay, idir)
      
      DO iz = myzb, myze
        DO ig = 1, ng
          phisum = ZERO
          
          DO izf = fmRange(iz, 1), fmRange(iz, 2)
            phisum = phisum + phis(isxy, izf, ig) * (hzfm(izf) / hz(iz))
          END DO
          
          fmult = phisum / PinXS(isxy, iz)%phi(ig)
          
          PhiAngIn(:, isv, ig, iz) = PhiAngIn(:, isv, ig, iz) * fmult
        END DO
      END DO
    END DO
  END DO
  !$OMP END DO
  !$OMP END PARALLEL
END IF
! ----------------------------------------------------
IF (nTracerCntl%lDomainDcmp) THEN
  !$OMP PARALLEL PRIVATE(iz, iAsy, idir, icnt, imxy, isurf, isxy, ingh, jsxy, slgh, ig, myphi, izf, nghphi, atil, ahat, surfphi, fmult)
  !$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
  DO iz = myzb, myze
    DO iAsy = 1, nAsy
      DO idir = 1, 2
        DO icnt = 1, DcmpAsyRayCount(iAsy)
          imxy  = DcmpAsyRayInCell(idir, icnt, iAsy) ! Global Idx. of MoC Pin
          isurf = DcmpAsyRayInSurf(idir, icnt, iAsy)
          isxy  = hPinInfo(imxy)%ihcPin              ! Global Idx. of Super-Pin
          ingh  = hPinInfo(imxy)%DcmpMP2slfSPngh(isurf)
          jsxy  = hPinInfo(imxy)%DcmpMP2nghSPidx(isurf)
                    
          slgh  = superPin(isxy)%BdLength(ingh)
          
          DO ig = 1, ng
            myphi = ZERO
            DO izf = fmRange(iz, 1), fmRange(iz, 2)
              myphi = myphi + phis(isxy, izf, ig) * (hzfm(izf) / hz(iz))
            END DO
            
            nghphi = ZERO
            IF (jsxy .GT. 0) THEN
              DO izf = fmRange(iz, 1), fmRange(iz, 2)
                nghphi = nghphi + phis(jsxy, izf, ig) * (hzfm(izf) / hz(iz))
              END DO
            END IF
            
            atil = PinXS(isxy, iz)%atil(ingh, ig)
            ahat = PinXS(isxy, iz)%ahat(ingh, ig)
                        
            surfphi = atil * myphi + (slgh - atil) * nghphi
            
            IF (ItrCntl%mocit .EQ. 0) THEN
              IF (nTracerCntl%lNodeMajor) THEN
                AsyPhiAngIn(1:nPolarAng, ig, idir, icnt, iAsy, iz) = surfphi / slgh
              ELSE
                AsyPhiAngIn(1:nPolarAng, idir, icnt, iAsy, ig, iz) = surfphi / slgh
              END IF
            ELSE
              surfphi = surfphi + ahat * (myphi + nghphi)
              fmult   = surfphi / superJout(3, ingh, isxy, iz, ig)
              
              IF (nTracerCntl%lNodeMajor) THEN
                AsyPhiAngIn(1:nPolarAng, ig, idir, icnt, iAsy, iz) = AsyPhiAngIn(:, ig, idir, icnt, iAsy, iz) * fmult
              ELSE
                AsyPhiAngIn(1:nPolarAng, idir, icnt, iAsy, ig, iz) = AsyPhiAngIn(:, idir, icnt, iAsy, ig, iz) * fmult
              END IF
            END IF
          END DO
        END DO
      END DO
    END DO
  END DO
  !$OMP END DO
  !$OMP END PARALLEL
END IF
! ----------------------------------------------------
! Dcmp.
NULLIFY (DcmpAsyRayCount)
NULLIFY (DcmpAsyRayInSurf)
NULLIFY (DcmpAsyRayInCell)

! Loc.
NULLIFY (RotRayInOutCell)
NULLIFY (PhiAngInSvIdx)
NULLIFY (phis)
NULLIFY (PhiAngIn)
NULLIFY (AsyPhiAngIn)
NULLIFY (superPin)

! Geo
NULLIFY (fmRange)
NULLIFY (hz)
NULLIFY (hzfm)
! ----------------------------------------------------

END SUBROUTINE HexSetMocPhiIn
#endif
! ------------------------------------------------------------------------------------------------------------

END MODULE HexCmfd