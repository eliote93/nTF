MODULE HexCmfd
  
IMPLICIT NONE

CONTAINS
! ------------------------------------------------------------------------------------------------------------
!                                     01. HEX : Super Pin Currnet
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSuperPinCurrent(Pin, Jout, superJout, ng, nxy, myzb, myze)

USE TYPEDEF_COMMON, ONLY : superPin_Type

IMPLICIT NONE

TYPE(superPin_Type), POINTER :: Pin(:)
REAL, POINTER :: Jout(:, :, :, :, :), superJout(:, :, :, :, :)
INTEGER :: ig, ix, iy, iz, ixy, jxy, ipin, iNgh, iBndy, jBndy
INTEGER :: ng, nxy
INTEGER :: myzb, myze
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
          
          superJout(:, iNgh, ixy, iz, ig) = superJout(:, iNgh, ixy, iz, ig) &
                                               + Jout(:, jBndy, iPin, iz, ig) * ratio
        END DO
      END DO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL
! ----------------------------------------------------

END SUBROUTINE HexSuperPinCurrent
! ------------------------------------------------------------------------------------------------------------
!                                     02. HEX : Ray Info 4 CMFD Gen
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexRayInfo4CmfdGen(RayInfo, RayInfo4Cmfd)

USE TYPEDEF, ONLY : RayInfo_type, RayInfo4Cmfd_Type
USE HexData, ONLY : hRotRay, hcRay, hAsy, hAsyTypInfo, haRay
USE HexUtil, ONLY : SetSgn_INT
USE ALLOCS

IMPLICIT NONE

TYPE(RayInfo_type) :: RayInfo
TYPE(RayInfo4Cmfd_type), POINTER :: RayInfo4Cmfd

INTEGER :: iRotRay, tDir, kDir, isv, icNum, imNum, ipNum, nRotRay, ncRay, nmRay, npRay
INTEGER :: jcRay, jmRay, jAsy, jaTyp, jGeo, jcBss, jPin1, jPin2, jPin3
! ----------------------------------------------------

nRotRay = RayInfo%nRotRay

ALLOCATE (RayInfo4CMFD)

RayInfo4CMFD%nRotRay   = nRotRay
RayInfo4CMFD%nPolAngle = RayInfo%nPolarAngle

CALL DMALLOC(RayInfo4CMFD%RotRayInOutCell, nRotRay, 2)

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

RayInfo%RayInfo4Cmfd       => RayInfo4CMfd
RayInfo4CMFD%PhiAngInSvIdx => RayInfo%PhiAngInSvIdx

END SUBROUTINE HexRayInfo4CmfdGen
! ------------------------------------------------------------------------------------------------------------

END MODULE HexCmfd