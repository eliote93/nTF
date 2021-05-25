SUBROUTINE HexSetAsyRay(icBss)

USE PARAM,   ONLY : TRUE, FALSE, MESG
USE PE_MOD,  ONLY : PE
USE IOUTIL,  ONLY : message, terminate
USE FILES,   ONLY : io8
USE HexType, ONLY : Type_HexAsyRay, Type_HexCelRay, Type_HexPinRay
USE HexData, ONLY : NumMray, hmRay, hCelBss, hAsyTypInfo, RayCel, RayPinInfo, haRay, nGeoTyp, hLgc

USE HexAsyRayConst
USE HexAsyRayIntSct

IMPLICIT NONE

INTEGER :: icBss
INTEGER :: imRay, iAng, iGeo, ihpRay, nhpRay, nMaxPinPt
REAL    :: RayEqn(3)
LOGICAL :: lGeo(7) = TRUE

TYPE(Type_HexPinRay), POINTER :: hpRay(:)
TYPE(Type_HexAsyRay), POINTER :: haRay_Loc
TYPE(Type_HexCelRay), POINTER :: CelRay_Loc
! ----------------------------------------------------

IF (hLgc%lSngCel) THEN
  CALL HexSetAsyRay_SngCel
  
  RETURN
ELSE IF (hLgc%lRadVac) THEN
  lGeo(2) = FALSE
  lGeo(5) = FALSE
  lGeo(6) = FALSE
END IF

WRITE(MESG, '(9X, I2, X, A3)') icBss, '...'
IF(PE%Master) CALL message(io8, TRUE, TRUE, MESG)
! ----------------------------------------------------
!               01. INIT
! ----------------------------------------------------
ALLOCATE (RayPinInfo (hAsyTypInfo(hCelBss(icBss)%iaTyp)%nTotPin(1)))

CALL HexSetRayPinDat(RayPinInfo, icBss)

ALLOCATE (RayCel (2))

CALL HexSetRayRodCel_12(RayCel, icBss)
CALL HexSetRayGapCel   (RayCel, icBss)
! ----------------------------------------------------
!               02. ITER
! ----------------------------------------------------
nMaxPinPt = 10 * hCelBss(icBss)%nPin ! ARBITRARY

!$OMP PARALLEL PRIVATE(imRay, hpRay, RayEqn, iGeo, nhpRay, haRay_Loc, ihpRay, CelRay_Loc)
!$OMP DO SCHEDULE(GUIDED)
DO imRay = 1, NumMray(0)
  ALLOCATE (hpRay (nMaxPinPt))
  
  RayEqn = hmRay(imRay)%Eq ! Origin : Asy Cnt
  
  DO iGeo = 1, nGeoTyp
    IF (.NOT. lGeo(iGeo)) CYCLE
    
    CALL HexSetRayIntSct_Pin(icBss, RayEqn, imRay, iGeo, nhpRay, hpRay)
    
    IF (nhpRay > nMaxPinPt) CALL terminate("SET : ASY RAY - NHPRAY")
    
    haRay_Loc => haRay(iGeo, icBss, imRay)
    
    haRay_Loc%nhpRay = nhpRay
    
    ALLOCATE (haRay_Loc%CelRay (nhpRay))
    
    DO ihpRay = 1, nhpRay
      CelRay_Loc => haRay_Loc%CelRay(ihpRay)
      
      CelRay_Loc%hPinIdx      = hpRay(ihpRay)%PinIdx
      CelRay_Loc%hSufIdx(1:2) = hpRay(ihpRay)%SurfIdx(1:2)
      
      CALL HexSetRayIntSct_Msh(CelRay_Loc, RayEqn, iGeo, CelRay_Loc%hPinIdx, hpRay(ihpRay)%PinPt)
    END DO
  END DO
  
  NULLIFY (hpRay)
END DO
!$OMP END DO
!$OMP END PARALLEL

NULLIFY (haRay_Loc, CelRay_Loc, RayPinInfo, RayCel)
! ----------------------------------------------------

END SUBROUTINE HexSetAsyRay