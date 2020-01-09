
SUBROUTINE SetCmfdLinearSystem(lcmfd, l3dim, AxSolver)
USE PARAM
USE TYPEDEF, ONLY : CMFDLS_TYPE
USE itrcntl_mod, ONLY : itrcntl
USE CMFD_MOD, ONLY : ng,        nxy,          myzb,         myze,                &
                     myzbf,     myzef,                                           &
                     hzfm,      hz,           SubPlaneMap,  SubPlaneRange,       &
                     CmfdPinXS, PinNeighIdx,  PinVol,       PinVolFm,            &
                     CmfdLs,    CmfdLs1g,                                        &
                     AxDhat,    AxDtil,       AxPdhat
USE BASICOPERATION, ONLY : CP_CA
IMPLICIT NONE
LOGICAL :: lcmfd, l3dim
INTEGER :: AxSolver
INTEGER :: ig, iz, iz0, ixy, ibd, iz1, iz2
INTEGER, PARAMETER :: nbd = 4
REAL, POINTER :: Diag(:, :), RadOffDiag(:, :, :), AxOffDiag(:, :, :)
REAL :: Vol, neighVol, lmnt, offlmnt, dhat, dtil, alpha, dhmax
CHARACTER(2) :: Step
dhmax = 0
DO IG = 1, ng
  !CALL CP_CA(Core)
  CMFDLS(ig)%AxialPlaneMap => SubPlaneMap
  Diag => CmfdLs(ig)%Diag
  RadOffDiag => CmfdLs(ig)%RadOffDiag
  AxOffDiag => CmfdLs(ig)%AxOffDiag
  CALL CP_CA(diag(:, myzbf:myzef), zero, nxy, myzef-myzbf + 1)
  CALL CP_CA(RadOffDiag(:, :, myzb:myze), zero, 4, nxy, myze-myzb + 1)
  CALL CP_CA(AxOffDiag(:, :, myzbf:myzef), zero, 2, nxy, myzef-myzbf + 1)
  DO iz = myzb, myze
    iz1 = SubPlaneRange(1, iz); iz2 = SubPlaneRange(2, iz); 
    DO ixy  = 1, nxy
      !Vol = PinVol(ixy, iz)
      Vol = PinVolFm(ixy, iz)
      !Diagonal Component
      lmnt = Vol*(CmfdPinXS(ixy, iz)%xsr(ig))
      DO ibd = 1, nbd
        dtil = CmfdPinXs(ixy, iz)%Dtil(ibd, ig)
        dhat = CmfdPinXs(ixy, iz)%Dhat(ibd, ig)
        lmnt = lmnt + (dtil - dhat) * hzfm(iz1)
        offlmnt = -(dtil + dhat) * hzfm(iz1)
        RadOffDiag(ibd, ixy, iz) = offlmnt  
      ENDDO
      Diag(ixy, iz1:iz2) = lmnt
    ENDDO
  ENDDO
  IF(.NOT. l3dim) CYCLE
  DO iz = myzbf, myzef
    iz0 = SubPlaneMap(iz)
    DO ixy = 1, nxy
      vol = PinVolFm(ixy, iz0)/hzfm(iz)
      dtil = AxDtil(1, ixy, iz, ig)
      dhat = AxDhat(1, ixy, iz, ig)
      !Under-Relaxation
      dhmax = max(abs(dhat), dhmax)
      !IF(abs(Dhat-AxFlx(iz, ixy)%Pdhat(1, ig))>10.*Dtil) THEN
      !  ALPHA = DTIL/(ABS(DHAT-AxFlx(iz, ixy)%Pdhat(1, ig))-DTIL)
        !DHAT = DHAT * ALPHA + AxFlx(iz, ixy)%Pdhat(1, ig)*(1.-ALPHA)
      !ENDIF
      AxPDhat(1, ixy, iz, ig) = Dhat
      lmnt = (dtil - dhat)
      offlmnt = -(dtil + dhat)*vol
      AxOffDiag(1, ixy, iz) = OffLmnt
      dhmax = max(abs(dhat), dhmax)
      dtil = AxDtil(2, ixy, iz, ig); dhat = AxDhat(2, ixy, iz, ig)
      
      !dtil = AxFlx(iz, ixy)%dtil(2, ig); dhat = AxFlx(iz, ixy)%dhat(2, ig)
     ! IF(abs(Dhat-AxFlx(iz, ixy)%Pdhat(2, ig))>10.*Dtil) THEN
      !  ALPHA = DTIL/(ABS(DHAT-AxFlx(iz, ixy)%Pdhat(2, ig))-DTIL)
        !DHAT = DHAT * ALPHA + AxFlx(iz, ixy)%Pdhat(2, ig)*(1.-ALPHA)
      !ENDIF
      AxPDhat(2, ixy, iz, ig) = Dhat
      !AxFlx(iz, ixy)%Pdhat(2, ig) = Dhat
      lmnt = lmnt + (dtil - dhat)
      offlmnt = -(dtil + dhat)*vol
      AxOffDiag(2, ixy, iz) = OffLmnt
      lmnt = lmnt * vol 
      Diag(ixy, iz) = Diag(ixy, iz) + lmnt
    ENDDO
  ENDDO
ENDDO
!print *, dhmax
NULLIFY(Diag, RadOffDiag, AxOffDiag)
!#define PRINT_LS
#ifdef PRINT_LS
WRITE(Step, '(I2)'), itrcntl%mocit
OPEN(66, FILE = 'CMFD ' // trim(Step) // '.txt')
DO ig = 1, ng
  WRITE(66, *), "Group ", ig
  DO iz = myzb, myze
    DO ixy = 1, nxy
      IF (l3dim) WRITE(66, *), CMFDLS(ig)%AxOffDiag(2, ixy, iz)
      WRITE(66, '(5F10.4)'), CMFDLS(ig)%RadOffDiag(NORTH, ixy, iz), CMFDLS(ig)%RadOffDiag(WEST, ixy, iz), CMFDLS(ig)%Diag(ixy, iz), &
                             CMFDLS(ig)%RadOffDiag(EAST, ixy, iz), CMFDLS(ig)%RadOffDIag(SOUTH, ixy, iz)
      IF (l3dim) WRITE(66, *), CMFDLS(ig)%AxOffDiag(1, ixy, iz)
    ENDDO
  ENDDO
ENDDO
CLOSE(66)
STOP
#endif

END SUBROUTINE
