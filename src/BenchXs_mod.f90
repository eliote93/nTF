module BenchXs
USE PARAM
USE TYPEDEF,     only : BenchXS_type, XsMac_Type
IMPLICIT NONE
SAVE
INTEGER :: ngben, nxsltype, benchxstype, scatod
INTEGER :: nxschangeben, nxsltype_perturb
!INTEGER, PRIVATE :: i, j, k, ig
TYPE(BenchXs_type), POINTER :: MacXsBen(:)
TYPE(BenchXs_type), PRIVATE, POINTER :: MacXsBen0(:)

LOGICAL :: lKinParam = .FALSE.
LOGICAL :: lTranInit = .FALSE. 
INTEGER :: DNP_NGRP
REAL, POINTER :: DNP_BETA(:)
REAL, POINTER :: DNP_Lambda(:)
REAL, POINTER :: Neut_Velo(:)
REAL, POINTER :: Dneut_Chi(:)

CONTAINS


SUBROUTINE XsBenChange_new(iso0, iso1, isonew, wt)
IMPLICIT NONE
INTEGER, INTENT(IN) :: iso0, iso1
INTEGER, INTENT(INOUT) :: isonew
REAL, INTENT(IN) :: wt

INTEGER :: ig, ig0
REAL :: wtbar 
IF(.NOT. lTranInit) CALL InitXsChange()
IF(isonew .EQ. 0) THEN
  nxsltype_perturb = nxsltype_perturb + 1
  isonew = nxsltype_perturb
  CALL MakenCopyXsBen(isonew, iso0, .TRUE.)
ENDIF

wtbar = 1._8 - wt
DO ig = 1, ngben
  MacXsBen(isonew)%xst(ig) =  MacXsBen(isonew)%xst(ig) * wtbar + MacXsBen(iso1)%xst(ig) * wt
  MacXsBen(isonew)%xstr(ig) =  MacXsBen(isonew)%xstr(ig) * wtbar + MacXsBen(iso1)%xstr(ig) * wt
  MacXsBen(isonew)%xskf(ig) =  MacXsBen(isonew)%xskf(ig) * wtbar + MacXsBen(iso1)%xskf(ig) * wt
  MacXsBen(isonew)%xsnf(ig) =  MacXsBen(isonew)%xsnf(ig) * wtbar + MacXsBen(iso1)%xsnf(ig) * wt
  MacXsBen(isonew)%chi(ig) =  MacXsBen(isonew)%chi(ig) * wtbar + MacXsBen(iso1)%chi(ig) * wt
  MacXsBen(isonew)%xsa(ig) =  MacXsBen(isonew)%xsa(ig) * wtbar + MacXsBen(iso1)%xsa(ig) * wt
  DO ig0 = 1, ngben
    MacXsBen(isonew)%xss(ig0, ig) =  MacXsBen(isonew)%xss(ig0, ig) * wtbar + MacXsBen(iso1)%xss(ig0, ig) * wt
  ENDDO
ENDDO

END SUBROUTINE 
SUBROUTINE XsBenChange(iso0, iso1, wt)
IMPLICIT NONE
INTEGER, INTENT(IN) :: iso0, iso1
REAL, INTENT(IN) :: wt

INTEGER :: ig, ig0
REAL :: wtbar 
IF(.NOT. lTranInit) CALL InitXsChange()

wtbar = 1._8 - wt
DO ig = 1, ngben
  MacXsBen(iso0)%xst(ig) =  MacXsBen(iso0)%xst(ig) * wtbar + MacXsBen0(iso1)%xst(ig) * wt
  MacXsBen(iso0)%xstr(ig) =  MacXsBen(iso0)%xstr(ig) * wtbar + MacXsBen0(iso1)%xstr(ig) * wt
  MacXsBen(iso0)%xskf(ig) =  MacXsBen(iso0)%xskf(ig) * wtbar + MacXsBen0(iso1)%xskf(ig) * wt
  MacXsBen(iso0)%xsnf(ig) =  MacXsBen(iso0)%xsnf(ig) * wtbar + MacXsBen0(iso1)%xsnf(ig) * wt
  MacXsBen(iso0)%chi(ig) =  MacXsBen(iso0)%chi(ig) * wtbar + MacXsBen0(iso1)%chi(ig) * wt
  MacXsBen(iso0)%xsa(ig) =  MacXsBen(iso0)%xsa(ig) * wtbar + MacXsBen0(iso1)%xsa(ig) * wt
  DO ig0 = 1, ngben
    MacXsBen(iso0)%xss(ig0, ig) =  MacXsBen(iso0)%xss(ig0, ig) * wtbar + MacXsBen0(iso1)%xss(ig0, ig) * wt
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE MakenCopyXsBen(iso0, iso1, lalloc)
USE ALLOCS
IMPLICIT NONE
INTEGER :: iso0, iso1
LOGICAL :: lalloc
IF(lalloc) THEN
  MacXsBen(iso0)%lempty = MacXsBen(iso1)%lempty;     MacXsBen(iso0)%lFuel = MacXsBen(iso1)%lFuel
  CALL Dmalloc(MacXsBen(iso0)%xskf, ngben); CALL Dmalloc(MacXsBen(iso0)%xsnf, ngben);
  CALL Dmalloc(MacXsBen(iso0)%chi, ngben);  CALL Dmalloc(MacXsBen(iso0)%xsa, ngben);
  CALL Dmalloc(MacXsBen(iso0)%xss, ngben, ngben);
  MacXsBen(iso0)%Xss => MacXsBen(iso0)%Xss
ENDIF

MacXsBen(iso0)%xst = MacXsBen(iso1)%xst;   MacXsBen(iso0)%xstr = MacXsBen(iso1)%xstr
MacXsBen(iso0)%xskf = MacXsBen(iso1)%xskf;   MacXsBen(iso0)%xsnf = MacXsBen(iso1)%xsnf
MacXsBen(iso0)%chi = MacXsBen(iso1)%chi;   MacXsBen(iso0)%xsa = MacXsBen(iso1)%xsa
MacXsBen(iso0)%xss = MacXsBen(iso1)%xss

END SUBROUTINE

SUBROUTINE InitXsChange()
USE ALLOCS 
IMPLICIT NONE
INTEGER :: ixsl

IF(lTranInit) RETURN
nxsltype_perturb = nXslType
ALLOCATE(MacXsBen0(nXslType))
DO ixsl = 1, nXslType
  MacXsBen0(ixsl)%lempty = MacXsBen(ixsl)%lempty;     MacXsBen0(ixsl)%lFuel = MacXsBen(ixsl)%lFuel
  CALL Dmalloc(MacXsBen0(ixsl)%xst, ngben);  CALL Dmalloc(MacXsBen0(ixsl)%xstr, ngben);    !Allocate 
  CALL Dmalloc(MacXsBen0(ixsl)%xskf, ngben); CALL Dmalloc(MacXsBen0(ixsl)%xsnf, ngben);
  CALL Dmalloc(MacXsBen0(ixsl)%chi, ngben);  CALL Dmalloc(MacXsBen0(ixsl)%xsa, ngben);
  CALL Dmalloc(MacXsBen0(ixsl)%xss, ngben, ngben);
  
  MacXsBen(ixsl)%Xss0 => MacXsBen(ixsl)%Xss
  
  MacXsBen0(ixsl)%xst = MacXsBen(ixsl)%xst;   MacXsBen0(ixsl)%xstr = MacXsBen(ixsl)%xstr
  MacXsBen0(ixsl)%xskf = MacXsBen(ixsl)%xskf;   MacXsBen0(ixsl)%xsnf = MacXsBen(ixsl)%xsnf
  MacXsBen0(ixsl)%chi = MacXsBen(ixsl)%chi;   MacXsBen0(ixsl)%xsa = MacXsBen(ixsl)%xsa
  MacXsBen0(ixsl)%xss = MacXsBen(ixsl)%xss
ENDDO
lTranInit = .TRUE.
END SUBROUTINE


SUBROUTINE xsbaseBen(itype, igb, ige, jsfr, jsto, lscat1, XsMac)
!Base Macroscopic XS
IMPLICIT NONE
INTEGER :: i, j, k, ig
INTEGER :: itype, igb, ige, jsfr, jsto
LOGICAL :: lscat1
LOGICAL :: lscat1_
!INTEGER :: scatod
TYPE(XsMac_Type) :: XsMac

lscat1_=lscat1
IF( scatod .GT. 0 )THEN
    lscat1_=.TRUE. .AND. lscat1_
ENDIF

XsMac%lFuel = MacXsBen(itype)%lFuel
DO i = igb, ige
  XsMac%xsmact(i) = MacXsBen(itype)%xst(i)
  XsMac%xsmaca(i) = MacXsBen(itype)%xsa(i)
  XsMac%xsmactr(i) = MacXsBen(itype)%xstr(i)
  XsMac%XsMacnf(i) = MacXsBen(itype)%xsnf(i)
  XsMac%XsMackf(i) = MacXsBen(itype)%xskf(i)
  XsMac%chi(i) = MacXsBen(itype)%chi(i)
ENDDO

IF( benchxstype .eq. 1 )THEN
DO i = igb, ige
    DO j = jsfr, jsto
      XsMac%xsmacsm(j, i) = MacXsBen(itype)%xss(j, i)
    ENDDO
    XsMac%xsmacs(i)=0
    DO j = igb, ige
        XsMac%xsmacs(i) = XsMac%xsmacs(i)+ MacXsBen(itype)%xss(i, j)
    ENDDO
    XsMac%xsmacstr(i)=XsMac%xsmacs(i)
    DO j = jsfr, jsto
      IF (ScatOd .GE. 1) THEN
        XsMac%xsmacp1sm(j, i) = MacXsBen(itype)%xss1(j, i)
        IF (ScatOd .GE. 2) THEN
          XsMac%xsmacp2sm(j, i) = MacXsBen(itype)%xss2(j, i)
          IF (ScatOd .GE. 3) THEN
            XsMac%xsmacp3sm(j, i) = MacXsBen(itype)%xss3(j, i)
          ENDIF
        ENDIF
      ENDIF
    ENDDO
    !IF( lscat1_ )THEN
    !    DO j = jsfr, jsto
    !        XsMac%xsmacsm(j,j) = XsMac%xsmacsm(j,j) + (MacXsBen(itype)%xs0sum(j)-MacXsBen(itype)%xs1sum(j))
    !    ENDDO
    !ENDIF
ENDDO
ELSEIF( benchxstype .EQ. 4 )THEN
    DO i = igb, ige
        XsMac%XsMacSm(:,i)=0
        DO j = MacXsBen(itype)%PnSM(1,i)%ib, MacXsBen(itype)%PnSM(1,i)%ie
            XsMac%XsMacSm(j,i) = MacXsBen(itype)%PnSM(1,i)%from(j)
        ENDDO
        IF( scatod .GE. 1)THEN
        XsMac%XsMacP1Sm(:,i) = 0
        DO j = MacXsBen(itype)%PnSM(2,i)%ib, MacXsBen(itype)%PnSM(2,i)%ie
          XsMac%XsMacP1Sm(j,i) = MacXsBen(itype)%PnSM(2,i)%from(j)
        ENDDO
        IF( scatod .GE. 2)THEN
          XsMac%XsMacP2Sm(:,i) = 0
          DO j = MacXsBen(itype)%PnSM(3,i)%ib, MacXsBen(itype)%PnSM(3,i)%ie
            XsMac%XsMacP2Sm(j,i) = MacXsBen(itype)%PnSM(3,i)%from(j)
          ENDDO
          IF( scatod .GE. 3 )THEN
            XsMac%XsMacP3Sm(:,i) = 0
            DO j = MacXsBen(itype)%PnSM(4,i)%ib, MacXsBen(itype)%PnSM(4,i)%ie
              XsMac%XsMacP3Sm(j,i) = MacXsBen(itype)%PnSM(4,i)%from(j)
            ENDDO
          ENDIF
        ENDIF
        ENDIF
    ENDDO
    DO i = igb, ige
        XsMac%xsmacs(i)=0
        DO j = igb, ige
            XsMac%xsmacs(i) = XsMac%xsmacs(i)+ XsMac%XsMacSm(i,j)
        ENDDO
        XsMac%xsmacstr(i)=XsMac%xsmacs(i) + MacXsBen(itype)%xstr(i)-MacXsBen(itype)%xst(i)
    ENDDO
    !IF( lscat1_ )THEN
    !    DO j = jsfr, jsto
    !        !XsMac%xsmacsm(j,j) = XsMac%xsmacsm(j,j) + (MacXsBen(itype)%xs0sum(j)-MacXsBen(itype)%xs1sum(j))
    !        XsMac%xsmacsm(j,j) = XsMac%xsmacsm(j,j) + MacXsBen(itype)%xs1sum(j)
    !        XsMac%xsmacs(j)    = XsMac%xsmacs(j)    + MacXsBen(itype)%xs1sum(j) ! sigs_tr > sigs_tot
    !    ENDDO
    !ENDIF
ENDIF

END SUBROUTINE

SUBROUTINE xsbaseBen_old(itype, igb, ige, jsfr, jsto, lscat1, XsMac)
!Base Macroscopic XS
IMPLICIT NONE
INTEGER :: i, j, k, ig
INTEGER :: itype, igb, ige, jsfr, jsto
LOGICAL :: lscat1
!INTEGER :: scatod
TYPE(XsMac_Type) :: XsMac

IF( scatod .GT. 0 )THEN
    lscat1=.TRUE.
ENDIF

XsMac%lFuel = MacXsBen(itype)%lFuel
DO i = igb, ige
  XsMac%xsmact(i) = MacXsBen(itype)%xst(i)
  XsMac%xsmaca(i) = MacXsBen(itype)%xsa(i)
  XsMac%xsmactr(i) = MacXsBen(itype)%xstr(i)
  XsMac%XsMacnf(i) = MacXsBen(itype)%xsnf(i)
  XsMac%XsMackf(i) = MacXsBen(itype)%xskf(i)
  XsMac%chi(i) = MacXsBen(itype)%chi(i)
  IF( benchxstype .eq. 1 )THEN
    DO j = jsfr, jsto
      XsMac%xsmacsm(j, i) = MacXsBen(itype)%xss(j, i)
    ENDDO
    IF( lscat1 )THEN
        XsMac%xsmacs(i)=0
        DO j = igb, ige
            XsMac%xsmacs(i) = XsMac%xsmacs(i)+ MacXsBen(itype)%xss(i, j)
        ENDDO
        XsMac%xsmacstr(i)=XsMac%xsmact(i)
        DO j = igb, ige
            XsMac%xsmacstr(i) = XsMac%xsmacstr(i)- MacXsBen(itype)%xss1(i, j)
        ENDDO
        DO j = jsfr, jsto
            XsMac%xsmacp1sm(j, i) = MacXsBen(itype)%xss1(j, i)
            IF( scatod .GE. 2)THEN
                XsMac%xsmacp2sm(j, i) = MacXsBen(itype)%xss2(j, i)
                IF( scatod .GE. 3 )THEN
                    XsMac%xsmacp3sm(j, i) = MacXsBen(itype)%xss3(j, i)
                ENDIF
            ENDIF
        ENDDO
    ENDIF
  ENDIF
ENDDO
IF( benchxstype .EQ. 4 )THEN
  DO i = igb, ige
      XsMac%XsMacSm(:,i)=0
      DO j = MacXsBen(itype)%PnSM(1,i)%ib, MacXsBen(itype)%PnSM(1,i)%ie
          XsMac%XsMacSm(j,i) = MacXsBen(itype)%PnSM(1,i)%from(j)
      ENDDO
      IF( lscat1 )THEN
          XsMac%XsMacP1Sm(:,i) = 0
        DO j = MacXsBen(itype)%PnSM(2,i)%ib, MacXsBen(itype)%PnSM(2,i)%ie
          XsMac%XsMacP1Sm(j,i) = MacXsBen(itype)%PnSM(2,i)%from(j)
        ENDDO
        IF( scatod .GE. 2)THEN
          XsMac%XsMacP2Sm(:,i) = 0
        DO j = MacXsBen(itype)%PnSM(3,i)%ib, MacXsBen(itype)%PnSM(3,i)%ie
          XsMac%XsMacP2Sm(j,i) = MacXsBen(itype)%PnSM(3,i)%from(j)
        ENDDO
        IF( scatod .GE. 3 )THEN
          XsMac%XsMacP3Sm(:,i) = 0
        DO j = MacXsBen(itype)%PnSM(4,i)%ib, MacXsBen(itype)%PnSM(4,i)%ie
          XsMac%XsMacP3Sm(j,i) = MacXsBen(itype)%PnSM(4,i)%from(j)
        ENDDO
        ENDIF
        ENDIF
      ENDIF
  ENDDO
  IF( lscat1 )THEN
  DO i = igb, ige
      XsMac%xsmacs(i)=0
      DO j = igb, ige
          XsMac%xsmacs(i) = XsMac%xsmacs(i)+ XsMac%XsMacSm(i,j)
      ENDDO
      XsMac%xsmacstr(i)=XsMac%xsmacs(i)
      DO j = igb, ige 
          XsMac%xsmacstr(i) = XsMac%xsmacstr(i)- XsMac%XsMacP1Sm(i, j)
      ENDDO
  ENDDO
  ENDIF
ENDIF

END SUBROUTINE
SUBROUTINE xstben(itype, igb, ige, xst)
!Transport XS
IMPLICIT NONE
INTEGER :: itype, igb, ige
REAL, POINTER :: xst(:)
INTEGER :: i, j, k, ig
DO ig = igb, ige, 1
  xst(ig) = MacXsBen(itype)%xst(ig)
ENDDO
END SUBROUTINE

SUBROUTINE xstrben(itype, igb, ige, xstr)
!Transport XS
IMPLICIT NONE
INTEGER :: itype, igb, ige
REAL, POINTER :: xstr(:)
INTEGER :: i, j, k, ig
DO ig = igb, ige, 1
  xstr(ig) = MacXsBen(itype)%xstr(ig)
ENDDO
END SUBROUTINE

FUNCTION GetXstrBen(itype, igb, ige)
IMPLICIT NONE
INTEGER :: itype, igb, ige, ig
REAL :: GetXstrBen(igb:ige)
INTEGER :: i, j, k
DO ig = igb, ige
  GetXsTrBen(ig) = MacXsBen(itype)%xstr(ig)
  !GetXsTrBen(ig) = MacXsBen(itype)%xst(ig)
ENDDO
!GetXstrBen1g = MacXsBen(itype)%xstr(ig)
END FUNCTION

FUNCTION GetXstBen(itype, igb, ige)
IMPLICIT NONE
INTEGER :: itype, igb, ige, ig
REAL :: GetXstBen(igb:ige)
INTEGER :: i, j, k
DO ig = igb, ige
  GetXsTBen(ig) = MacXsBen(itype)%xst(ig)
ENDDO
END FUNCTION

SUBROUTINE xsnfBen(itype, igb, ige, xsnf)
!Nu-FS scattering
IMPLICIT NONE
INTEGER :: itype, igb, ige
REAL, POINTER :: xsnf(:)
INTEGER :: i, j, k, ig
DO ig = igb, ige, 1
  xsnf(ig) = MacXsBen(itype)%xsnf(ig)
ENDDO
END SUBROUTINE

SUBROUTINE xskfBen(itype, igb, ige, xskf)
!In Scattering
IMPLICIT NONE
INTEGER :: itype, igb, ige
REAL, POINTER :: xskf(:)
INTEGER :: i, j, k, ig
DO ig = igb, ige
  xskf(ig) = MacXsBen(itype)%xskf(ig)
ENDDO
END SUBROUTINE

FUNCTION GetChiBen(itype, igb, ige)
IMPLICIT NONE
INTEGER :: itype, igb, ige
REAL :: GetChiBen(igb: ige)

INTEGER :: i, j, k, ig
DO ig = igb, ige
  GetChiBen(ig) = MacXsBen(itype)%chi(ig)
ENDDO

ENDFUNCTION

!SUBROUTINE xssben(itype, jgrp, jsfr1g, jsto1g, xss)
!IMPLICIT NONE
!INTEGER :: itype, jgrp, jsfr1g, jsto1g
!REAL, POINTER :: XSS(:,:)
!
!INTEGER :: i, j, k, ig
!DO j = jsfr1g, jsto1g
!  xss(j, jgrp) = MacXsBen(itype)%xss(j, jgrp)
!ENDDO
!ENDSUBROUTINE
!
!SUBROUTINE xssmben(itype, igb, ige, jsfr, jsto, xss)
!!Scattering Matrix
!IMPLICIT NONE
!INTEGER :: itype, igb, ige, jsfr, jsto
!REAL, POINTER :: XSS(:, :)
!
!INTEGER :: i, j, k, ig
!DO ig = igb, ige
!  DO j = jsfr, jsto
!    xss(j, ig) = MacXsBen(itype)%xss(j, ig)
!  ENDDO
!ENDDO
!END SUBROUTINE
!
!
!SUBROUTINE xssm1ben(itype, igb, ige, jsfr, jsto, xss)
!!Scattering Matrix
!IMPLICIT NONE
!INTEGER :: itype, igb, ige, jsfr, jsto
!REAL, POINTER :: XSS(:, :)
!
!INTEGER :: i, j, k, ig
!DO ig = igb, ige
!  DO j = jsfr, jsto
!    xss(j, ig) = MacXsBen(itype)%xss1(j, ig)
!  ENDDO
!ENDDO
!END SUBROUTINE
!SUBROUTINE xssm2ben(itype, igb, ige, jsfr, jsto, xss)
!!Scattering Matrix
!IMPLICIT NONE
!INTEGER :: itype, igb, ige, jsfr, jsto
!REAL, POINTER :: XSS(:, :)
!
!INTEGER :: i, j, k, ig
!DO ig = igb, ige
!  DO j = jsfr, jsto
!    xss(j, ig) = MacXsBen(itype)%xss2(j, ig)
!  ENDDO
!ENDDO
!END SUBROUTINE
!SUBROUTINE xssm3ben(itype, igb, ige, jsfr, jsto, xss)
!!Scattering Matrix
!IMPLICIT NONE
!INTEGER :: itype, igb, ige, jsfr, jsto
!REAL, POINTER :: XSS(:, :)
!
!INTEGER :: i, j, k, ig
!DO ig = igb, ige
!  DO j = jsfr, jsto
!    xss(j, ig) = MacXsBen(itype)%xss3(j, ig)
!  ENDDO
!ENDDO
!END SUBROUTINE


SUBROUTINE xssben(itype, jgrp, jsfr1g, jsto1g, xss, lscat1)
IMPLICIT NONE
INTEGER :: itype, jgrp, jsfr1g, jsto1g
REAL, POINTER :: XSS(:,:)
INTEGER :: i, j, k, ig
LOGICAL :: lscat1
SELECT CASE(benchxstype)
    case(1)
        DO j = jsfr1g, jsto1g
            xss(j, jgrp) = MacXsBen(itype)%xss(j, jgrp)
        ENDDO
        !DO j = jsfr1g, jsto1g
        !    IF( lscat1 ) xss(j, j) = xss(j, j) + MacXsBen(itype)%Xs1sum(j) 
        !ENDDO
    case(4)
        DO j = jsfr1g, jsto1g
            xss(j, jgrp) = 0
        ENDDO
        DO j =  MacXsBen(itype)%PnSM(1, jgrp)%ib, MacXsBen(itype)%PnSM(1, jgrp)%ie
            xss(j, jgrp) = MacXsBen(itype)%PnSM(1, jgrp)%from(j)
        ENDDO
        !DO j = jsfr1g, jsto1g
        !    IF( lscat1 ) xss(j, j) = xss(j, j) + MacXsBen(itype)%Xs1sum(j) 
        !ENDDO
ENDSELECT
    
        
ENDSUBROUTINE

SUBROUTINE xssmben(itype, igb, ige, jsfr, jsto, xss)
!Scattering Matrix
IMPLICIT NONE
INTEGER :: itype, igb, ige, jsfr, jsto
INTEGER :: i, j, k, ig
REAL, POINTER :: XSS(:, :)

SELECT CASE(benchxstype)
    case(1)
        DO ig = igb, ige
            DO j = jsfr, jsto
            xss(j, ig) = MacXsBen(itype)%xss(j, ig)
            ENDDO            
        ENDDO
    case(4)
        DO ig = igb, ige
            DO j = 1, igb, ige
                xss(j, ig) = 0
            ENDDO            
            DO j =  MacXsBen(itype)%PnSM(1, ig)%ib, MacXsBen(itype)%PnSM(1, ig)%ie
                xss(j, ig) = MacXsBen(itype)%PnSM(1, ig)%from(j)
            ENDDO
        ENDDO
ENDSELECT    
END SUBROUTINE


!
SUBROUTINE xssm1ben(itype, igb, ige, jsfr, jsto, xss)
!Scattering Matrix
IMPLICIT NONE
INTEGER :: itype, igb, ige, jsfr, jsto
REAL, POINTER :: XSS(:, :)

INTEGER :: i, j, k, ig

SELECT CASE(benchxstype)
    case(1)
        DO ig = igb, ige
            DO j = jsfr, jsto
            xss(j, ig) = MacXsBen(itype)%xss1(j, ig)
            ENDDO            
        ENDDO
    case(4)
        DO ig = igb, ige
            DO j = igb, ige
                xss(j, ig) = 0
            ENDDO
            DO j =  MacXsBen(itype)%PnSM(2, ig)%ib, MacXsBen(itype)%PnSM(2, ig)%ie
                xss(j, ig) = MacXsBen(itype)%PnSM(2, ig)%from(j)
            ENDDO
        ENDDO
ENDSELECT   
    
END SUBROUTINE
SUBROUTINE xssm2ben(itype, igb, ige, jsfr, jsto, xss)
!Scattering Matrix
IMPLICIT NONE
INTEGER :: itype, igb, ige, jsfr, jsto
REAL, POINTER :: XSS(:, :)

INTEGER :: i, j, k, ig

SELECT CASE(benchxstype)
    case(1)
        DO ig = igb, ige
            DO j = jsfr, jsto
            xss(j, ig) = MacXsBen(itype)%xss2(j, ig)
            ENDDO            
        ENDDO
    case(4)
        DO ig = igb, ige
            DO j = igb, ige
                xss(j, ig) = 0
            ENDDO
            DO j =  MacXsBen(itype)%PnSM(3, ig)%ib, MacXsBen(itype)%PnSM(3, ig)%ie
                xss(j, ig) = MacXsBen(itype)%PnSM(3, ig)%from(j)
            ENDDO
        ENDDO
ENDSELECT   
    
END SUBROUTINE
SUBROUTINE xssm3ben(itype, igb, ige, jsfr, jsto, xss)
!Scattering Matrix
IMPLICIT NONE
INTEGER :: itype, igb, ige, jsfr, jsto
REAL, POINTER :: XSS(:, :)

INTEGER :: i, j, k, ig


SELECT CASE(benchxstype)
    case(1)
        DO ig = igb, ige
            DO j = jsfr, jsto
            xss(j, ig) = MacXsBen(itype)%xss3(j, ig)
            ENDDO            
        ENDDO
    case(4)
        DO ig = igb, ige
            DO j = igb, ige
                xss(j, ig) = 0
            ENDDO
            DO j =  MacXsBen(itype)%PnSM(4, ig)%ib, MacXsBen(itype)%PnSM(4, ig)%ie
                xss(j, ig) = MacXsBen(itype)%PnSM(4, ig)%from(j)
            ENDDO
        ENDDO
ENDSELECT   
END SUBROUTINE



!Kinetic Parameter

SUBROUTINE DnpBetaBen(itype, beta)
IMPLICIT NONE
REAL :: Beta(DNP_NGRP)
INTEGER :: itype
INTEGER :: ig 

DO ig= 1, DNP_NGRP
  Beta(ig) = DNP_Beta(ig)
ENDDO
END SUBROUTINE

SUBROUTINE DnpLambdaBen(Lambda)
IMPLICIT NONE
REAL :: Lambda(DNP_NGRP)
INTEGER :: ig 

DO ig= 1, DNP_NGRP
  Lambda(ig) = DNP_Lambda(ig)
ENDDO
END SUBROUTINE

SUBROUTINE NeutVeloBen(itype, velo)
REAL :: Velo(ngben)
INTEGER :: itype
INTEGER :: ig
DO ig = 1, ngben
  Velo(ig) = Neut_Velo(ig)
ENDDO
END SUBROUTINE

SUBROUTINE ChidBen(Chid)
REAL :: Chid(ngben)
INTEGER :: ig
DO ig = 1, ngben
  Chid(ig) = Dneut_chi(ig)
ENDDO

END SUBROUTINE

END module 