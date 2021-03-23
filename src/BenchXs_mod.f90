module BenchXs
USE PARAM
USE TYPEDEF,     only : BenchXS_type, XsMac_Type, BenchKinParam_Type, DynBenchXs_Type, &
                        NeaCrpBenchXs_type
USE CNTL,        ONLY : nTracerCntl
IMPLICIT NONE
SAVE
INTEGER :: ngben, nxsltype, benchxstype, scatod
INTEGER :: nxschangeben, nxsltype_perturb
INTEGER :: nxsltype_total
!INTEGER, PRIVATE :: i, j, k, ig
TYPE(BenchXs_type), POINTER :: MacXsBen(:)
TYPE(BenchXs_type), PRIVATE, POINTER :: MacXsBen0(:)

TYPE(DynBenchXs_type), POINTER :: DynMacXsBen(:)

TYPE(NeaCrpBenchXs_type), POINTER :: NeaCrpXsBen(:)
TYPE(NeaCrpBenchXs_type), POINTER :: NeaCrpCAben(:)

LOGICAL :: lKinParam = .FALSE.
LOGICAL :: lTranInit = .FALSE. 
INTEGER :: DNP_NGRP
REAL, POINTER :: DNP_BETA(:)
REAL, POINTER :: DNP_Lambda(:)
REAL, POINTER :: Neut_Velo(:)
REAL, POINTER :: Dneut_Chi(:)

TYPE(BenchKinParam_Type), POINTER :: KinParamXsBen(:)

CONTAINS
SUBROUTINE XsBenChange_new(iso0, iso1, isonew, wt, lCusping, lCuspingDirection)
IMPLICIT NONE
INTEGER, INTENT(IN) :: iso0, iso1
INTEGER, INTENT(INOUT) :: isonew
REAL, INTENT(IN) :: wt
LOGICAL, INTENT(IN) :: lCusping, lCuspingDirection

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
  MacXsBen(isonew)%xst(ig) =  MacXsBen(iso0)%xst(ig) * wtbar + MacXsBen(iso1)%xst(ig) * wt
  MacXsBen(isonew)%xstr(ig) =  MacXsBen(iso0)%xstr(ig) * wtbar + MacXsBen(iso1)%xstr(ig) * wt
  MacXsBen(isonew)%xskf(ig) =  MacXsBen(iso0)%xskf(ig) * wtbar + MacXsBen(iso1)%xskf(ig) * wt
  MacXsBen(isonew)%xsnf(ig) =  MacXsBen(iso0)%xsnf(ig) * wtbar + MacXsBen(iso1)%xsnf(ig) * wt
  MacXsBen(isonew)%chi(ig) =  MacXsBen(iso0)%chi(ig) * wtbar + MacXsBen(iso1)%chi(ig) * wt
  MacXsBen(isonew)%xsa(ig) =  MacXsBen(iso0)%xsa(ig) * wtbar + MacXsBen(iso1)%xsa(ig) * wt
  DO ig0 = 1, ngben
    MacXsBen(isonew)%xss(ig0, ig) =  MacXsBen(iso0)%xss(ig0, ig) * wtbar + MacXsBen(iso1)%xss(ig0, ig) * wt
  ENDDO
ENDDO

IF(nTracerCntl%lKineticBen) THEN
  KinParamXsBen(isonew)%Beta = KinParamXsBen(iso0)%Beta * wtbar + KinParamXsBen(iso1)%Beta * wt
  KinParamXsBen(isonew)%Velo = KinParamXsBen(iso0)%Velo * wtbar + KinParamXsBen(iso1)%Velo * wt
ENDIF

IF(lCusping) THEN
  MacXsBen(isonew)%lCusping = .TRUE.
  IF(lCuspingDirection) THEN !Insertion
    MacXsBen(isonew)%wt = wt
    MacXsBen(isonew)%iso0 = iso0
    MacXsBen(isonew)%iso1 = iso1
  ELSE !Withdrawal
    MacXsBen(isonew)%wt = wtbar
    MacXsBen(isonew)%iso0 = iso1
    MacXsBen(isonew)%iso1 = iso0
  ENDIF
ENDIF

IF(nTracerCntl%libtyp .EQ. 11) THEN
  NeaCrpXsBen(isonew)%basexsl = iso0
END IF


END SUBROUTINE

SUBROUTINE XsBenNoise(XsNoise, nowtime)
USE TYPEDEF,          ONLY : XsNoise_Type
IMPLICIT NONE
TYPE(XsNoise_Type) :: XsNoise
REAL :: nowtime

REAL :: amp, freq, phase
REAL :: ampxs, nowphase, delxs, xsf_ratio
REAL :: sumXss
INTEGER :: itype, iso0, isonew
INTEGER :: ig, jg

itype = XsNoise%itype
iso0 = XsNoise%iso0

IF(.NOT. lTranInit) CALL InitXsChange()
IF(XsNoise%isonew .EQ. 0) THEN
  nxsltype_perturb = nxsltype_perturb + 1
  XsNoise%isonew = nxsltype_perturb
  CALL MakenCopyXsBen(XsNoise%isonew, iso0, .TRUE.)
END IF
isonew = XsNoise%isonew

amp = XsNoise%amp
freq = XsNoise%freq
phase = XsNoise%phase

nowphase = 2. * PI * freq * nowtime + phase

DO ig = 1, ngben
  SELECT CASE(itype)
  CASE(1)
    ampxs = amp * (MacXsBen(iso0)%xsa(ig) - MacXsBen(iso0)%xskf(ig))
    delxs = ampxs * sin(nowphase)
    MacXsBen(isonew)%xst(ig) = MacXsBen(iso0)%xst(ig) + delxs
    MacXsBen(isonew)%xstr(ig) = MacXsBen(iso0)%xstr(ig) + delxs
    MacXsBen(isonew)%xsa(ig) = MacXsBen(iso0)%xsa(ig) + delxs
  CASE(2)
    ampxs = amp * MacXsBen(iso0)%xskf(ig)
    delxs = ampxs * sin(nowphase)
    MacXsBen(isonew)%xst(ig) = MacXsBen(iso0)%xst(ig) + delxs
    MacXsBen(isonew)%xstr(ig) = MacXsBen(iso0)%xstr(ig) + delxs
    MacXsBen(isonew)%xsa(ig) = MacXsBen(iso0)%xsa(ig) + delxs
    MacXsBen(isonew)%xskf(ig) = MacXsBen(iso0)%xskf(ig) + delxs
    xsf_ratio = MacXsBen(isonew)%xskf(ig) / MacXsBen(iso0)%xskf(ig)
    MacXsBen(isonew)%xsnf(ig) = MacXsBen(iso0)%xsnf(ig) * xsf_ratio
  CASE(3)
    DO jg = 1, ngben
      ampxs = amp * MacXsBen(iso0)%xss(ig, jg)
      delxs = ampxs * sin(nowphase)
      MacXsBen(isonew)%xss(ig, jg) = MacXsBen(iso0)%xss(ig, jg) + delxs
    END DO 
    sumXss = sum(MacXsBen(isonew)%xss(ig,1:ngben))
    MacXsBen(isonew)%xst(ig) = MacXsBen(isonew)%xsa(ig) + sumXss
    MacXsBen(isonew)%xs0sum(ig) = sumXss
    MacXsBen(isonew)%xstr(ig) = MacXsBen(isonew)%xst(ig)

  END SELECT
END DO


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

IF(nTracerCntl%lKineticBen) THEN
  KinParamXsBen(iso0)%Beta(:) = KinParamXsBen(iso0)%Beta(:) * wtbar + KinParamXsBen(iso1)%Beta(:) * wt
  KinParamXsBen(iso0)%Velo(:) = KinParamXsBen(iso0)%Velo(:) * wtbar + KinParamXsBen(iso1)%Velo(:) * wt
ENDIF

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
  IF(nTracerCntl%lKineticBen) THEN
    CALL Dmalloc(KinParamXsBen(iso0)%Beta, DNP_nGRP)
    CALL Dmalloc(KinParamXsBen(iso0)%Velo, ngben)
  ENDIF
ENDIF

MacXsBen(iso0)%xst = MacXsBen(iso1)%xst;   MacXsBen(iso0)%xstr = MacXsBen(iso1)%xstr
MacXsBen(iso0)%xskf = MacXsBen(iso1)%xskf;   MacXsBen(iso0)%xsnf = MacXsBen(iso1)%xsnf
MacXsBen(iso0)%chi = MacXsBen(iso1)%chi;   MacXsBen(iso0)%xsa = MacXsBen(iso1)%xsa
MacXsBen(iso0)%xss = MacXsBen(iso1)%xss

IF(nTracerCntl%lKineticBen) THEN
  KinParamXsBen(iso0)%Beta = KinParamXsBen(iso1)%Beta
  KinParamXsBen(iso0)%Velo = KinParamXsBen(iso1)%Velo
ENDIF

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

SUBROUTINE xsbaseBen_cusping(itype, igb, ige, jsfr, jsto, lscat1, XsMac, phis, phisl, phisu, hz, hzl, hzu)
!Base Macroscopic XS
IMPLICIT NONE
INTEGER :: i, j, k, ig
INTEGER :: itype, igb, ige, jsfr, jsto
LOGICAL :: lscat1
LOGICAL :: lscat1_
!INTEGER :: scatod
TYPE(XsMac_Type) :: XsMac
REAL :: phis(:), phisl(:), phisu(:)
REAL :: hz, hzl, hzu

REAL :: wt, wtbar, wtg, wtgbar
REAL :: lflux, uflux
INTEGER :: iso0, iso1

lscat1_=lscat1
IF( scatod .GT. 0 )THEN
    lscat1_=.TRUE. .AND. lscat1_
ENDIF

iso0 = MacXsBen(itype)%iso0
iso1 = MacXsBen(itype)%iso1
wt = MacXsBen(itype)%wt
wtbar = 1._8 -wt

XsMac%lFuel = MacXsBen(itype)%lFuel
DO i = igb, ige
  uflux = (hzu * phisu(i) + wt * hz * phis(i)) / (hzu + wt * hz)
  lflux = (hzl * phisl(i) + wtbar * hz * phis(i)) / (hzl + wtbar * hz)
  wtg = wt * uflux / (wt * uflux + wtbar * lflux)
  wtgbar = 1._8 - wtg
  
  XsMac%xsmact(i) = MacXsBen(iso0)%xst(i) * wtgbar + MacXsBen(iso1)%xst(i) * wtg 
  XsMac%xsmaca(i) = MacXsBen(iso0)%xsa(i) * wtgbar + MacXsBen(iso1)%xsa(i) * wtg 
  XsMac%xsmactr(i) = MacXsBen(iso0)%xstr(i) * wtgbar + MacXsBen(iso1)%xstr(i) * wtg 
  XsMac%xsmacnf(i) = MacXsBen(iso0)%xsnf(i) * wtgbar + MacXsBen(iso1)%xsnf(i) * wtg 
  XsMac%xsmackf(i) = MacXsBen(iso0)%xskf(i) * wtgbar + MacXsBen(iso1)%xskf(i) * wtg 
  XsMac%chi(i) = MacXsBen(itype)%chi(i)
ENDDO

IF( benchxstype .eq. 1 )THEN
  DO i = igb, ige
    DO j = jsfr, jsto
      uflux = (hzu * phisu(j) + wt * hz * phis(j)) / (hzu + wt * hz)
      lflux = (hzl * phisl(j) + wtbar * hz * phis(j)) / (hzl + wtbar * hz)
      wtg = wt * uflux / (wt * uflux + wtbar * lflux)
      wtgbar = 1._8 - wtg

      XsMac%xsmacsm(j, i) = MacXsBen(iso0)%xss(j, i) * wtgbar + MacXsBen(iso1)%xss(j, i) * wtg
      IF (ScatOd .GE. 1) THEN
        XsMac%xsmacp1sm(j, i) = MacXsBen(iso0)%xss1(j, i) * wtgbar + MacXsBen(iso1)%xss1(j, i) * wtg
        IF (ScatOd .GE. 2) THEN
          XsMac%xsmacp2sm(j, i) = MacXsBen(iso0)%xss2(j, i) * wtgbar + MacXsBen(iso1)%xss2(j, i) * wtg
          IF (ScatOd .GE. 3) THEN
            XsMac%xsmacp3sm(j, i) = MacXsBen(iso0)%xss3(j, i) * wtgbar + MacXsBen(iso1)%xss3(j, i) * wtg
          ENDIF
        ENDIF
      ENDIF
    ENDDO
    XsMac%xsmacs(i)=0
    DO j = igb, ige
      XsMac%xsmacs(i) = XsMac%xsmacs(i)+ XsMac%xsmacsm(i, j)
    ENDDO
    XsMac%xsmacstr(i)=XsMac%xsmacs(i)
  ENDDO
ELSEIF( benchxstype .EQ. 4 )THEN
  PRINT*, 'AFW Decusping for benchmark XS type 4 is not supported'
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

FUNCTION GetXstrBen_Cusping(itype, igb, ige, phis, phisl, phisu, hz, hzl, hzu)
IMPLICIT NONE
INTEGER :: itype, igb, ige, ig
REAL :: GetXstrBen_Cusping(igb:ige)
REAL :: phis(igb:ige), phisl(igb:ige), phisu(igb:ige)
REAL :: hz, hzl, hzu
INTEGER :: i, j, k
INTEGER :: iso0, iso1
REAL :: wt, wtbar
REAL :: wtg, wtgbar
REAL :: lflux, uflux

iso0 = MacXsBen(itype)%iso0
iso1 = MacXsBen(itype)%iso1
wt = MacXsBen(itype)%wt
wtbar = 1._8 - wt

DO ig = igb, ige
  uflux = (hzu*phisu(ig) + wt*hz*phis(ig)) / (hzu + wt*hz)
  lflux = (hzl*phisl(ig) + wtbar*hz*phis(ig)) / (hzl + wtbar*hz)
  wtg = wt*uflux / (wt*uflux + wtbar*lflux)
  wtgbar = 1._8 - wtg
  GetXsTrBen_Cusping(ig) = MacXsBen(iso0)%xstr(ig) * wtgbar +  MacXsBen(iso1)%xstr(ig) * wtg
ENDDO
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

FUNCTION GetXstBen_Cusping(itype, igb, ige, phis, phisl, phisu, hz, hzl, hzu)
IMPLICIT NONE
INTEGER :: itype, igb, ige, ig
REAL :: GetXstBen_Cusping(igb:ige)
REAL :: phis(igb:ige), phisl(igb:ige), phisu(igb:ige)
REAL hz, hzl, hzu
INTEGER :: i, j, k
INTEGER :: iso0, iso1
REAL :: wt, wtbar
REAL :: wtg, wtgbar
REAL :: lflux, uflux

iso0 = MacXsBen(itype)%iso0
iso1 = MacXsBen(itype)%iso1
wt = MacXsBen(itype)%wt
wtbar = 1._8 - wt

DO ig = igb, ige
  uflux = (hzu*phisu(ig) + wt*hz*phis(ig)) / (hzu + wt*hz)
  lflux = (hzl*phisl(ig) + wtbar*hz*phis(ig)) / (hzl + wtbar*hz)
  wtg = wt*uflux / (wt*uflux + wtbar*lflux)
  wtgbar = 1._8 - wtg
  GetXsTBen_Cusping(ig) = MacXsBen(iso0)%xst(ig) * wtgbar +  MacXsBen(iso1)%xst(ig) * wtg
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

SUBROUTINE xssben_Cusping(itype, jgrp, jsfr1g, jsto1g, xss, lscat1, &
                          phis, phisl, phisu, hz, hzl, hzu)
IMPLICIT NONE
INTEGER :: itype, jgrp, jsfr1g, jsto1g
REAL, POINTER :: XSS(:,:)
LOGICAL :: lscat1
REAL :: phis(:), phisl(:), phisu(:)
REAL :: hz, hzl, hzu

REAL :: wt, wtbar, wtg, wtgbar
REAL :: lflux, uflux
INTEGER :: iso0, iso1
INTEGER :: i, j, k, ig

iso0 = MacXsBen(itype)%iso0
iso1 = MacXsBen(itype)%iso1
wt = MacXsBen(itype)%wt
wtbar = 1._8 - wt

SELECT CASE(benchxstype)
    case(1)
      DO j = jsfr1g, jsto1g
        uflux = (hzu*phisu(j) + wt*hz*phis(j)) / (hzu + wt*hz)
        lflux = (hzl*phisl(j) + wtbar*hz*phis(j)) / (hzl + wtbar*hz)
        wtg = wt*uflux / (wt*uflux + wtbar*lflux)
        wtgbar = 1._8 - wtg

        xss(j, jgrp) = MacXsBen(iso0)%xss(j, jgrp) * wtgbar + MacXsBen(iso1)%xss(j, jgrp) * wtg
      ENDDO
      !DO j = jsfr1g, jsto1g
      !    IF( lscat1 ) xss(j, j) = xss(j, j) + MacXsBen(itype)%Xs1sum(j)
      !ENDDO
    case(4)
      PRINT*, 'AFW Decusping for benchmark XS type 4 is not supported'
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

IF(nTracerCntl%lKineticBen) THEN
  Beta = KinParamXsBen(itype)%Beta
ELSE
  Beta = DNP_Beta
END IF

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

IF(nTracerCntl%lKineticBen) THEN
  Velo = KinParamXsBen(itype)%Velo
ELSE
  Velo = Neut_Velo
END IF

END SUBROUTINE

SUBROUTINE ChidBen(Chid)
REAL :: Chid(ngben)
INTEGER :: ig
DO ig = 1, ngben
  Chid(ig) = Dneut_chi(ig)
ENDDO

END SUBROUTINE

SUBROUTINE alloc_DynBenchXS(ixsl, ntemp, nprec)
USE PARAM
USE allocs
IMPLICIT NONE
INTEGER :: ixsl, ntemp, nprec

DynMacXsBen(ixsl)%lFuel = .FALSE.

DynMacXsBen(ixsl)%ntemp = ntemp
ALLOCATE(DynMacXsBen(ixsl)%temp(ntemp))

ALLOCATE(DynMacXsBen(ixsl)%xst(ntemp, ngben))
ALLOCATE(DynMacXsBen(ixsl)%xstr(ntemp, ngben))
ALLOCATE(DynMacXsBen(ixsl)%xskf(ntemp, ngben))
ALLOCATE(DynMacXsBen(ixsl)%xsnf(ntemp, ngben))
ALLOCATE(DynMacXsBen(ixsl)%chi(ntemp, ngben))
ALLOCATE(DynMacXsBen(ixsl)%xsa(ntemp, ngben))
ALLOCATE(DynMacXsBen(ixsl)%xs0sum(ntemp, ngben))
ALLOCATE(DynMacXsBen(ixsl)%xs1sum(ntemp, ngben))

ALLOCATE(DynMacXsBen(ixsl)%xss(ntemp, ngben, ngben))
DynMacXsBen(ixsl)%xss0 => DynMacXsBen(ixsl)%xss
IF(scatod .GE. 1) THEN
  ALLOCATE(DynMacXsBen(ixsl)%xss1(ntemp, ngben, ngben))
  IF(scatod .GE. 2) THEN
    ALLOCATE(DynMacXsBen(ixsl)%xss2(ntemp, ngben, ngben))
    IF(scatod .GE. 3) THEN
      ALLOCATE(DynMacXsBen(ixsl)%xss3(ntemp, ngben, ngben))
    END IF
  END IF
END IF

IF(lKinParam) THEN
  ALLOCATE(DynMacXsBen(ixsl)%Beta(ntemp, nprec))
  ALLOCATE(DynMacXsBen(ixsl)%Lambda(ntemp, nprec))
  ALLOCATE(DynMacXsBen(ixsl)%Velo(ntemp, ngben))
  ALLOCATE(DynMacXsBen(ixsl)%Chid(ntemp, ngben))
END IF

END SUBROUTINE

SUBROUTINE tempInterpolation(itype, temp, wt1, temp1, temp2)
IMPLICIT NONE
INTEGER :: itype
REAL :: temp, wt1
INTEGER :: temp1, temp2

INTEGER :: ntemp
INTEGER :: itemp

ntemp = DynMacXsBen(itype)%ntemp


IF(temp .LE. DynMacXsBen(itype)%temp(1)) THEN
  temp1 = 1
  temp2 = 2
  wt1 = (DynMacXsBen(itype)%temp(temp2) - temp) / (DynMacXsBen(itype)%temp(temp2) - DynMacXsBen(itype)%temp(temp1))
  RETURN
END IF
DO itemp = 2, ntemp
  IF(temp .LE. DynMacXsBen(itype)%temp(itemp)) THEN
    temp2 = itemp
    temp1 = itemp-1
    wt1 = (DynMacXsBen(itype)%temp(temp2) - temp) / (DynMacXsBen(itype)%temp(temp2) - DynMacXsBen(itype)%temp(temp1))
    RETURN
  END IF
END DO 

temp2 = ntemp
temp1 = ntemp - 1
wt1 = (DynMacXsBen(itype)%temp(temp2) - temp) / (DynMacXsBen(itype)%temp(temp2) - DynMacXsBen(itype)%temp(temp1))

END SUBROUTINE

SUBROUTINE xsDynBenChange(iso0, iso1, isonew, wt, lCusping, lCuspingDirection)
IMPLICIT NONE
INTEGER, INTENT(IN) :: iso0, iso1
INTEGER, INTENT(INOUT) :: isonew
REAL, INTENT(IN) :: wt
LOGICAL, INTENT(IN) :: lCusping, lCuspingDirection

REAL :: wtbar
REAL :: twt, twtbar
REAL :: temp0
INTEGER :: ig, jg, itemp, jtemp1, jtemp2

IF(.NOT. lTranInit) THEN
  nxsltype_perturb = nxsltype
  lTranInit = .TRUE.
END IF
IF(isonew .EQ. 0) THEN
  nxsltype_perturb = nxsltype_perturb + 1
  isonew = nxsltype_perturb
  CALL  alloc_DynBenchXS(isonew, DynMacXsBen(iso0)%ntemp, DNP_nGRP)
  DynMacXsBen(isonew)%temp = DynMacXsBen(iso0)%temp
END IF

wtbar = 1. - wt
DO itemp = 1, DynMacXsBen(iso0)%ntemp
  temp0 = DynMacXsBen(iso0)%temp(itemp)
  CALL tempInterpolation(iso1, temp0, twt, jtemp1, jtemp2)
  twtbar = wt * (1. - twt)
  twt = twt * wt
  DO ig = 1, ngben
    DynMacXsBen(isonew)%xst(itemp, ig) = DynMacXsBen(iso0)%xst(itemp, ig) * wtbar + DynMacXsBen(iso1)%xst(jtemp1, ig) * twt + DynMacXsBen(iso1)%xst(jtemp2, ig) * twtbar
    DynMacXsBen(isonew)%xstr(itemp, ig) = DynMacXsBen(iso0)%xstr(itemp, ig) * wtbar + DynMacXsBen(iso1)%xstr(jtemp1, ig) * twt + DynMacXsBen(iso1)%xstr(jtemp2, ig) * twtbar
    DynMacXsBen(isonew)%xskf(itemp, ig) = DynMacXsBen(iso0)%xskf(itemp, ig) * wtbar + DynMacXsBen(iso1)%xskf(jtemp1, ig) * twt + DynMacXsBen(iso1)%xskf(jtemp2, ig) * twtbar
    DynMacXsBen(isonew)%xsnf(itemp, ig) = DynMacXsBen(iso0)%xsnf(itemp, ig) * wtbar + DynMacXsBen(iso1)%xsnf(jtemp1, ig) * twt + DynMacXsBen(iso1)%xsnf(jtemp2, ig) * twtbar
    DynMacXsBen(isonew)%chi(itemp, ig) = DynMacXsBen(iso0)%chi(itemp, ig) * wtbar + DynMacXsBen(iso1)%chi(jtemp1, ig) * twt + DynMacXsBen(iso1)%chi(jtemp2, ig) * twtbar
    DynMacXsBen(isonew)%xsa(itemp, ig) = DynMacXsBen(iso0)%xsa(itemp, ig) * wtbar + DynMacXsBen(iso1)%xsa(jtemp1, ig) * twt + DynMacXsBen(iso1)%xsa(jtemp2, ig) * twtbar
    DO jg = 1, ngben
      DynMacXsBen(isonew)%xss(itemp, jg, ig) = DynMacXsBen(iso0)%xss(itemp, jg, ig) * wtbar + &
        DynMacXsBen(iso1)%xss(jtemp1, jg, ig) * twt + DynMacXsBen(iso1)%xss(jtemp2, jg, ig) * twtbar
    END DO
    IF(scatod .GE. 1) THEN
      DO jg = 1, ngben
        DynMacXsBen(isonew)%xss1(itemp, jg, ig) = DynMacXsBen(iso0)%xss1(itemp, jg, ig) * wtbar + &
          DynMacXsBen(iso1)%xss1(jtemp1, jg, ig) * twt + DynMacXsBen(iso1)%xss1(jtemp2, jg, ig) * twtbar
      END DO
      IF(scatod .GE. 2) THEN
        DO jg = 1, ngben
          DynMacXsBen(isonew)%xss2(itemp, jg, ig) = DynMacXsBen(iso0)%xss2(itemp, jg, ig) * wtbar + &
            DynMacXsBen(iso1)%xss2(jtemp1, jg, ig) * twt + DynMacXsBen(iso1)%xss2(jtemp2, jg, ig) * twtbar
        END DO
        IF(scatod .GE. 3) THEN
          DO jg = 1, ngben
            DynMacXsBen(isonew)%xss3(itemp, jg, ig) = DynMacXsBen(iso0)%xss3(itemp, jg, ig) * wtbar + &
              DynMacXsBen(iso1)%xss3(jtemp1, jg, ig) * twt + DynMacXsBen(iso1)%xss3(jtemp2, jg, ig) * twtbar
          END DO
        END IF
      END IF
    END IF
    IF(lKinParam) THEN
      DynMacXsBen(isonew)%velo(itemp, ig) = DynMacXsBen(iso0)%velo(itemp, ig) * wtbar + DynMacXsBen(iso1)%velo(jtemp1, ig) * twt + DynMacXsBen(iso1)%velo(jtemp2, ig) * twtbar
    END IF
  END DO
  IF(lKinParam) THEN
    DO ig = 1, DNP_nGRP
      DynMacXsBen(isonew)%beta(itemp, ig) = DynMacXsBen(iso0)%beta(itemp, ig) * wtbar + DynMacXsBen(iso1)%beta(jtemp1, ig) * twt + DynMacXsBen(iso1)%beta(jtemp2, ig) * twtbar
    END DO
  END IF
END DO 

IF(lCusping) THEN
  DynMacXsBen(isonew)%lCusping = .TRUE.
  IF(lCuspingDirection) THEN !Insertion
    DynMacXsBen(isonew)%wt = wt
    DynMacXsBen(isonew)%iso0 = iso0
    DynMacXsBen(isonew)%iso1 = iso1
  ELSE !Withdrawal
    DynMacXsBen(isonew)%wt = wtbar
    DynMacXsBen(isonew)%iso0 = iso1
    DynMacXsBen(isonew)%iso1 = iso0
  END IF
END IF

END SUBROUTINE

SUBROUTINE xsbaseDynBen(itype, temp, igb, ige, jsfr, jsto, lscat1, XsMac)
IMPLICIT NONE
INTEGER :: itype
REAL :: temp
INTEGER :: igb, ige, jsfr, jsto
LOGICAL :: lscat1
TYPE(XsMac_Type) :: XsMac

REAL :: wt, wtbar
INTEGER :: temp1, temp2
INTEGER :: ig, jg, itemp
LOGICAL :: lscat1_

lscat1_ = lscat1
IF(scatod .GT. 0) THEN
  lscat1_ = .TRUE. .AND. lscat1_
END IF

CALL tempInterpolation(itype, temp, wt, temp1, temp2)
wtbar = 1. - wt

XsMac%lFuel = DynMacXsBen(itype)%lFuel
DO ig = igb, ige
  XsMac%xsmact(ig) = wt * DynMacXsBen(itype)%xst(temp1, ig)  + wtbar * DynMacXsBen(itype)%xst(temp2, ig)
  XsMac%xsmaca(ig) = wt * DynMacXsBen(itype)%xsa(temp1, ig)  + wtbar * DynMacXsBen(itype)%xsa(temp2, ig)
  XsMac%xsmactr(ig) = wt * DynMacXsBen(itype)%xstr(temp1, ig)  + wtbar * DynMacXsBen(itype)%xstr(temp2, ig)
  XsMac%xsmacnf(ig) = wt * DynMacXsBen(itype)%xsnf(temp1, ig)  + wtbar * DynMacXsBen(itype)%xsnf(temp2, ig)
  XsMac%xsmackf(ig) = wt * DynMacXsBen(itype)%xskf(temp1, ig)  + wtbar * DynMacXsBen(itype)%xskf(temp2, ig)
  XsMac%chi(ig) = wt * DynMacXsBen(itype)%chi(temp1, ig)  + wtbar * DynMacXsBen(itype)%chi(temp2, ig)
END DO 

DO ig = igb, ige
  DO jg = jsfr, jsto
    XsMac%xsmacsm(jg, ig) = wt * DynMacXsBen(itype)%xss(temp1, jg, ig) + wtbar * DynMacXsBen(itype)%xss(temp2, jg, ig)
  END DO
  XsMac%xsmacs(ig) = 0.
  DO jg = igb, ige
    XsMac%xsmacs(ig) = XsMac%xsmacs(ig) + wt * DynMacXSBen(itype)%xss(temp1, ig, jg) + wtbar * DynMacXsBen(itype)%xss(temp2, ig, jg)
  END DO
  XsMac%xsmacstr(ig) = XsMac%xsmacs(ig)
  DO jg = jsfr, jsto
    IF(scatod .GE. 1) THEN
      XsMac%xsmacp1sm(jg, ig) = wt * DynMacXsBen(itype)%xss1(temp1, jg, ig) + wtbar * DynMacXsBen(itype)%xss1(temp2, jg, ig)
      IF(scatod .GE. 2) THEN
        XsMac%xsmacp2sm(jg, ig) = wt * DynMacXsBen(itype)%xss2(temp1, jg, ig) + wtbar * DynMacXsBen(itype)%xss2(temp2, jg, ig)
        IF(scatod .GE. 3) THEN
          XsMac%xsmacp3sm(jg, ig) = wt * DynMacXsBen(itype)%xss3(temp1, jg, ig) + wtbar * DynMacXsBen(itype)%xss3(temp2, jg, ig)
        END IF
      END IF
    END IF
  END DO
END DO 

END SUBROUTINE 

SUBROUTINE xsbaseDynBen_cusping(itype, temp, igb, ige, jsfr, jsto, lscat1, XsMac, phis, phisl, phisu, hz, hzl, hzu) 
IMPLICIT NONE
INTEGER :: itype
REAL :: temp
INTEGER :: igb, ige, jsfr, jsto
LOGICAL :: lscat1
TYPE(XsMac_Type) :: XsMac
REAL :: phis(:), phisl(:), phisu(:)
REAL :: hz, hzl, hzu

REAL :: wt, wtbar, wtg, wtgbar
REAL :: twt0, twtbar0, twt1, twtbar1
REAL :: twtg0, twtgbar0, twtg1, twtgbar1
REAL :: lflux, uflux
INTEGER :: iso0, iso1, itemp1, itemp2, jtemp1, jtemp2
INTEGER :: ig, jg


iso0 = DynMacXsBen(itype)%iso0
iso1 = DynMacXsBen(itype)%iso1
wt = DynMacXsBen(itype)%wt
wtbar = 1. - wt

CALL tempInterpolation(iso0, temp, twt0, itemp1, itemp2)
twtbar0 = 1. - twt0
CALL tempInterpolation(iso1, temp, twt1, jtemp1, jtemp2)
twtbar1 = 1. - twt1

XsMac%lFuel = DynMacXsBen(itype)%lFuel
DO ig = igb, ige
  uflux = (hzu * phisu(ig) + wt * hz * phis(ig)) / (hzu + wt * hz)
  lflux = (hzl * phisl(ig) + wtbar * hz * phis(ig)) / (hzl + wtbar * hz)
  wtg = wt * uflux / (wt * uflux + wtbar * lflux)
  wtgbar = 1. - wtg

  twtg0 = wtgbar * twt0
  twtgbar0 = wtgbar * twtbar0
  twtg1 = wtg * twt1
  twtgbar1 = wtg * twtbar1

  XsMac%xsmact(ig) = DynMacXsBen(iso0)%xst(itemp1, ig) * twtg0 + DynMacXsBen(iso0)%xst(itemp2, ig) * twtgbar0 + &
                     DynMacXsBen(iso1)%xst(jtemp1, ig) * twtg1 + DynMacXsBen(iso1)%xst(jtemp2, ig) * twtgbar1
  XsMac%xsmaca(ig) = DynMacXsBen(iso0)%xsa(itemp1, ig) * twtg0 + DynMacXsBen(iso0)%xsa(itemp2, ig) * twtgbar0 + &
                     DynMacXsBen(iso1)%xsa(jtemp1, ig) * twtg1 + DynMacXsBen(iso1)%xsa(jtemp2, ig) * twtgbar1
  XsMac%xsmactr(ig) = DynMacXsBen(iso0)%xstr(itemp1, ig) * twtg0 + DynMacXsBen(iso0)%xstr(itemp2, ig) * twtgbar0 + &
                      DynMacXsBen(iso1)%xstr(jtemp1, ig) * twtg1 + DynMacXsBen(iso1)%xstr(jtemp2, ig) * twtgbar1
  XsMac%xsmacnf(ig) = DynMacXsBen(iso0)%xsnf(itemp1, ig) * twtg0 + DynMacXsBen(iso0)%xsnf(itemp2, ig) * twtgbar0 + &
                      DynMacXsBen(iso1)%xsnf(jtemp1, ig) * twtg1 + DynMacXsBen(iso1)%xsnf(jtemp2, ig) * twtgbar1
  XsMac%xsmackf(ig) = DynMacXsBen(iso0)%xskf(itemp1, ig) * twtg0 + DynMacXsBen(iso0)%xskf(itemp2, ig) * twtgbar0 + &
                      DynMacXsBen(iso1)%xskf(jtemp1, ig) * twtg1 + DynMacXsBen(iso1)%xskf(jtemp2, ig) * twtgbar1
  XsMac%chi(ig) = DynMacXsBen(iso0)%chi(itemp1, ig) * twtg0 + DynMacXsBen(iso0)%chi(itemp2, ig) * twtgbar0 + &
                  DynMacXsBen(iso1)%chi(jtemp1, ig) * twtg1 + DynMacXsBen(iso1)%chi(jtemp2, ig) * twtgbar1
  
END DO

DO ig = igb, ige
  DO jg = jsfr, jsto
    uflux = (hzu * phisu(jg) + wt * hz * phis(jg)) / (hzu + wt * hz)
    lflux = (hzl * phisl(jg) + wtbar * hz * phis(jg)) / (hzl + wtbar * hz)
    wtg = wt * uflux / (wt * uflux + wtbar * lflux)
    wtgbar = 1. - wtg

    twtg0 = wtgbar * twt0
    twtgbar0 = wtgbar * twtbar0
    twtg1 = wtg * twt1
    twtgbar1 = wtg * twtbar1

    
    XsMac%xsmacsm(jg, ig) = DynMacXsBen(iso0)%xss(itemp1, jg, ig) * twtg0 + DynMacXsBen(iso0)%xss(itemp2, jg, ig) * twtgbar0 + &
                            DynMacXsBen(iso1)%xss(jtemp1, jg, ig) * twtg1 + DynMacXsBen(iso1)%xss(jtemp2, jg, ig) * twtgbar1
    IF(scatod .GE. 1) THEN
      XsMac%xsmacp1sm(jg, ig) = DynMacXsBen(iso0)%xss1(itemp1, jg, ig) * twtg0 + DynMacXsBen(iso0)%xss1(itemp2, jg, ig) * twtgbar0 + &
                                DynMacXsBen(iso1)%xss1(jtemp1, jg, ig) * twtg1 + DynMacXsBen(iso1)%xss1(jtemp2, jg, ig) * twtgbar1
      IF(scatod .GE. 2) THEN
        XsMac%xsmacp2sm(jg, ig) = DynMacXsBen(iso0)%xss2(itemp1, jg, ig) * twtg0 + DynMacXsBen(iso0)%xss2(itemp2, jg, ig) * twtgbar0 + &
                                  DynMacXsBen(iso1)%xss2(jtemp1, jg, ig) * twtg1 + DynMacXsBen(iso1)%xss2(jtemp2, jg, ig) * twtgbar1
        IF(scatod .GE. 3) THEN
          XsMac%xsmacp3sm(jg, ig) = DynMacXsBen(iso0)%xss3(itemp1, jg, ig) * twtg0 + DynMacXsBen(iso0)%xss3(itemp2, jg, ig) * twtgbar0 + &
                                    DynMacXsBen(iso1)%xss3(jtemp1, jg, ig) * twtg1 + DynMacXsBen(iso1)%xss3(jtemp2, jg, ig) * twtgbar1
        END IF
      END IF
    END IF
  END DO
  XsMac%xsmacs(ig) = 0.
  DO jg = igb, ige
    XsMac%xsmacs(ig) = XsMac%xsmacs(ig) + XsMac%xsmacsm(ig, jg)
  END DO
  XsMac%xsmacstr(ig) = XsMac%xsmacs(ig)
END DO

END SUBROUTINE

FUNCTION GetXstrDynBen(itype, temp, igb, ige)
IMPLICIT NONE
REAL :: GetXstrDynBen(igb:ige)
INTEGER :: itype
REAL :: temp
INTEGER :: igb, ige

REAL :: wt, wtbar
INTEGER :: temp1, temp2
INTEGER :: ig

CALL tempInterpolation(itype, temp, wt, temp1, temp2)
wtbar = 1. - wt
DO ig = igb, ige
  GetXsTrDynBen(ig) = DynMacXsBen(itype)%xstr(temp1, ig) * wt + DynMacXsBen(itype)%xstr(temp2, ig) * wtbar
END DO 

END FUNCTION 

FUNCTION GetXstrDynBen_cusping(itype, temp, igb, ige, phis, phisl, phisu, hz, hzl, hzu)
IMPLICIT NONE
REAL :: GetXsTrDynBen_cusping(igb:ige)
INTEGER :: itype
REAL :: temp
INTEGER :: igb, ige
REAL :: phis(igb:ige), phisl(igb:ige), phisu(igb:ige)
REAL :: hz, hzl, hzu

REAL :: wt, wtbar, wtg, wtgbar
REAL :: twt0, twtbar0, twt1, twtbar1
REAL :: twtg0, twtgbar0, twtg1, twtgbar1
REAL :: lflux, uflux
INTEGER :: iso0, iso1, itemp1, itemp2, jtemp1, jtemp2
INTEGER :: ig

iso0 = DynMacXsBen(itype)%iso0
iso1 = DynMacXsBen(itype)%iso1
wt = DynMacXsBen(itype)%wt
wtbar = 1. - wt

CALL tempInterpolation(iso0, temp, twt0, itemp1, itemp2)
twtbar0 = 1. - twt0
CALL tempInterpolation(iso1, temp, twt1, jtemp1, jtemp2)
twtbar1 = 1. - twt1

DO ig = igb, ige
  uflux = (hzu * phisu(ig) + wt * hz * phis(ig)) / (hzu + wt * hz)
  lflux = (hzl * phisl(ig) + wtbar * hz * phis(ig)) / (hzl + wtbar * hz)
  wtg = wt * uflux / (wt * uflux + wtbar * lflux)
  wtgbar = 1. - wtg

  twtg0 = wtgbar * twt0
  twtgbar0 = wtgbar * twtbar0
  twtg1 = wtg * twt1
  twtgbar1 = wtg * twtbar1

  GetXsTrDynBen_cusping(ig) = DynMacXsBen(iso0)%xstr(itemp1, ig) * twtg0 + DynMacXsBen(iso0)%xstr(itemp2, ig) * twtgbar0 + &
                              DynMacXsBen(iso1)%xstr(jtemp1, ig) * twtg1 + DynMacXsBen(iso1)%xstr(jtemp2, ig) * twtgbar1
END DO
END FUNCTION

FUNCTION GetXstDynBen(itype, temp, igb, ige)
IMPLICIT NONE
REAL :: GetXstDynBen(igb:ige)
INTEGER :: itype
REAL :: temp
INTEGER :: igb, ige

REAL :: wt, wtbar
INTEGER :: temp1, temp2
INTEGER :: ig

CALL tempInterpolation(itype, temp, wt, temp1, temp2)
wtbar = 1. - wt
DO ig = igb, ige
  GetXsTDynBen(ig) = DynMacXsBen(itype)%xst(temp1, ig) * wt + DynMacXsBen(itype)%xst(temp2, ig) * wtbar
END DO 

END FUNCTION 

FUNCTION GetXstDynBen_cusping(itype, temp, igb, ige, phis, phisl, phisu, hz, hzl, hzu)
IMPLICIT NONE
REAL :: GetXsTDynBen_cusping(igb:ige)
INTEGER :: itype
REAL :: temp
INTEGER :: igb, ige
REAL :: phis(igb:ige), phisl(igb:ige), phisu(igb:ige)
REAL :: hz, hzl, hzu

REAL :: wt, wtbar, wtg, wtgbar
REAL :: twt0, twtbar0, twt1, twtbar1
REAL :: twtg0, twtgbar0, twtg1, twtgbar1
REAL :: lflux, uflux
INTEGER :: iso0, iso1, itemp1, itemp2, jtemp1, jtemp2
INTEGER :: ig

iso0 = DynMacXsBen(itype)%iso0
iso1 = DynMacXsBen(itype)%iso1
wt = DynMacXsBen(itype)%wt
wtbar = 1. - wt

CALL tempInterpolation(iso0, temp, twt0, itemp1, itemp2)
twtbar0 = 1. - twt0
CALL tempInterpolation(iso1, temp, twt1, jtemp1, jtemp2)
twtbar1 = 1. - twt1

DO ig = igb, ige
  uflux = (hzu * phisu(ig) + wt * hz * phis(ig)) / (hzu + wt * hz)
  lflux = (hzl * phisl(ig) + wtbar * hz * phis(ig)) / (hzl + wtbar * hz)
  wtg = wt * uflux / (wt * uflux + wtbar * lflux)
  wtgbar = 1. - wtg

  twtg0 = wtgbar * twt0
  twtgbar0 = wtgbar * twtbar0
  twtg1 = wtg * twt1
  twtgbar1 = wtg * twtbar1

  GetXsTDynBen_cusping(ig) = DynMacXsBen(iso0)%xst(itemp1, ig) * twtg0 + DynMacXsBen(iso0)%xst(itemp2, ig) * twtgbar0 + &
                             DynMacXsBen(iso1)%xst(jtemp1, ig) * twtg1 + DynMacXsBen(iso1)%xst(jtemp2, ig) * twtgbar1
END DO
END FUNCTION

SUBROUTINE xsnfDynBen(itype, temp, igb, ige, xsnf)
IMPLICIT NONE
INTEGER :: itype 
REAL :: temp
INTEGER :: igb, ige
REAL, POINTER :: xsnf(:)

REAL :: wt, wtbar
INTEGER :: temp1, temp2
INTEGER :: ig

CALL tempInterpolation(itype, temp, wt, temp1, temp2)
wtbar = 1. - wt
DO ig = igb, ige
  xsnf(ig) = DynMacXsBen(itype)%xsnf(temp1, ig) * wt + DynMacXsBen(itype)%xsnf(temp2, ig) * wtbar
END DO 
END SUBROUTINE

SUBROUTINE xskfDynBen(itype, temp, igb, ige, xskf)
IMPLICIT NONE
INTEGER :: itype 
REAL :: temp
INTEGER :: igb, ige
REAL, POINTER :: xskf(:)

REAL :: wt, wtbar
INTEGER :: temp1, temp2
INTEGER :: ig

CALL tempInterpolation(itype, temp, wt, temp1, temp2)
wtbar = 1. - wt
DO ig = igb, ige
  xskf(ig) = DynMacXsBen(itype)%xskf(temp1, ig) * wt + DynMacXsBen(itype)%xskf(temp2, ig) * wtbar
END DO 
END SUBROUTINE

FUNCTION GetChiDynBen(itype, temp, igb, ige)
IMPLICIT NONE
REAL :: GetChiDynBen(igb:ige)
INTEGER :: itype
REAL :: temp
INTEGER :: igb, ige

REAL :: wt, wtbar
INTEGER :: temp1, temp2
INTEGER :: ig

CALL tempInterpolation(itype, temp, wt, temp1, temp2)
wtbar = 1. - wt
DO ig = igb, ige
  GetChiDynBen(ig) = DynMacXsBen(itype)%Chi(temp1, ig) * wt + DynMacXsBen(itype)%Chi(temp2, ig) * wtbar
END DO 

END FUNCTION 

SUBROUTINE xssDynBen(itype, temp, jgrp, jsfr1g, jsto1g, xss, lscat1)
IMPLICIT NONE
INTEGER :: itype
REAL :: temp
INTEGER:: jgrp, jsfr1g, jsto1g
REAL, POINTER :: xss(:, :)
LOGICAL :: lscat1

REAL :: wt, wtbar
INTEGER :: temp1, temp2
INTEGER :: ig

CALL tempInterpolation(itype, temp, wt, temp1, temp2)
wtbar = 1. - wt

DO ig = jsfr1g, jsto1g
  xss(ig, jgrp) = DynMacXsBen(itype)%xss(temp1, ig, jgrp) * wt + DynMacXsBen(itype)%xss(temp2, ig, jgrp) * wtbar
END DO

END SUBROUTINE

SUBROUTINE xssDynBen_cusping(itype, temp, jgrp, jsfr1g, jsto1g, xss, lscat1, &
                             phis, phisl, phisu, hz, hzl, hzu)
IMPLICIT NONE
INTEGER :: itype
REAL :: temp
INTEGER :: jgrp, jsfr1g, jsto1g
REAL, POINTER :: xss(:, :)
LOGICAL :: lscat1
REAL :: phis(:), phisl(:), phisu(:)
REAL :: hz, hzl, hzu

REAL :: wt, wtbar, wtg, wtgbar
REAL :: twt0, twtbar0, twt1, twtbar1
REAL :: twtg0, twtgbar0, twtg1, twtgbar1
REAL :: lflux, uflux
INTEGER :: iso0, iso1, itemp1, itemp2, jtemp1, jtemp2
INTEGER :: ig

iso0 = DynMacXsBen(itype)%iso0
iso1 = DynMacXsBen(itype)%iso1
wt = DynMacXsBen(itype)%wt
wtbar = 1. - wt

CALL tempInterpolation(iso0, temp, twt0, itemp1, itemp2)
twtbar0 = 1. - twt0
CALL tempInterpolation(iso1, temp, twt1, jtemp1, jtemp2)
twtbar1 = 1. - twt1

DO ig = jsfr1g, jsto1g
  uflux = (hzu * phisu(ig) + wt * hz * phis(ig)) / (hzu + wt * hz)
  lflux = (hzl * phisl(ig) + wtbar * hz * phis(ig)) / (hzl + wtbar * hz)
  wtg = wt * uflux / (wt * uflux + wtbar * lflux)
  wtgbar = 1. - wtg

  twtg0 = wtgbar * twt0
  twtgbar0 = wtgbar * twtbar0
  twtg1 = wtg * twt1
  twtgbar1 = wtg * twtbar1

  xss(ig, jgrp) = DynMacXsBen(iso0)%xss(itemp1, ig, jgrp) * twtg0 + DynMacXsBen(iso0)%xss(itemp2, ig, jgrp) * twtgbar0 + &
                  DynMacXsBen(iso1)%xss(jtemp1, ig, jgrp) * twtg1 + DynMacXsBen(iso1)%xss(jtemp2, ig, jgrp) * twtgbar1
END DO
END SUBROUTINE

SUBROUTINE xssmDynBen(itype, temp, igb, ige, jsfr, jsto, xss)
IMPLICIT NONE
INTEGER :: itype
REAL :: temp
INTEGER :: igb, ige, jsfr, jsto
REAL, POINTER :: xss(:, :)

REAL :: wt, wtbar
INTEGER :: temp1, temp2
INTEGER :: ig, jg

CALL tempInterpolation(itype, temp, wt, temp1, temp2)
wtbar = 1. - wt
DO ig = igb, ige
  DO jg = jsfr, jsto
    xss(jg, ig) = DynMacXsBen(itype)%xss(temp1, jg, ig) * wt + DynMacXsBen(itype)%xss(temp2, jg, ig) * wtbar
  END DO
END DO
END SUBROUTINE

SUBROUTINE xssm1DynBen(itype, temp, igb, ige, jsfr, jsto, xss)
IMPLICIT NONE
INTEGER :: itype
REAL :: temp
INTEGER :: igb, ige, jsfr, jsto
REAL, POINTER :: xss(:, :)

REAL :: wt, wtbar
INTEGER :: temp1, temp2
INTEGER :: ig, jg

CALL tempInterpolation(itype, temp, wt, temp1, temp2)
wtbar = 1. - wt
DO ig = igb, ige
  DO jg = jsfr, jsto
    xss(jg, ig) = DynMacXsBen(itype)%xss1(temp1, jg, ig) * wt + DynMacXsBen(itype)%xss1(temp2, jg, ig) * wtbar
  END DO
END DO
END SUBROUTINE

SUBROUTINE xssm2DynBen(itype, temp, igb, ige, jsfr, jsto, xss)
IMPLICIT NONE
INTEGER :: itype
REAL :: temp
INTEGER :: igb, ige, jsfr, jsto
REAL, POINTER :: xss(:, :)

REAL :: wt, wtbar
INTEGER :: temp1, temp2
INTEGER :: ig, jg

CALL tempInterpolation(itype, temp, wt, temp1, temp2)
wtbar = 1. - wt
DO ig = igb, ige
  DO jg = jsfr, jsto
    xss(jg, ig) = DynMacXsBen(itype)%xss2(temp1, jg, ig) * wt + DynMacXsBen(itype)%xss2(temp2, jg, ig) * wtbar
  END DO
END DO
END SUBROUTINE

SUBROUTINE xssm3DynBen(itype, temp, igb, ige, jsfr, jsto, xss)
IMPLICIT NONE
INTEGER :: itype
REAL :: temp
INTEGER :: igb, ige, jsfr, jsto
REAL, POINTER :: xss(:, :)

REAL :: wt, wtbar
INTEGER :: temp1, temp2
INTEGER :: ig, jg

CALL tempInterpolation(itype, temp, wt, temp1, temp2)
wtbar = 1. - wt
DO ig = igb, ige
  DO jg = jsfr, jsto
    xss(jg, ig) = DynMacXsBen(itype)%xss3(temp1, jg, ig) * wt + DynMacXsBen(itype)%xss3(temp2, jg, ig) * wtbar
  END DO
END DO
END SUBROUTINE

SUBROUTINE DnpBetaDynBen(itype, temp, beta)
IMPLICIT NONE
INTEGER :: itype
REAL :: temp, beta(DNP_nGRP)

REAL :: wt, wtbar
INTEGER :: temp1, temp2
INTEGER :: iprec

CALL tempInterpolation(itype, temp, wt, temp1, temp2)
wtbar = 1. - wt

DO iprec = 1, DNP_nGRP
  Beta(iprec) = DynMacXsBen(itype)%Beta(temp1, iprec) * wt + DynMacXsBen(itype)%Beta(temp2, iprec) * wtbar
END DO
END SUBROUTINE

SUBROUTINE NeutVeloDynBen(itype, temp, velo)
IMPLICIT NONE
INTEGER :: itype
REAL :: temp, velo(ngben)

REAL :: wt, wtbar
INTEGER :: temp1, temp2
INTEGER :: ig

CALL tempInterpolation(itype, temp, wt, temp1, temp2)
wtbar = 1. - wt

DO ig = 1, ngben
  velo(ig) = DynMacXsBen(itype)%velo(temp1, ig) * wt + DynMacXsBen(itype)%velo(temp2, ig) * wtbar
END DO

END SUBROUTINE

SUBROUTINE SetBoronCoolant_NEACRP(boronppm)
IMPLICIT NONE
REAL :: boronppm

REAL :: prevppm, delppm
INTEGER :: ixsl, ig

DO ixsl = 1, nxsltype_total
  IF(MacXsBen(ixsl)%lempty) CYCLE
  prevppm = MacXsBen(ixsl)%boronppm
  MacXsBen(ixsl)%boronppm = boronppm

  delppm = boronppm - prevppm
  MacXsBen(ixsl)%xstr(1) = MacXsBen(ixsl)%xstr(1) + delppm * NeaCrpXsBen(ixsl)%gradxs_c(1)
  MacXsBen(ixsl)%xstr(2) = MacXsBen(ixsl)%xstr(2) + delppm * NeaCrpXsBen(ixsl)%gradxs_c(6)
  MacXsBen(ixsl)%xsa(1) = MacXsBen(ixsl)%xsa(1) + delppm * NeaCrpXsBen(ixsl)%gradxs_c(3)
  MacXsBen(ixsl)%xsa(2) = MacXsBen(ixsl)%xsa(2) + delppm * NeaCrpXsBen(ixsl)%gradxs_c(7)
  MacXsBen(ixsl)%xsnf(1) = MacXsBen(ixsl)%xsnf(1) + delppm * NeaCrpXsBen(ixsl)%gradxs_c(4)
  MacXsBen(ixsl)%xsnf(2) = MacXsBen(ixsl)%xsnf(2) + delppm * NeaCrpXsBen(ixsl)%gradxs_c(8)
  MacXsBen(ixsl)%xskf(1) = MacXsBen(ixsl)%xskf(1) + delppm * NeaCrpXsBen(ixsl)%gradxs_c(5)
  MacXsBen(ixsl)%xskf(2) = MacXsBen(ixsl)%xskf(2) + delppm * NeaCrpXsBen(ixsl)%gradxs_c(9)
  MacXsBen(ixsl)%xss0(1,2) = MacXsBen(ixsl)%xss0(1,2) + delppm * NeaCrpXsBen(ixsl)%gradxs_c(2)
  MacXsBen(ixsl)%xss0(1,1) = MacXsBen(ixsl)%xstr(1) - MacXsBen(ixsl)%xsa(1) - MacXsBen(ixsl)%xss0(1,2) 
  MacXsBen(ixsl)%xss0(2,2) = MacXsBen(ixsl)%xstr(2) - MacXsBen(ixsl)%xsa(2)
  MacXsBen(ixsl)%xst(1) = MacXsBen(ixsl)%xstr(1)
  MacXsBen(ixsl)%xst(2) = MacXsBen(ixsl)%xstr(2)
  DO ig = 1, 2
    MacXsBen(ixsl)%Xs0sum(ig) = sum(MacXsBen(ixsl)%Xss(ig,1:ngben))
  END DO 
END DO

END SUBROUTINE

SUBROUTINE xsbaseben_NEACRP(ixsl_, rho, tm, td, xsmac)
IMPLICIT NONE
INTEGER :: ixsl_
REAL :: rho, tm, td
INTEGER :: igb, ige, jsfr, jsto
LOGICAL :: lscat1
TYPE(XsMac_Type) :: xsmac

REAL :: xstr(2), xsa(2), xsnf(2), xskf(2), xss0(2,2), xst(2), xs0sum(2)
REAL :: delrho, deltm, deltd
INTEGER :: ixsl
INTEGER :: i, ig, j

IF(ixsl_ .GT. nxsltype) THEN
  ixsl = NeacrpXsBen(ixsl_)%basexsl
ELSE
  ixsl = ixsl_
END IF

delrho = rho - NeaCrpXsBen(ixsl)%rho0
deltm = tm - NeaCrpXsBen(ixsl)%tm0
deltd = sqrt(td) - sqrt(NeaCrpXsBen(ixsl)%td0)

xstr(1) = MacXsBen(ixsl_)%xstr(1) + delrho * NeaCrpXsBen(ixsl)%gradxs_rho(1) + &
  deltm * NeaCrpXsBen(ixsl)%gradxs_tm(1) + deltd * NeaCrpXsBen(ixsl)%gradxs_td(1)
xstr(2) = MacXsBen(ixsl_)%xstr(2) + delrho * NeaCrpXsBen(ixsl)%gradxs_rho(6) + &
  deltm * NeaCrpXsBen(ixsl)%gradxs_tm(6) + deltd * NeaCrpXsBen(ixsl)%gradxs_td(6)
xsa(1) = MacXsBen(ixsl_)%xsa(1)   + delrho * NeaCrpXsBen(ixsl)%gradxs_rho(3) + &
  deltm * NeaCrpXsBen(ixsl)%gradxs_tm(3) + deltd * NeaCrpXsBen(ixsl)%gradxs_td(3)
xsa(2) = MacXsBen(ixsl_)%xsa(2)   + delrho * NeaCrpXsBen(ixsl)%gradxs_rho(7) + &
  deltm * NeaCrpXsBen(ixsl)%gradxs_tm(7) + deltd * NeaCrpXsBen(ixsl)%gradxs_td(7)
xsnf(1) = MacXsBen(ixsl_)%xsnf(1) + delrho * NeaCrpXsBen(ixsl)%gradxs_rho(4) + &
  deltm * NeaCrpXsBen(ixsl)%gradxs_tm(4) + deltd * NeaCrpXsBen(ixsl)%gradxs_td(4)
xsnf(2) = MacXsBen(ixsl_)%xsnf(2) + delrho * NeaCrpXsBen(ixsl)%gradxs_rho(8) + &
  deltm * NeaCrpXsBen(ixsl)%gradxs_tm(8) + deltd * NeaCrpXsBen(ixsl)%gradxs_td(8)
xskf(1) = MacXsBen(ixsl_)%xskf(1) + delrho * NeaCrpXsBen(ixsl)%gradxs_rho(5) + &
  deltm * NeaCrpXsBen(ixsl)%gradxs_tm(5) + deltd * NeaCrpXsBen(ixsl)%gradxs_td(5)
xskf(2) = MacXsBen(ixsl_)%xskf(2) + delrho * NeaCrpXsBen(ixsl)%gradxs_rho(9) + &
  deltm * NeaCrpXsBen(ixsl)%gradxs_tm(9) + deltd * NeaCrpXsBen(ixsl)%gradxs_td(9)
xss0(1,2) = MacXsBen(ixsl_)%xss0(1,2) + delrho * NeaCrpXsBen(ixsl)%gradxs_rho(2) + &
  deltm * NeaCrpXsBen(ixsl)%gradxs_tm(2) + deltd * NeaCrpXsBen(ixsl)%gradxs_td(2)
xss0(1,1) = xstr(1) - xsa(1) - xss0(1,2)
xss0(2,2) = xstr(2) - xsa(2)
xst(1) = xstr(1)
xst(2) = xstr(2)

DO ig = 1, 2
  xs0sum(ig) = sum(xss0(ig,1:ngben))
END DO

XsMac%lFuel = MacXsBen(ixsl_)%lFuel
DO i = 1, 2
  XsMac%xsmact(i) = xst(i)
  XsMac%xsmaca(i) = xsa(i)
  XsMac%xsmactr(i) = xstr(i)
  XsMac%xsmacnf(i) = xsnf(i)
  XsMac%xsmackf(i) = xskf(i)
  XsMac%chi(i) = MacXsBen(ixsl_)%chi(i)
END DO 

DO i = 1, 2
  DO j = 1, 2
      XsMac%xsmacsm(j, i) = xss0(j, i)
  END DO
  XsMac%xsmacs(i)=0
  DO j = 1, 2
    XsMac%xsmacs(i) = XsMac%xsmacs(i) + xss0(i, j)
  END DO 
  XsMac%xsmacstr(i) = XsMac%xsmacs(i)
  
END DO

END SUBROUTINE

END MODULE 