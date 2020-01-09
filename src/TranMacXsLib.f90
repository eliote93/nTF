MODULE TranMacXsLib_Mod
USE PARAM
USE TYPEDEF,         ONLY : XsMac_Type,     GroupInfo_Type,      Fxrinfo_type
USE XSLIB_MOD
USE XsUtil_mod,      ONLY : AllocXsMac,          AllocMacIsoXs,  &
                            LineIntPol,          XsTempInterpolation, P1XsTempInterpolation
USE BasicOperation,  ONLY : CP_CA, CP_VA,   MULTI_CA,            DotProduct 
USE IOUTIL,          ONLY : TERMINATE
!Default Kinetic Parameter for XS Library Calculation Mode
INTEGER, PARAMETER :: nprec_LibBase = 6
INTEGER, PARAMETER :: ng_LibBase = 47;

REAL :: DecayConst(6), DelaySpectrum(ng_LibBase)
DATA DecayConst / 0.0128_8, 0.0318_8, 0.119_8, 0.3181_8, 1.4027_8, 3.9286_8/
DATA DelaySpectrum / 4*0._8, 0.005_8, 0.021_8, 0.269_8, 0.247_8, 0.429_8, 0.029_8, 37*0._8/

CONTAINS

SUBROUTINE SetTranXsLibInfo(GroupInfo)
TYPE(GroupInfo_Type) :: GroupInfo
GroupInfo%nPrec = nprec_LibBase
IF(GroupInfo%ng .NE. ng_LibBase) THEN
  CALL TERMINATE('For Transient Calculation w/ XS Lib., the 47 Group Structure is allowed only')
ENDIF
END SUBROUTINE

SUBROUTINE InitChidLib(Chid, ng)
REAL :: Chid(ng)
INTEGER :: ng

Chid(1:ng) = DelaySpectrum
END SUBROUTINE


SUBROUTINE InitLambdaLib(Lambda)
REAL :: Lambda(nPrec_LibBase)
Lambda = DecayConst
END SUBROUTINE

SUBROUTINE FxrBeta(XsMac, Fxr, FxrPhi, ng)
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac
TYPE(FxrINfo_Type) :: Fxr
REAL :: FxrPhi(ng)
INTEGER :: ng
!EffMacXs and FxrAvgPhi should be called prior to this routine
INTEGER, POINTER :: IdIso(:)
REAL, POINTER :: pnum(:)

REAL :: beta(nprec_LibBase)

REAL :: FisIso                !Fission Rate of One Isotope
REAL :: FisSum

INTEGER :: niso
INTEGER :: ig, iso, id, iprec, ifis

niso = Fxr%niso
IdIso => Fxr%IdIso; pnum => Fxr%pnum
FisSum = 0; Beta = 0
DO iso = 1, niso
  id = IdIso(iso);  ifis = MapFis(id)
  IF(ifis .EQ. 0) CYCLE
  id = MapNucl(idiso(iso)); 
  FisIso = 0;
  DO ig = 1, ng
    !FisIso = FisIso + XsMac%IsoXsMacNf(iso, ig) * FxrPhi(ig)  
    FisIso = FisIso + XsMac%IsoXsMacNf(iso, ig) * Fxr%fresoFIso(iso,ig) * FxrPhi(ig)
  ENDDO
  FisSum = FisSum + FisIso
  
  DO iprec = 1, nprec_LibBase
    Beta(iprec) = Beta(iprec) + FisIso * ldiso(id)%Beta(iprec)
  ENDDO
ENDDO
IF(FisSum .GT. 0) THEN
  Fissum = 1._8 / FisSum
  DO iprec = 1, nprec_LibBase
    Beta(iprec) = Beta(iprec) * FisSum 
  ENDDO
ENDIF
Fxr%Beta = Beta
END SUBROUTINE

SUBROUTINE FxrVelo(Fxr, Temp, ng)
IMPLICIT NONE
TYPE(FxrInfo_Type) :: Fxr
REAL :: Temp
INTEGER :: ng

REAL :: sigab10, velot, ubar, ebar
REAL :: wt1, wt2
INTEGER :: it1, it2
INTEGER :: ig, id

id = MapNucl(5010)
CALL  XsTempInterpolation(id, ldiso(id), temp, wt1, wt2, it1, it2)

DO ig = 1, ng
  IF(enbhel(ig) .LT. 0.1e+6) THEN  !less than 0.1 MeV 
    sigab10 = (wt2 * ldiso(id)%siga(ig, it2) + wt1 * ldiso(id)%siga(ig, it1))
    velot = 3837._8 * 2.2e5_8/sigab10
  ELSE
    ubar=0.5_8*(uhel(ig)+uhel(ig+1))
    ebar=1.0e7_8/exp(ubar)
    velot=2.2e5_8*sqrt(ebar/0.0253_8)    
  ENDIF
  Fxr%Velo(ig) = Velot
ENDDO

END SUBROUTINE


END MODULE