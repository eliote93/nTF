MODULE itrcntl_mod
USE param
IMPLICIT NONE
TYPE MOCItrCntl_Type
  INTEGER :: nitermax = 1
  INTEGER :: nitermin = 1
  REAL :: ResErr = 1.
  REAL :: eigverr = 1.
  REAL :: psierr = 1.
END TYPE

TYPE InSolverItrCntl_TYPE
  INTEGER :: ninmax = 20
  INTEGER :: ninmin = 3
  REAL :: convcrit = epsm5

END TYPE

TYPE CMFDItrCntl_TYPE
  INTEGER :: nitermax = 40
  INTEGER :: nitermin = 7
  INTEGER :: nCMFDpNodal = 5
  LOGICAL :: lfirstnodal = .TRUE.
  REAL :: convcrit = 0.1_8
  REAL :: AxDhatConvCrit = 1.0E-3_8
  REAL :: Res3dCmfdConv = epsm4
  REAL :: eigv3dCmfdconv = epsm5
END TYPE


TYPE GcCMFDItrCntl_TYPE
  INTEGER :: GcCmfdIt = 0
  INTEGER :: nitermax = 20
  INTEGER :: nitermin = 6
  REAL :: conv_rescrit = 0.1_8
  REAL :: conv_eigvcrit = 1.e-4
  INTEGER :: convmod = 1   !1: Residual Check,  2: eigchk 3: psichek
  LOGICAL :: lLogOut = .TRUE.
END TYPE
TYPE DcplXsGenCntl_Type
  REAL :: eigconv = epsm5
  REAL :: psiconv = epsm4
  REAL :: ResConv = epsm4
  INTEGER :: iRefTemp, iRefPln
  INTEGER :: nXsGenIterMax = 5
  INTEGER :: nXsGenFbIterMax(0:3)
  LOGICAL :: lMocConv = .FALSE.
END TYPE

TYPE DbgCntl_TYPE
  LOGICAL :: lPhiCNegChk = .FALSE.
  LOGICAL :: lNegFix = .FALSE.
END TYPE

TYPE ItrCntl_TYPE
  TYPE(MOCItrCntl_Type) :: MocItrCntl
  TYPE(InSolverItrCntl_TYPE) :: InSolverItrCntl
  TYPE(CMFDItrCntl_TYPE) :: CMFDItrCntl
  TYPE(GcCMFDItrCntl_TYPE) :: GcCMFDItrCntl
  TYPE(DcplXsGenCntl_Type) :: DcplXsGenCntl
  TYPE(CMFDItrCntl_TYPE) :: DcplCMFD3dItrCntl
  INTEGER :: SrcIt = 0
  INTEGER :: cmfdit = 0
  INTEGER :: GcCmfdIt = 0
  INTEGER :: mocit = 0
  INTEGER :: AxIt = 0
  INTEGER :: innerit = 0

  INTEGER :: SrcIt0 = 0
  INTEGER :: cmfdit0 = 0
  INTEGER :: GcCmfdIt0 = 0
  INTEGER :: mocit0 = 0
  INTEGER :: AxIt0 = 0

  INTEGER :: OuterMax = 20
  INTEGER :: DcplItrData(0:3,5)
  INTEGER :: nThCondiGen
  REAL :: eigconv = epsm5
  REAL :: psiconv = epsm4
  REAL :: resconv = epsm4
  REAL :: ThConv = epsm4
  REAL :: decuspconv = epsm3
  LOGICAL :: lconv = .FALSE.
  LOGICAL :: lnodal = .TRUE.
  LOGICAL :: lThConv = .FALSE.
  LOGICAL :: lRadDhatUpdt
  LOGICAL :: lGammaConv = .FALSE.   !--- CNJ Edit : Gamma Transport Convergence Flag
  LOGICAL :: lDecuspConv = .FALSE.
END TYPE

TYPE ConvItrCntl_TYPE
  INTEGER :: OuterMax = 15
  REAL :: eigconv = epsm5
  REAL :: psiconv = epsm4
  REAL :: resconv = epsm4
  REAL :: decuspconv = epsm3
END TYPE

save
LOGICAL :: IFcmfd = .TRUE.
INTEGER :: noutmax = 100
TYPE(ItrCntl_TYPE) :: ItrCntl
TYPE(ConvItrCntl_TYPE) :: ConvItrCntl
TYPE(ItrCntl_TYPE) :: DcplItrCntl(100)

#ifndef __GFORTRAN__
DATA ItrCntl%DcplItrData(0:3, 1) /1, 3, 2, 2/
DATA ItrCntl%DcplItrData(1:2, 2) /200, 40 /
DATA ItrCntl%DcplItrData(1:2, 3) /20, 10 /
DATA ItrCntl%DcplItrData(1:2, 4) /2, 2 /
DATA ItrCntl%DcplItrData(1:3, 5) /15, 30, 3 /
#endif

END module