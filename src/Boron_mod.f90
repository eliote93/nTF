MODULE Boron_mod
REAL, SAVE :: ppmd, ppmd2
REAL, SAVE :: eigvd = 1._8
INTEGER, SAVE :: iter = 0 
REAL, SAVE :: b10frac=0.198_8
REAL, SAVE :: b10frac0=0.198_8
REAL :: sigc1g_b10(0:500), phi1g_mod(0:500), DeplB10frac(0:500)
REAL :: BoronPPM(0:500), DeplBoronPPM(0:500)
INTERFACE 

SUBROUTINE SetBoronCoolant(Core, Fxr, boronppm, myzb, myze)
USE PARAM
USE TypeDef,       ONLY : CoreInfo_type, Fxrinfo_type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
REAL :: boronppm
INTEGER :: myzb, myze
END SUBROUTINE

SUBROUTINE SetBoronCoolantCTF(Core, Fxr,THinfo, myzb, myze)
USE PARAM
USE TypeDef,       ONLY : CoreInfo_type,THInfo_type, Fxrinfo_type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(THinfo_Type) :: THinfo
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
INTEGER :: myzb, myze
END SUBROUTINE

SUBROUTINE UpdtBoronPPM(target_eigv, eigv, ppm, lreset, master)
IMPLICIT NONE
REAL :: target_eigv, eigv, ppm
LOGICAL :: lreset, master
END SUBROUTINE


SUBROUTINE UpdtBoronCmfdXS(Core, Fxr, Phis, PinXS, boronppm, myzb, myze, ng)
USE PARAM
USE TypeDef,       ONLY : CoreInfo_type, Fxrinfo_type, PinXS_Type,   &
                          Pin_Type,      Cell_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(PinXS_Type), POINTER :: PinXS(:, :)
REAL, POINTER :: Phis(:, :, :)

REAL :: boronppm
INTEGER :: myzb, myze, ng 

END SUBROUTINE

SUBROUTINE CalB10XS(Core, FmInfo, istep, ng, GroupInfo, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,       ONLY : CoreInfo_Type,    FmInfo_Type,     GroupInfo_Type,     &
                          PE_Type
USE CNTL,          ONLY : nTracerCntl_Type                          
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE
INTEGER :: ng
INTEGER :: istep

END SUBROUTINE

SUBROUTINE PostB10Depl(DeplCntl, PE)
USE PARAM
USE TYPEDEF,    ONLY : PE_TYPE
USE DeplType,   ONLY : DeplCntl_Type
IMPLICIT NONE

TYPE(DeplCntl_Type) :: DeplCntl
TYPE(PE_TYPE) :: PE
END SUBROUTINE

SUBROUTINE PostB10Depl0(nowstep, DeplCntl, PE)
USE PARAM
USE TYPEDEF,    ONLY : PE_TYPE
USE DeplType,   ONLY : DeplCntl_Type
IMPLICIT NONE
INTEGER :: nowstep
TYPE(DeplCntl_Type) :: DeplCntl
TYPE(PE_TYPE) :: PE
END SUBROUTINE


END INTERFACE

CONTAINS

SUBROUTINE SaveBoronPPM(istep, ppm)
IMPLICIT NONE
REAL :: ppm
INTEGER :: istep
BoronPPM(istep) = ppm
END SUBROUTINE

FUNCTION GetBoronPPM(istep)
REAL :: GetBoronPPM
INTEGER :: istep
GetBoronPPM = BoronPPM(istep)
END FUNCTION

END MODULE