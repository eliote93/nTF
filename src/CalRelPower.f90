#include <defines.h>
SUBROUTINE CalRelPower(Core, CmInfo, RelPower, ng, nTracerCntl, PE, lGather)
! Update Normalized Power for T/H

USE PARAM,          ONLY : ZERO, FALSE
USE TypeDef,        ONLY : CoreInfo_Type, CmInfo_Type, PE_Type, PinXS_Type, Pin_Type
USE CNTL,           ONLY : nTracerCntl_Type
USE BasicOperation, ONLY : CP_CA, CP_VA, DotProduct
USE TH_Mod,         ONLY : GatherRelPower

#ifdef __PGI
USE CUDA_MASTER,    ONLY : cuCntl, cuPwDist
USE CUDA_PWDIST,    ONLY : cuSetCMRelPower, cuUpdtPinPw
#endif

IMPLICIT NONE

TYPE(CoreInfo_Type)    :: Core
TYPE(CmInfo_Type)      :: CmInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE)          :: PE

INTEGER :: ng
REAL, POINTER, DIMENSION(:,:) :: RelPower
LOGICAL :: lGather
! ----------------------------------------------------
TYPE(PinXs_Type), POINTER, DIMENSION(:,:) :: PinXS
TYPE(Pin_Type),   POINTER, DIMENSION(:)   :: Pin

REAL, POINTER, DIMENSION(:)     :: hz
REAL, POINTER, DIMENSION(:,:)   :: PinVol
REAL, POINTER, DIMENSION(:,:,:) :: PhiC

REAL :: hzeff, inv_hzeff, Volsum, hactive, TotPow, localPow, vol, sum0, AvgPow, fnorm
INTEGER :: nch, icel, ig, iz, ixy, nxy, nz, myzb, myze
! ----------------------------------------------------

Pin    => Core%pin
PinVol => Core%PinVol
hz     => Core%hz
nxy     = Core%nxy
nz      = Core%nz

PinXs => CmInfo%PinXS
PhiC  => CmInfo%PhiC

myzb = PE%myzb
myze = PE%myze

! Obtain Active Height
hactive = ZERO
DO iz = 1, nz
  IF (Core%lFuelPlane(iz)) Hactive = Hactive + Core%hz(iz)
END DO
! ----------------------------------------------------
#ifdef __PGI
IF (PE%lCUDACMFD .AND. cuCntl%lPwDist) THEN
  IF (.NOT. cuPwDist%lTran) CALL cuUpdtPinPw(cuPwDist, FALSE, FALSE)
  
  CALL cuSetCMRelPower(cuPwDist, RelPower, PE, nz, lGather)
  
  hzeff = ZERO
  nch   = 0
  
  DO ixy = 1, nxy
    IF (.NOT. Pin(ixy)%lfuel) CYCLE
    
    nch = nch + 1
    
    DO iz = 1, nz
      icel = Pin(ixy)%cell(iz)
      
      IF (.NOT. Core%CellInfo(icel)%lfuel) CYCLE
      
      hzeff = hzeff + Core%hz(iz)
    END DO
  END DO
  
  hzeff     = hzeff / nch
  inv_hzeff = 1. / hzeff
  
  DO ixy = 1, nxy
    IF (.NOT. Pin(ixy)%lfuel) THEN
      RelPower(0, ixy) = ZERO
      
      CYCLE
    END IF
    
    sum0 = ZERO
    DO iz = 1, nz
      icel = Pin(ixy)%cell(iz)
      
      IF (.NOT. Core%CellInfo(icel)%lfuel) CYCLE
      
      sum0 = sum0 + RelPower(iz, ixy) * Core%hz(iz)
    END DO
    
    RelPower(0, ixy) = sum0 * inv_hzeff
  END DO
  
  RETURN
END IF
#endif
! ----------------------------------------------------
!--- CNJ Edit : MKL, CUDA Reorder Pin Cross Section
#ifdef __INTEL_MKL
IF (PE%lMKL) CALL MKL_ReorderPinXS(Core, CmInfo)
#endif

#ifdef __PGI
IF (PE%lCUDACMFD) CALL CUDAReorderPinXS(Core, CmInfo)
#endif
! ----------------------------------------------------
TotPow = ZERO
Volsum = ZERO

RelPower(1:nz, 1:nxy) = ZERO

DO iz = myzb, myze
  IF (.NOT. Core%lFuelPlane(iz)) CYCLE
  
  DO ixy = 1, nxy
    localPow = DotProduct(PinXS(ixy, iz)%xskf(1:ng), PhiC(ixy, iz, 1:ng), ng)
    
    IF (localPow .LE. ZERO) CYCLE
    
    vol      = PinVol(ixy, iz)
    localPow = localPow * vol
    
    RelPower(iz, ixy) = localPow
    
    Volsum = Volsum + vol
    TotPow = TotPow + LocalPow
  END DO
END DO

#ifdef MPI_ENV
IF (lGather) THEN
  CALL GatherRelPower(Core, RelPower, PE)
  
  DO iz = 1, nz
    IF (.NOT. Core%lFuelPlane(iz))  CYCLE
    IF ((iz-myzb)*(iz-myze) .LE. 0) CYCLE
    
    DO ixy = 1, nxy
      IF (RelPower(iz, ixy) .LE. 0_8) CYCLE
      
      Volsum = Volsum + PinVol(ixy, iz)
      TotPow = TotPow + RelPower(iz, ixy)
    END DO
  END DO
END IF
#endif
! ----------------------------------------------------
! Adj : Flux Lv.
AVGPOW = TotPow / VolSUM
fnorm  = 1._8 / AvgPow

DO ixy = 1, nxy
  IF (.NOT. Pin(ixy)%lfuel) CYCLE
  
  sum0 = ZERO
  DO iz =1 ,nz
    RelPower(iz, ixy) = RelPower(iz, ixy) * fnorm / PinVol(ixy, iz)
    
    SUM0 = SUM0 + RelPower(iz, ixy) * hz(iz)
  END DO
  
  RelPower(0, ixy) = SUM0 / Hactive
END DO

DO iz = 1, nz
  IF (.NOT. Core%lFuelPlane(iz)) CYCLE
  
  SUM0   = ZERO
  volsum = ZERO
  DO ixy = 1, nxy
    IF (.NOT. Pin(ixy)%lfuel) CYCLE
    
    sum0   = sum0   + RelPower(iz, ixy) * PinVol(ixy, iz)
    volsum = volsum + PinVol(ixy, iz)
  END DO
  
  RelPower(iz, 0) = Sum0 / VolSum
END DO
! ----------------------------------------------------
!#define relpowerout
#ifdef relpowerout
IF (PE%CMFDMASTER) THEN
  OPEN(89, FILE='REL.pw', status='replace')
  DO ixy = 1, nxy
    write(89, '(2I8, F15.5, 500e20.5)') Pin(ixy)%ix, Pin(ixy)%iy, sum(RelPower(:, ixy)),(RelPower(iz, ixy), iz=1, nz)
  END DO
  CLOSE(89)
END IF
#endif

NULLIFY (Pin)
NULLIFY (PinXS)
NULLIFY (PinVol)
NULLIFY (hz)
NULLIFY (PhiC)
! ----------------------------------------------------

END SUBROUTINE CalRelPower
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SubCellPowerProfile(Core, FmInfo, Profile, ixy, npr, nz, ng, GroupInfo, PE)

USE PARAM,          ONLY : ZERO, FALSE
USE TYPEDEF,        ONLY : CoreInfo_Type, FmInfo_Type, PE_TYPE, GroupInfo_Type, FxrInfo_Type, Cell_Type, Pin_Type, XsMac_Type
USE MacXsLib_Mod,   ONLY : MacXskf
USE BasicOperation, ONLY : CP_CA, MULTI_VA

#ifdef MPI_ENV
USE MPICOMM_MOD,  ONLY : BCAST
#endif

IMPLICIT NONE

TYPE (CoreInfo_Type)  :: Core
TYPE (FmInfo_Type)    :: FmInfo
TYPE (GroupInfo_Type) :: GroupInfo
TYPE (PE_Type)        :: PE
REAL, POINTER, DIMENSION(:,:) :: Profile
INTEGER :: ixy, ng, nz, npr
! ----------------------------------------------------
TYPE(FxrInfo_Type), POINTER, DIMENSION(:,:) :: Fxr
TYPE(Cell_Type),    POINTER, DIMENSION(:)   :: Cell
TYPE(Pin_Type),     POINTER, DIMENSION(:)   :: Pin

REAL, POINTER, DIMENSION(:)     :: xsmackf
REAL, POINTER, DIMENSION(:,:,:) :: phis

TYPE(XsMac_Type), SAVE :: XsMac

REAL :: psum, volsum, fxrprofile(1:200)
INTEGER :: myzb, myze, FsrIdxSt, FxrIdxSt, nlocalFxr, nlocalFsr, nFsrInFxr, norg
INTEGER :: icel, iz, iz1, iz2, ig, ifxr, ifsr, ifsrlocal, iResoGrpBeg, iResoGrpEnd, i, j
! ----------------------------------------------------

myze = PE%myze
myzb = PE%myzb

iResoGrpBeg = GroupInfo%nofg + 1
iResoGrpEnd = GroupInfo%nofg + GroupInfo%norg
norg        = GroupInfo%norg

Fxr  => FmInfo%Fxr
Phis => FmInfo%Phis

Pin  => Core%Pin
Cell => Core%CellInfo

FsrIdxSt = Pin(ixy)%FsrIdxSt
FxrIdxSt = Pin(ixy)%FxrIdxSt

DO iz = myzb, myze
  icel = Pin(ixy)%Cell(iz)
  
  nlocalFxr = Cell(icel)%nFxr
  nlocalFsr = Cell(icel)%nFsr
  
  volsum = ZERO
  
  DO j = 1, nLocalFxr
    fxrprofile(j) = 0
    ifxr = FxrIdxSt + j -1
    
    IF (.NOT. Fxr(ifxr,iz)%lfuel) CYCLE
    
    nFsrInFxr = Cell(icel)%nFsrInFxr(j)
    
    CALL MacXsKf(XsMac, Fxr(ifxr, iz), 1, ng, ng, 1._8, FALSE)
    
    xsmackf => XsMac%XsMackf
    
    DO ig = iResoGrpBeg,iResoGrpEnd
      XsMackf(ig) = XsMackf(ig) * Fxr(ifxr, iz)%fresokF(ig)
    END DO
    
    DO i = 1, nFsrInFxr
      ifsrlocal = Cell(icel)%MapFxr2FsrIdx(i, j)
      ifsr      = FsrIdxSt + ifsrlocal - 1
      
      DO ig = 1, ng
        fxrprofile(j) = fxrprofile(j) + Cell(icel)%vol(ifsrlocal) * phis(ifsr, iz, ig) * XsMacKf(ig)
      END DO
    END DO
    
    IF (fxrprofile(j) .GT. ZERO) volsum = volsum + Fxr(ifxr, iz)%area
  END DO
  
  psum = volsum / SUM(fxrprofile(1:nLocalFxr))
  
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j - 1
    
    fxrprofile(j) = fxrprofile(j) * psum / Fxr(ifxr, iz)%area
  END DO
  
  DO j = 1, npr
    Profile(j,     iz) = 0
    Profile(j+npr, iz) = 0
    DO i = 1, nlocalFxr
      Profile(j,     iz) = Profile(j,     iz) + Cell(icel)%ThCell%InvFrac(i, j) * FxrProfile(i)
      Profile(j+npr, iz) = Profile(j+npr, iz) + Cell(icel)%ThCell%InvFrac(i, j)
    END DO
    
    CONTINUE
  END DO
END DO

#ifdef MPI_ENV
DO i = 0, PE%nproc-1
  iz1 = PE%AxDomRange(1, i)
  iz2 = PE%AxDomRange(2, i)
  
  CALL BCAST(Profile(1:2*npr, iz1:iz2), 2*npr, iz2-iz1 + 1, PE%MPI_CMFD_COMM, i)
END DO
#endif
! ----------------------------------------------------

END SUBROUTINE SubCellPowerProfile
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE GatherRelPower(Core, RelPow, PE)

USE TYPEDEF, ONLY : CoreInfo_Type, PE_TYPE

#ifdef MPI_ENV
USE MPIComm_mod, ONLY : BCAST, MPI_SYNC
#endif

IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(PE_TYPE)       :: PE
REAL, POINTER, DIMENSION(:,:) :: RelPow
! ----------------------------------------------------
REAL, POINTER :: Buf(:, :)

#ifdef MPI_ENV
INTEGER :: myzb, myze, nxy, nz, comm, nproc, i, j, iz1, iz2
! ----------------------------------------------------

nxy = Core%nxy
nz  = Core%nz

myzb  = PE%myzb
myze  = PE%myze
nproc = PE%nCmfdProc
comm  = PE%MPI_CMFD_COMM

ALLOCATE(Buf(nxy, nz))

DO j = 1, nz
  DO i = 1, nxy
    Buf(i, j) = RelPow(j, i)
  END DO
END DO

DO i = 0, nproc-1
  iz1 = PE%AxDomRange(1, i)
  iz2 = PE%AxDomRange(2, i)
  
  CALL BCAST(Buf(1:nxy, iz1:iz2), nxy, iz2 - iz1 + 1, comm, i)
END DO

DO j = 1, nz
  DO i = 1, nxy
    RelPow(j, i) = Buf(i, j)
  END DO
END DO

DEALLOCATE(Buf)
#endif
! ----------------------------------------------------

  END SUBROUTINE GatherRelPower
! ------------------------------------------------------------------------------------------------------------
!SUBROUTINE Grp_RelPower(Core, CmInfo, RelPower, RelPower_Grp, ng, nTracerCntl, PE, lGather)
SUBROUTINE Grp_RelPower(Core, CmInfo, RelPower, ng, nTracerCntl, PE, lGather)
!Update Normalized Power for T/H

USE PARAM,          ONLY : ZERO
USE TypeDef,        ONLY : CoreInfo_Type, CmInfo_Type, PE_Type, PinXS_Type, Pin_Type, Asy_Type, AsyInfo_Type
USE CNTL,           ONLY : nTracerCntl_Type
USE BasicOperation, ONLY : CP_CA, CP_VA, DotProduct
USE TH_Mod,         ONLY : GatherRelPower

IMPLICIT NONE

TYPE(CoreInfo_Type)    :: Core
TYPE(CmInfo_Type)      :: CmInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE)          :: PE
INTEGER :: ng
REAL, POINTER, DIMENSION(:,:) :: RelPower!, RelPower_Grp
LOGICAL :: lGather
! ----------------------------------------------------
TYPE(Pin_Type),     POINTER, DIMENSION(:) :: Pin
TYPE(AsyInfo_Type), POINTER, DIMENSION(:) :: AsyInfo
TYPE(Asy_Type),     POINTER, DIMENSION(:) :: Asy

REAL, POINTER, DIMENSION(:)   :: hz
REAL, POINTER, DIMENSION(:,:) :: PinVol

REAL, DIMENSION(1000) :: PW_ThChGrp, Vol_ThChGrp
INTEGER :: ixy, ixy0, iasy, iz, iasytype, igrp, nxy, nz, nasytype, nasy
! ----------------------------------------------------

Pin     => Core%Pin
PinVol  => Core%PinVol
Hz      => Core%hz
AsyInfo => Core%AsyInfo
Asy     => Core%Asy
nxy      = Core%nxy
nz       = Core%nz
nasytype = Core%nAsyType
nasy     = Core%nxya

DO iz = 1, nz
  DO iasy = 1, nasy
    iasytype = Asy(iasy)%AsyType
    
    PW_ThChGrp  = ZERO
    Vol_ThChGrp = ZERO
    DO ixy0 = 1, AsyInfo(iasytype)%nxy
      igrp = AsyInfo(iasytype)%ThChGrp(ixy0)
      ixy  = Asy(iasy)%GlobalPinIdx(ixy0)
      
      PW_ThChGrp (igrp) = PW_ThChGrp (igrp) + PinVol(ixy, iz) * RelPower(iz, ixy)
      Vol_ThChGrp(igrp) = Vol_ThChGrp(igrp) + PinVol(ixy, iz)
    END DO
    
    DO igrp = 1, AsyInfo(iasytype)%nThChGrp
      IF (Vol_ThChGrp(igrp) .LT. 1.0E-5) CYCLE
      
      PW_ThChGrp(igrp) = PW_ThChGrp(igrp) / Vol_ThChGrp(igrp)
    END DO
    
    DO ixy0 = 1, AsyInfo(iasytype)%nxy
      igrp = AsyInfo(iasytype)%ThChGrp(ixy0)
      ixy  = Asy(iasy)%GlobalPinIdx(ixy0)
      
      !RelPower_Grp(iz, ixy) = Pw_ThChGrp(igrp)
      RelPower(iz, ixy) = Pw_ThChGrp(igrp)
    END DO
  END DO
END DO

NULLIFY(Pin, PinVol, Hz)
NULLIFY(AsyInfo, Asy)
! ----------------------------------------------------

END SUBROUTINE Grp_RelPower
! ------------------------------------------------------------------------------------------------------------