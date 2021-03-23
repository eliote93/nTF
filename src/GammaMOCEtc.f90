#include <defines.h>    
!--- CNJ Edit : Auxiliary Routines for Gamma MOC Calculation
#ifdef __GAMMA_TRANSPORT    
FUNCTION GammaMocResidual(Core, FmInfo, GroupInfo, PE, nTracerCntl)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,           FmInfo_Type,            GroupInfo_Type,         &
                           Cell_Type,               Pin_Type,               FxrInfo_Type,           &
                           PE_Type
USE GammaCore_mod,  ONLY : gphis,                   gJout,                  gAxSrc,                 &
                           gAxPXS
USE GamMOC_MOD,     ONLY : SetGamSrc,               SetGamSrcNM,            SetGamMacXs,            &
                           SetGamMacXsNM,           GamPseudoAbsorption,    GamPseudoAbsorptionNM,  &
                           gxst1g,                  gsrc1g,                 gAxSrc1g,               &
                           gAxPXS1g,                gphisnm,                gxstnm,                 &
                           gsrcnm
USE CNTL,           ONLY : nTracerCntl_Type
#ifdef MPI_ENV
USE MPICOMM_MOD,    ONLY : REDUCE
#endif
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(PE_TYPE) :: PE
TYPE(nTracerCntl_Type) :: nTracerCntl
REAL :: GammaMocResidual

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
REAL, POINTER :: phis(:, :, :), phisnm(:, :)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
INTEGER :: ng, ngg, nxy, nFsr, nlocalFsr, FsrIdxSt, myzb, myze
INTEGER :: i, ig, ireg, ipin, icel, iz, ifsrlocal, nbd
LOGICAL :: l3dim, lscat1, lTrCorrection
REAL :: SrcSum, LocalResidual, localsrc, vol, temp
REAL :: Leakage(100), Collision(100), Source(100)

Fxr => FmInfo%Fxr
phis => FmInfo%phis

Pin => Core%Pin
Cell => Core%CellInfo
l3dim = nTracerCntl%l3dim
lscat1 = nTracerCntl%lGammaScat1
lTrCorrection = .NOT. lscat1
myzb = PE%myzb; myze = PE%myze
nFsr = Core%nCoreFsr; nxy = Core%nxy; nbd = 4; ng = GroupInfo%ng; ngg = GroupInfo%ngg
GammaMocResidual = 0; SrcSum = 0
IF (lscat1) THEN

DO ig = 1, ngg
  DO iz = myzb, myze
    IF (.NOT. Core%lFuelPlane(iz)) CYCLE
    gAxSrc1g => gAxSrc(:, ig, iz)
    gAxPXS1g => gAxPXS(:, ig, iz)
    CALL SetGamMacXs(Core, Fxr(:, iz), gxst1g, iz, ig, ngg, lTrCorrection, PE)
#ifdef LkgSplit
    CALL GamPseudoAbsorption(Core, Fxr(:, iz), gAxPXS1g, gxst1g, iz, l3dim)
#endif            
    CALL SetGamSrc(Core, Fxr(:, iz), gsrc1g, phis, gphis, gAxSrc1g, gxst1g, iz, ig, ng, ngg,                   &
                   l3dim, lscat1, PE, GroupInfo)
    DO ipin = 1, nxy
      FsrIdxSt = Pin(ipin)%FsrIdxSt
      icel = Pin(ipin)%Cell(iz)
      nlocalFsr = Cell(icel)%nFsr
      LocalResidual = 0; localsrc = 0
      DO i = 1, nbd
        LocalResidual = LocalResidual + (gJout(2, i, ipin, ig, iz) - gJout(1, i, ipin, ig, iz))
      ENDDO      
      DO i = 1, nlocalFsr
        ireg = FsrIdxSt + i - 1
        vol = Cell(icel)%vol(i)
        LocalResidual = LocalResidual + vol * gphis(ireg, ig, iz) * gxst1g(ireg)
        localsrc = localsrc + gsrc1g(ireg) * vol * gxst1g(ireg)
      ENDDO
      LocalResidual = localsrc - LocalResidual
      SrcSum = SrcSum + localsrc * localsrc
      GammaMocResidual = GammaMocResidual + LocalResidual * LocalResidual
    ENDDO
  ENDDO
ENDDO

ELSE

ALLOCATE(phisnm(ng, nFsr))

DO iz = myzb, myze
  IF (.NOT. Core%lFuelPlane(iz)) CYCLE
  DO ig = 1, ng
    phisnm(ig, :) = phis(:, iz, ig)
  ENDDO
  gphisnm => gphis(:, :, iz)
  CALL SetGamMacXsNM(Core, Fxr(:, iz), gxstnm, iz, ngg, lTrCorrection, PE)
#ifdef LkgSplit
  CALL GamPseudoAbsorptionNM(Core, Fxr(:, iz), gAxPXS, gxstnm, iz, ngg, l3dim)
#endif       
  CALL SetGamSrcNM(Core, Fxr(:, iz), gsrcnm, phisnm, gphisnm, gAxSrc, gxstnm, iz,                               &
                   1, ngg, ng, ngg, l3dim, lscat1, PE)
  DO ipin = 1, nxy
    FsrIdxSt = Pin(ipin)%FsrIdxSt
    icel = Pin(ipin)%Cell(iz)
    nlocalFsr = Cell(icel)%nFsr
    Collision = 0; Source = 0; Leakage = 0
    DO i = 1, nbd
      DO ig = 1, ngg
        Leakage(ig) = Leakage(ig) + gJout(2, ig, i, ipin, iz) - gJout(1, ig, i, ipin, iz)
      ENDDO
    ENDDO      
    DO i = 1, nlocalFsr
      ireg = FsrIdxSt + i - 1
      vol = Cell(icel)%vol(i)
      DO ig = 1, ngg
        Collision(ig) = Collision(ig) + vol * gphis(ig, ireg, iz) * gxstnm(ig, ireg)
        Source(ig) = Source(ig) + gsrcnm(ig, ireg) * vol * gxstnm(ig, ireg)
      ENDDO
    ENDDO
    DO ig = 1, ngg
      LocalResidual = Source(ig) - Collision(ig) - Leakage(ig)
      SrcSum = SrcSum + Source(ig) ** 2
      GammaMocResidual = GammaMocResidual + LocalResidual ** 2
    ENDDO
  ENDDO
ENDDO

DEALLOCATE(phisnm)

ENDIF

#ifdef MPI_ENV
CALL REDUCE(GammaMocResidual, temp, PE%MPI_RTMASTER_COMM, .TRUE.)
GammaMocResidual = temp
CALL REDUCE(SrcSum, temp, PE%MPI_RTMASTER_COMM, .TRUE.)
SrcSum = temp
#endif
GammaMocResidual = GammaMocResidual / SrcSum
GammaMocResidual = SQRT(GammaMocResidual)

END FUNCTION

!--- JSU Edit : Subroutine to Calculate Gamma Power 
SUBROUTINE GammaPowerUpdate(Core, Fxr, gphis, gpower, myzb, myze, ng, PE)
USE PARAM
USE TYPEDEF,      ONLY : coreinfo_type,       Fxrinfo_type,       Cell_Type,     pin_Type, &
                         PE_TYPE
USE CNTL,         ONLY : nTracerCntl
USE GamXsLib_Mod, ONLY : GamXsBase
#ifdef MPI_ENV
USE MPIComm_Mod, ONLY : BCAST
#endif
USE GammaTYPEDEF, ONLY : GamMacXS_TYPE

IMPLICIT NONE

TYPE(coreinfo_type) :: CORE
TYPE(Fxrinfo_type),POINTER :: Fxr(:, :)
TYPE(PE_Type) :: PE
REAL, POINTER :: gphis(:, :, :)
REAL, POINTER :: gpower(:, :)
INTEGER :: myzb, myze, ng


TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(Fxrinfo_type),POINTER :: myFxr
TYPE(GamMacXS_TYPE), SAVE :: GamXsMac

INTEGER :: nxy, nFsr, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr
INTEGER :: i, j, ipin, icel, ifsr, ifxr, iz, ig

REAL, POINTER :: xsmackerma(:)

Pin => Core%Pin
CellInfo => Core%CellInfo
nFsr = Core%nCoreFsr
nxy = Core%nxy
gPower =0.
IF (nTracerCntl%lGammaScat1) THEN
  DO iz = myzb, myze
    DO ig = 1, ng
      DO ipin = 1, nxy
        FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
        icel = Pin(ipin)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr  
        DO j = 1, nLocalFxr
          ifxr = FxrIdxSt + j -1
          nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)    
          myFxr => Fxr(ifxr, iz)
          CALL GamXsBase(GamXsMac, myFxr, ig, ig, ng, FALSE)
          xsmackerma => GamXsMac%MacKERMA
          
          DO i = 1, nFsrInFxr
            ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1  !Global FSR Index
            gPower(ifsr, iz) = gPower(ifsr, iz) + xsmackerma(ig) * gphis(ifsr, ig, iz)
          ENDDO
        ENDDO !Fsr Sweep    
      ENDDO
    ENDDO
  ENDDO
ELSE
  DO iz = myzb, myze
    DO ipin = 1, nxy
      FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
      icel = Pin(ipin)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr  
      DO j = 1, nLocalFxr
        ifxr = FxrIdxSt + j -1
        nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)    
        myFxr => Fxr(ifxr, iz)
        CALL GamXsBase(GamXsMac, myFxr, 1, ng, ng, FALSE)
        xsmackerma => GamXsMac%MacKERMA
        
        DO i = 1, nFsrInFxr
          ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1  !Global FSR Index
          DO ig = 1, ng
            gPower(ifsr, iz) = gPower(ifsr, iz) + xsmackerma(ig) * gphis(ig, ifsr, iz)
          ENDDO
        ENDDO !Fsr Sweep    
      ENDDO
    ENDDO
  ENDDO
ENDIF    

#ifdef MPI_ENV
DO iz = 1, Core%nz
  CALL BCAST(gPower(:, iz), nFsr, PE%MPI_RTMASTER_COMM, PE%AxDomList(iz))
ENDDO
#endif

END SUBROUTINE
  
SUBROUTINE NeutronLocalQUpdate(Core, Fxr, phis, localpower, GPOWERGEN, myzb, myze, ng, PE)
USE PARAM
USE TYPEDEF,        ONLY : coreinfo_type,       Fxrinfo_type,       Cell_Type,     pin_Type, &
                           PE_TYPE,             XsMac_TYPE
USE CNTL,           ONLY : nTracerCntl
USE CORE_MOD,       ONLY : GroupInfo
USE GamXsLib_Mod,   ONLY : GetLocalQMat
#ifdef MPI_ENV      
USE MPIComm_Mod,    ONLY : BCAST
#endif              
IMPLICIT NONE

TYPE(coreinfo_type) :: CORE
TYPE(Fxrinfo_type),POINTER :: Fxr(:, :)
TYPE(PE_Type) :: PE
REAL, POINTER :: phis(:, :, :)
REAL, POINTER :: localpower(:, :), GPOWERGEN(:, :)
INTEGER :: myzb, myze, ng, ngg

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(Fxrinfo_type),POINTER :: myFxr

TYPE(XsMac_TYPE) :: XsMac

INTEGER :: nxy, nCoreFsr, nCoreFxr, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr
INTEGER :: ipin, icel, ifsr, ifxr, iz, ig
INTEGER :: i, j

REAL, POINTER :: MacLocalQ(:)
LOGICAL :: lfis

Pin => Core%Pin
CellInfo => Core%CellInfo; nCoreFsr = Core%nCoreFsr
nCoreFxr = Core%nCoreFxr; nxy = Core%nxy

ngg = GroupInfo%ngg

localpower = 0.
GPOWERGEN = 0.
DO iz = myzb, myze
  DO ipin = 1, nxy
    FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
    icel = Pin(ipin)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr  
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j -1
      nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)    
      myFxr => Fxr(ifxr, iz)
      CALL GetLocalQMat(XSMac, myFxr, 1, ng, ng, .TRUE.)
      DO i = 1, nFsrInFxr
        ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1  !Global FSR Index
        DO ig = 1, ng
          localpower(ifsr, iz) = localpower(ifsr, iz) + XsMac%MacKERMA_t(ig) * phis(ifsr, iz, ig)
          GPOWERGEN(ifsr, iz)        = GPOWERGEN(ifsr, iz)        + XsMac%MacKERMA_P(ig) * phis(ifsr, iz, ig)
        ENDDO
        CONTINUE
      ENDDO !Fsr Sweep    
    ENDDO
  ENDDO
ENDDO
#ifdef MPI_ENV
DO iz = 1, Core%nz
  CALL BCAST(localpower(1:nCoreFsr, iz), nCoreFsr, PE%MPI_RTMASTER_COMM, PE%AxDomList(iz))
ENDDO
#endif
END SUBROUTINE
  
FUNCTION FxrAvgGPhi(Core, Fxr, GPhis, ipin, iLocalfxr, iz, ng, PE)
USE PARAM
USE TYPEDEF,     ONLY : CoreInfo_Type,    PE_Type,     FxrInfo_Type,               &
                        Cell_Type,    Pin_Type
USE CNTL,        ONLY : nTRACERCNTL
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(PE_Type) :: PE
REAL, POINTER :: GPhis(:, :, :)
INTEGER :: iLocalfxr, ipin, iz, ng
REAL :: FxrAvgGPhi(ng)

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)

REAL :: area, areasum
INTEGER :: FsrIdxSt, FxrIdxSt, nFsrInFxr
INTEGER :: ifxr, ifsr, icell
INTEGER :: i, j

CellInfo => Core%CellInfo
Pin => Core%Pin

FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt

ifxr = FxrIdxSt + iLocalFxr - 1; icell = Pin(ipin)%Cell(iz)
nFsrInFxr = Fxr(ifxr, iz)%nFsrInFxr
FxrAvgGPhi = 0; areasum = 0
IF (nTracerCntl%lGammaScat1) THEN
  DO i = 1, nFsrInFxr
    j = CellInfo(icell)%MapFxr2FsrIdx(i, iLocalFxr)
    iFsr = FsrIdxSt + j - 1
    Area = CellInfo(icell)%vol(j); AreaSum = AreaSum + Area
    FxrAvgGPhi = FxrAvgGPhi + Area * GPhis(ifsr, 1:ng, iz)
  ENDDO
ELSE
  DO i = 1, nFsrInFxr
    j = CellInfo(icell)%MapFxr2FsrIdx(i, iLocalFxr)
    iFsr = FsrIdxSt + j - 1
    Area = CellInfo(icell)%vol(j); AreaSum = AreaSum + Area
    FxrAvgGPhi = FxrAvgGPhi + Area * GPhis(1:ng, ifsr, iz)
  ENDDO
END IF
FxrAvgGPhi = FxrAvgGPhi / AreaSum

NULLIFY(Pin, CellInfo)
END FUNCTION  

SUBROUTINE CompensateGPower(Core,Fxr,GPowerGen,GPower, myzb, myze,PE)
USE PARAM
USE TYPEDEF,      ONLY : coreinfo_type,       Fxrinfo_type,       Cell_Type,     pin_Type, &
                         PE_TYPE
USE CNTL,         ONLY : nTracerCntl
USE GamXsLib_Mod, ONLY : GamXsBase
#ifdef MPI_ENV
USE MPIComm_Mod, ONLY : BCAST
#endif
USE GammaTYPEDEF, ONLY : GamMacXS_TYPE

IMPLICIT NONE

TYPE(coreinfo_type) :: CORE
TYPE(Fxrinfo_type),POINTER :: Fxr(:, :)
TYPE(PE_Type) :: PE
REAL, POINTER, DIMENSION(:, :) :: GPowerGen, gpower
INTEGER :: myzb, myze


TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(Fxrinfo_type),POINTER :: myFxr
TYPE(GamMacXS_TYPE), SAVE :: GamXsMac

INTEGER :: nxy, nFsr, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr
INTEGER :: i, j, ipin, icel, ifsr, ifxr, iz, ig

REAL :: TotGPower, TotGPowerGen, RatioGPower

REAL, POINTER :: xsmackerma(:)

Pin => Core%Pin
CellInfo => Core%CellInfo
nFsr = Core%nCoreFsr
nxy = Core%nxy
! Core%hz, FXR%area
TotGPower= 0.
TotGPowerGen = 0.
DO iz = myzb, myze
    DO ipin = 1, nxy
      FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
      icel = Pin(ipin)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr  
      DO j = 1, nLocalFxr
        ifxr = FxrIdxSt + j -1
        nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
        myFxr => Fxr(ifxr, iz)
        DO i = 1, nFsrInFxr
          ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1  !Global FSR Index
          TotGPower = TotGPower + GPower(ifsr,iz) * Core%hz(iz) * CellInfo(icel)%Vol(Cellinfo(icel)%MapFxr2FsrIdx(i, j))
          TotGPowerGen = TotGPowerGen + GPowerGen(ifsr,iz) * Core%hz(iz) * CellInfo(icel)%Vol(Cellinfo(icel)%MapFxr2FsrIdx(i, j))
        ENDDO
      ENDDO !Fsr Sweep
    ENDDO
ENDDO   
RatioGPower = (TotGPowerGen / TotGPower)
GPower = GPower * RatioGPower

PRINT '(a,es14.6)', 'RatioGPower: ', RatioGPower

#ifdef MPI_ENV
DO iz = 1, Core%nz
  CALL BCAST(gPower(:, iz), nFsr, PE%MPI_RTMASTER_COMM, PE%AxDomList(iz))
ENDDO
#endif

END SUBROUTINE  
  
  
#endif