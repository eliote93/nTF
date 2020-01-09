#include <defines.h>    
!--- CNJ Edit : Auxiliary Routines for Gamma MOC Calculation
#ifdef __GAMMA_TRANSPORT    
FUNCTION GammaMocResidual(Core, FmInfo, GroupInfo, PE, nTracerCntl)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,       FmInfo_Type,        GroupInfo_Type,     &
                           Cell_Type,           Pin_Type,           FxrInfo_Type,       &
                           PE_Type
USE GammaCore_mod,  ONLY : gphis,               gJout,              GamGroupInfo
USE GamMOC_MOD,     ONLY : SetGamSrc,           SetGamSrcNM,        SetGamMacXs,        &
                           SetGamMacXsNM,       gxst1g,             gsrc1g,             &
                           gphisnm,             gxstnm,             gsrcnm
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
REAL, POINTER :: AxSrc(:, :, :), AxSrc1g(:)
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
nFsr = Core%nCoreFsr; nxy = Core%nxy; nbd = 4; ng = GroupInfo%ng; ngg = GamGroupInfo%ngg
GammaMocResidual = 0; SrcSum = 0

IF (lscat1) THEN

DO ig = 1, ngg
  DO iz = myzb, myze
    IF (.NOT. Core%lFuelPlane(iz)) CYCLE
    CALL SetGamMacXs(Core, Fxr(:, iz), gxst1g, iz, ig, ngg, lTrCorrection, PE)
    CALL SetGamSrc(Core, Fxr(:, iz), gsrc1g, phis, gphis, AxSrc1g, gxst1g, iz, ig, ng, ngg,                   &
                   l3dim, lscat1, PE, GamGroupInfo)
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
  CALL SetGamSrcNM(Core, Fxr(:, iz), gsrcnm, phisnm, gphisnm, AxSrc, gxstnm, iz,                               &
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
  
SUBROUTINE NeutronLocalQUpdate(Core, Fxr, phis, localpower, myzb, myze, ng, PE)
USE PARAM
USE TYPEDEF,        ONLY : coreinfo_type,       Fxrinfo_type,       Cell_Type,     pin_Type, &
                           PE_TYPE
USE CNTL,           ONLY : nTracerCntl
USE CORE_MOD,       ONLY : GroupInfo
USE GamXsLib_Mod,   ONLY : GetLocalQMat
#ifdef MPI_ENV      
USE MPIComm_Mod,    ONLY : BCAST
#endif              
USE GammaTYPEDEF,   ONLY : GamMacXS_TYPE
USE GammaCore_mod,  ONLY : GamGroupInfo
IMPLICIT NONE

TYPE(coreinfo_type) :: CORE
TYPE(Fxrinfo_type),POINTER :: Fxr(:, :)
TYPE(PE_Type) :: PE
REAL, POINTER :: phis(:, :, :)
REAL, POINTER :: localpower(:, :)
INTEGER :: myzb, myze, ng, ngg

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(Fxrinfo_type),POINTER :: myFxr

TYPE(GamMacXS_TYPE) :: XsMac

INTEGER :: nxy, nCoreFsr, nCoreFxr, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr
INTEGER :: ipin, icel, ifsr, ifxr, iz, ig
INTEGER :: i, j

REAL, POINTER :: MacLocalQ(:)
LOGICAL :: lfis

Pin => Core%Pin
CellInfo => Core%CellInfo; nCoreFsr = Core%nCoreFsr
nCoreFxr = Core%nCoreFxr; nxy = Core%nxy

ngg = GamGroupInfo%ngg

localpower = 0.

DO iz = myzb, myze
  DO ipin = 1, nxy
    FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
    icel = Pin(ipin)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr  
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j -1
      nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)    
      myFxr => Fxr(ifxr, iz)
      lfis = myFxr%lres .and. Core%lFuelPlane(iz)
      CALL GetLocalQMat(XSMac, myFxr, 1, ng, ng, ngg, lfis)
      MacLocalQ => XsMac%LocalQ
      DO i = 1, nFsrInFxr
        ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1  !Global FSR Index
        DO ig = 1, ng
          localpower(ifsr, iz) = localpower(ifsr, iz) + MacLocalQ(ig) * phis(ifsr, iz, ig)
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
#endif