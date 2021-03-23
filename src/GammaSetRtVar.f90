#include <defines.h>
#ifdef __GAMMA_TRANSPORT
! Subroutines to produce raytracing source and cross section of photon transport calculation
SUBROUTINE SetGamMacXs(core, Fxr, xstr, iz, igg, ngg, lTrCorrection, PE)
! Photo-Atomic Reaction Cross Section Generation Routine
!   Only 0 K  case is treated with Photo-atomic reaction -> No resonance treatment
!   Only Outflow Correction is provided for PHOTON TRANSPORT CROSS SECTION
USE PARAM
USE TYPEDEF,  ONLY : coreinfo_type, Fxrinfo_type, Cell_Type, pin_Type, PE_TYPE
USE OMP_LIB
USE GammaTYPEDEF, ONLY : GamMacXS_TYPE
USE GamXsLib_Mod, ONLY : GamXsBase
USE CNTL,             ONLY : nTracerCntl

IMPLICIT NONE

TYPE(coreinfo_type) :: Core
TYPE(Fxrinfo_type) :: Fxr(:)
TYPE(PE_TYPE) :: PE
REAL, POINTER :: xstr(:)

INTEGER :: iz, igg, ngg
logical :: lGamma, lTrCorrection

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
INTEGER :: nCoreFxr, nCoreFsr, nxy, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr
INTEGER :: icel, ipin, ifxr, ifsr, tid
INTEGER :: i, j
INTEGER :: nofg, norg
REAL :: xsmactr(ngg)

TYPE(GamMacXS_TYPE), SAVE :: GamMacXs(nThreadMax)

Pin => Core%Pin
CellInfo => Core%CellInfo
nCoreFsr = Core%nCoreFsr
nCoreFxr = Core%nCoreFxr
nxy = Core%nxy
tid = 1

!$  call omp_set_dynamic(.FALSE.)
!$  call omp_set_num_threads(PE%nThread)
!$OMP PARALLEL DEFAULT(SHARED)                                                                      &
!$OMP PRIVATE(i, j, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr, icel, ipin, ifxr, ifsr,            &
!$OMP           tid, xsmactr)
!$  tid = omp_get_thread_num()+1
!$OMP DO
DO ipin = 1, nxy
  FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
  icel = Pin(ipin)%Cell(iz)
  nlocalFxr = CellInfo(icel)%nFxr
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j - 1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    !Get Transport XS
      !Obtain Absorption Xs
      CALL GamXsBase(GamMacXs(tid), Fxr(ifxr), igg, igg, ngg, TRUE)
!********** Meaningless in photo-atomic reaction***************************************
!      !Obtain Outgoingg total XS(trasport Corrected)                                 |
!      CALL GamTotScatXs(GamMacXs(tid), Fxr(ifxr), igg, ngg, TRUE)                    |
!**************************************************************************************
      IF(lTrCorrection) THEN
        xsmactr(igg) = GamMacXs(tid)%XsMacA(igg) + GamMacXs(tid)%XsMacSTR(igg)
      ELSE
        xsmactr(igg) = GamMacXs(tid)%XsMacA(igg) + GamMacXs(tid)%XsMacS(igg)
      ENDIF
    !Assgin
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      xstr(ifsr) = xsmactr(igg)
    ENDDO
  ENDDO !End of Fxr Sweep
ENDDO !End of Pin Sweep
!$OMP END DO
!$OMP END PARALLEL
NULLIFY(Pin)
NULLIFY(CellInfo)
END SUBROUTINE


SUBROUTINE SetGamMacXsNM(core, Fxr, xstnm, iz, ngg, lTrCorrection, PE)
! Photo-Atomic Reaction Cross Section Generation Routine
!   Only 0 K  case is treated with Photo-atomic reaction => No resonance treatment
!   Only Outflow Correction is provided for PHOTON TRANSPORT CROSS SECTION
!   For Node Major
USE PARAM
USE TYPEDEF,          ONLY : coreinfo_type, Fxrinfo_type, Cell_Type, pin_Type, PE_TYPE
USE GammaTYPEDEF,     ONLY : GamMacXS_TYPE
USE GamXsLib_Mod,     ONLY : GamXsBase
USE OMP_LIB
USE CNTL,             ONLY : nTracerCntl

IMPLICIT NONE

TYPE(CoreInfo_TYPE) :: Core
TYPE(FxrInfo_TYPE) :: Fxr(:)
TYPE(PE_TYPE) :: PE
REAL, POINTER :: xstnm(:, :)
INTEGER :: iz, ngg
LOGICAL :: lGamma, lTrCorrection

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
INTEGER :: nCoreFxr, nCoreFsr, nxy, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr
INTEGER :: icel, ipin, ifxr, ifsr
INTEGER :: i, j, tid!, igg
INTEGER :: xyb, xye   !--- CNJ Edit : Domain Decomposition + MPI
REAL :: xsmactr(ngg)

TYPE(GamMacXS_TYPE), SAVE :: GamMacXs(nThreadMax)

Pin => Core%Pin
Cellinfo => Core%CellInfo
nCoreFsr = Core%nCoreFsr
nCoreFxr = Core%nCoreFxr
nxy = Core%nxy
xyb = PE%myPinBeg; xye = PE%myPinEnd

!$  call omp_set_dynamic(.FALSE.)
!$  call omp_set_num_threads(PE%nThread) 
!$OMP PARALLEL DEFAULT(SHARED)   &
!$OMP PRIVATE(i, j, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr,icel, ipin, ifxr, ifsr, tid, xsmactr)
!$  tid = omp_get_thread_num() + 1
!$OMP DO
DO ipin = xyb, xye
  FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
  icel = Pin(ipin)%Cell(iz)
  nlocalFxr = CellInfo(icel)%nFxr
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j - 1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    CALL GamXsBase(GamMacXs(tid), Fxr(ifxr), 1, ngg, ngg, TRUE)
!    DO igg = 1, ngg
!      CALL MacTotScatXs(XsMac(tid), Fxr(ifxr), igg, ngg, TRUE)   ! Is this meaningful in photon?
      IF(lTrCorrection) THEN
        xsmactr = GamMacXs(tid)%XsMacA + GamMacXs(tid)%XsMacSTR
      ELSE
        xsmactr = GamMacXs(tid)%XsMacA + GamMacXs(tid)%XsMacS
      ENDIF
!    ENDDO ! gamma group loop end
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      xstnm(:, ifsr) = xsmactr
    ENDDO
  ENDDO ! FXR loop (j) end
ENDDO ! pin loop
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE


! Photon Source Generation Routine
SUBROUTINE SetGamSrc(Core, Fxr, src, phis, gphis, gaxsrc1g, xstr1g, iz, igg, ng, ngg,                  &
                        l3dim, lscat1, PE, GroupInfo)
! Photon Source Generation Routine in Photon Transport Equation
!         Production by Neutron-Isotope Reaction and Photo-Atomic Reaction
USE PARAM
USE TYPEDEF,        ONLY : coreinfo_type, Fxrinfo_type, Cell_Type, pin_Type, PE_TYPE, GROUPINFO_TYPE
USE BasicOperation, ONLY : CP_CA
USE CNTL,           ONLY : nTracerCntl
USE GammaTYPEDEF,   ONLY : GamMacXS_TYPE
USE GamXsLib_Mod,   ONLY : GamProdMatrix, GamScatMatrix
USE OMP_LIB
IMPLICIT NONE
! INPUT VARIABLES
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type) :: Fxr(:)
TYPE(PE_TYPE) :: PE
TYPE(GROUPINFO_TYPE) :: GroupInfo

REAL, POINTER :: src(:), phis(:, :, :), gphis(:, :, :), gAxSrc1g(:), xstr1g(:)
INTEGER :: igg, ng, iz, ifsr, ifxr, ngg
LOGICAL :: lscat1, l3dim

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(GamMacXS_TYPE), SAVE :: GamMacXS(nThreadMax)
REAL, POINTER :: xsmacs(:,:)
REAL, POINTER :: prodmat(:,:)

INTEGER :: nxy, nCoreFsr, nCoreFxr, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr
INTEGER :: ipin, icel, ig, tid
INTEGER :: i, j
INTEGER :: iso, niso

!ALLOCATE(GamMacXS(PE%nThread))

Pin => Core%Pin
CellInfo => Core%CellInfo
nCoreFsr = Core%nCoreFsr
nCoreFxr = Core%nCoreFxr
nxy = Core%nxy
!lNegSrcFix = FALSE   <=== what is this??

src = zero

tid = 1

!$ call omp_set_dynamic(.FALSE.)
!$ call omp_set_num_threads(PE%nThread)
!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(i, j, ifsr, ifxr, ipin, icel, ig, tid, FsrIdxSt, FxrIdxSt,                         &
!$OMP              nlocalFxr, nFsrInFxr, xsmacs, prodmat)
!$  tid = omp_get_thread_num()+1
!$OMP DO
DO ipin = 1, nxy ! Do loop begin
  FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
  icel = Pin(ipin)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr

!   Neutron Induced Photon Source (resonance dependence)
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j -1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    CALL GamProdMatrix(GamMacXS(tid), Fxr(ifxr), igg, igg, ng, ngg, GroupInfo, FALSE)
    prodmat => GamMacXS(tid)%GProdTot
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      DO ig = 1, ng
          src(ifsr) = src(ifsr) + prodmat(ig, igg) * phis(ifsr, iz, ig)
      ENDDO
    ENDDO !Fsr Sweep
  ENDDO ! Local Fxr loop

!   Photo-atomic Reaction Source (scattering source)
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j -1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)

    CALL GamScatMatrix(GamMacXS(tid), Fxr(ifxr), igg, igg, ngg, lscat1,.FALSE.)
    XsMacS => GamMacXS(tid)%XsMacSM

    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      DO ig = 1, ngg
        src(ifsr) = src(ifsr) + xsmacs(ig, igg) * gphis(ifsr, ig, iz)
      ENDDO
    ENDDO !Fsr Sweep
  ENDDO !End of Fxr Sweep
!  Axial Source Substitue
  IF (l3dim) THEN
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j - 1
      nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
      DO i = 1, nFsrInFxr
        ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
#ifndef LkgSplit
        src(ifsr) = src(ifsr) - gAxSrc1g(ipin)
#else
        IF(gAxSrc1g(ipin) .LT. 0 .AND. .NOT. Fxr(ifxr)%lvoid) THEN
          src(ifsr) = src(ifsr) - gAxSrc1g(ipin)
        ENDIF
#endif
      END DO
    END DO
  END IF
  DO j = 1, CellInfo(icel)%nFsr
    ifsr = FsrIdxSt + j - 1
    src(ifsr) = src(ifsr)/xstr1g(ifsr)
  END DO

  NULLIFY(XsMacS, prodmat)
END DO   ! pin loop
!$OMP END DO
!$OMP END PARALLEL
!DEALLOCATE(GamMacXS)
NULLIFY(Pin)
NULLIFY(CellInfo)
END SUBROUTINE

SUBROUTINE SetGamSrcNM(Core, Fxr, srcNM, phisNM, gphisNM, gAxSrc, xstnm, iz, igb, ige,               &
                       ng, ngg, l3dim, lscat1, PE)
! Photon Source Generation Routine in Photon Transport Equation
!      FOR NODE MAJOR
USE PARAM
USE TYPEDEF,          ONLY : coreinfo_type, Fxrinfo_type, Cell_Type, pin_Type, PE_TYPE
USE BasicOperation,   ONLY : CP_CA
USE CNTL,             ONLY : nTracerCntl
USE GammaTYPEDEF,     ONLY : GamMacXS_TYPE
USE Core_mod,         ONLY : GroupInfo
USE GamXsLib_Mod,     ONLY : GamProdMatrix, GamScatMatrix
USE OMP_LIB
IMPLICIT NONE
! INPUT VARIABLES
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type) :: Fxr(:)
REAL, POINTER :: srcNM(:,:), phisNM(:, :), gphisNM(:, :), gAxSrc(:, :, :), xstnm(:, :)
REAL :: eigv
INTEGER :: iz, igb, ige, ng, ngg, ifsr, ifxr, fsridx
LOGICAL :: l3dim, lscat1
TYPE(PE_TYPE) :: PE

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(GamMacXS_TYPE), SAVE :: GamMacXS(nThreadMax)
REAL, POINTER :: prodmat(:,:)
REAL, POINTER :: xsmacs(:,:)

INTEGER :: nxy, nCoreFsr, nCoreFxr, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr
INTEGER :: ipin, icel, igg, tid, ig, igg2
INTEGER :: i, j
INTEGER :: iso, niso

Pin => Core%Pin
CellInfo => Core%CellInfo
nCoreFsr = Core%nCoreFsr
nCoreFxr = Core%nCoreFxr
nxy = Core%nxy

SrcNM = zero

tid = 1
!$  call omp_set_dynamic(.FALSE.)
!$  call omp_set_num_threads(PE%nThread)
!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(i, j, ifsr, ifxr, ipin, icel, ig, igg, igg2, tid, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr, xsmacs, prodmat)
tid = omp_get_thread_num() + 1
!$OMP DO
DO ipin = 1, nxy
  FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
  icel = Pin(ipin)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
! Neutron Induced Photon Source (resonance dependence)
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j - 1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    CALL GamProdMatrix(GamMacXS(tid), Fxr(ifxr), igb, ige, ng, ngg, GroupInfo, FALSE)
    prodmat => GamMacXS(tid)%GProdTot
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      ! Neutron Induced Source  (ig -> igg)
      DO igg = igb, ige
        DO ig = 1, ng
          srcnm(igg, ifsr) = srcnm(igg, ifsr) + prodmat(ig, igg) * phisNM(ig, ifsr)
        END DO ! Departure Neutron Group loop
      ENDDO ! Arrival Photon Group loop
    ENDDO ! Fsr loop
  ENDDO ! Local Fxr loop (j)
! Photo-atomic Reaction Source (scattering source)
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j -1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    CALL GamScatMatrix(GamMacXS(tid), Fxr(ifxr), igb, ige, ngg, lscat1,.FALSE.)
    XsMacS => GamMacXS(tid)%XsMacSM

    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      ! Photon Induced Source  (igg2 -> igg)
      DO igg = igb, ige
        DO igg2 = 1, ngg
          srcnm(igg, ifsr) = srcnm(igg, ifsr) + xsmacs(igg2, igg) * gphisNM(igg2, ifsr)
        ENDDO  ! Departure Photon Group loop
      ENDDO  ! Arrival Photon Group loop
    ENDDO  ! Fsr loop
  ENDDO  ! Local Fxr loop (j)
!  Axial Source Substitue
  IF(l3dim) THEN
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j -1
      nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
      DO i = 1, nFsrInFxr
        ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
        DO igg = igb, ige
#ifndef LkgSplit
          srcnm(igg, ifsr) = srcnm(igg, ifsr) - gAxSrc(ipin, igg, iz)
#else
          IF(gAxSrc(ipin, igg, iz) .LT. 0 .AND. .NOT. Fxr(ifxr)%lvoid) THEN
            srcnm(igg, ifsr) = srcnm(igg, ifsr) - gAxSrc(ipin, igg, iz)
          ENDIF
#endif
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  DO j = 1, CellInfo(icel)%nFsr
    ifsr = FsrIdxSt + j - 1
    srcnm(igb:ige, ifsr) = srcnm(igb:ige, ifsr) / xstnm(igb:ige, ifsr)
  ENDDO
  NULLIFY(Xsmacs, prodmat)
ENDDO ! pin loop
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE SetGamP1Src(Core, Fxr, srcm, phim, xstr1g, iz, igg, ngg, lscat1, ScatOd, PE)
USE PARAM
USE TYPEDEF,        ONLY : coreinfo_type,          Fxrinfo_type,          Cell_Type,                &
                           pin_Type,               PE_Type
USE BasicOperation, ONLY : CP_CA
USE GammaTYPEDEF,   ONLY : GamMacXS_TYPE
USE GamXsLib_Mod,   ONLY : GamP1XsScatMatrix,     GamP2XsScatMatrix,     GamP3XsScatMatrix
USE OMP_LIB
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type) :: Fxr(:)
REAL, POINTER :: srcm(:, :)
REAL, POINTER :: phim(:, :, :, :)
REAL, POINTER :: xstr1g(:)
INTEGER :: myzb, myze, igg, ngg, iz, ifsr, ifxr
LOGICAL :: lscat1
INTEGER :: ScatOd
TYPE(PE_Type) :: PE

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(GamMacXS_TYPE), SAVE :: GamMacXS(nThreadMax)
INTEGER :: nxy, nCoreFsr, nCoreFxr, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr, nchi
INTEGER :: ipin, icel, ifsrlocal, itype, igg2, tid
INTEGER :: i, j

REAL, POINTER :: XsMacP1sm(:,:), XsMacP2sm(:,:), XsMacP3sm(:,:)

LOGICAL :: lscat1sum


Pin => Core%Pin
CellInfo => Core%CellInfo
nCoreFsr = Core%nCoreFsr
nCoreFxr = Core%nCoreFxr
nxy = Core%nxy

lscat1sum = lscat1

IF (ScatOd .EQ. 1) THEN
  CALL CP_CA(srcm(:, :), zero, 2, nCoreFsr)
ELSEIF (ScatOd .EQ. 2) THEN
  CALL CP_CA(srcm(:, :), zero, 5, nCoreFsr)
ELSEIF (ScatOd .EQ. 3) THEN
  CALL CP_CA(srcm(:, :), zero, 9, nCoreFsr)
ENDIF

tid = 1

!Scattering Source Update
!$  call omp_set_dynamic(.FALSE.)
!$  call omp_set_num_threads(PE%nThread)
!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(i, j, ifsr, ifxr, ipin, icel, ifsrlocal, itype, igg2, tid, FsrIdxSt,                  &
!$OMP                   FxrIdxSt, nlocalFxr, nFsrInFxr, XsMacP1Sm, XsMacP2Sm, XsMacP3Sm)
!$  tid = omp_get_thread_num()+1
!$OMP DO
DO ipin = 1, nxy
  FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
  icel = Pin(ipin)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j -1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    CALL GamP1XsScatMatrix(GamMacXS(tid), Fxr(ifxr), igg, igg, ngg)
    IF(ScatOd .GE. 2) CALL GamP2XsScatMatrix(GamMacXS(tid), Fxr(ifxr), igg, igg, ngg)
    IF(ScatOd .EQ. 3) CALL GamP3XsScatMatrix(GamMacXS(tid), Fxr(ifxr), igg, igg, ngg)
    XsMacP1Sm => GamMacXS(tid)%MacGSM1
    XsMacP2Sm => GamMacXS(tid)%MacGSM2
    XsMacP3Sm => GamMacXS(tid)%MacGSM3
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      IF(ScatOd .EQ. 1) THEN
        DO igg2 = 1, ngg
           srcm(1, ifsr) = srcm(1, ifsr) + XsMacP1Sm(igg2, igg) * phim(1, ifsr, igg2, iz)
           srcm(2, ifsr) = srcm(2, ifsr) + XsMacP1Sm(igg2, igg) * phim(2, ifsr, igg2, iz)
        ENDDO
      ELSEIF(ScatOd .EQ. 2) THEN
        DO igg2 = 1, ngg
          srcm(1:2, ifsr) = srcm(1:2, ifsr) + XsMacP1Sm(igg2, igg) * phim(1:2, ifsr, igg2, iz)
          srcm(3:5, ifsr) = srcm(3:5, ifsr) + XsMacP2Sm(igg2, igg) * phim(3:5, ifsr, igg2, iz)
        ENDDO
      ELSEIF(ScatOd .EQ. 3) THEN
        DO igg2 = 1, ngg
          srcm(1:2, ifsr) = srcm(1:2, ifsr) + XsMacP1Sm(igg2, igg) * phim(1:2, ifsr, igg2, iz)
          srcm(3:5, ifsr) = srcm(3:5, ifsr) + XsMacP2Sm(igg2, igg) * phim(3:5, ifsr, igg2, iz)
          srcm(6:9, ifsr) = srcm(6:9, ifsr) + XsMacP3Sm(igg2, igg) * phim(6:9, ifsr, igg2, iz)
        ENDDO
      ENDIF

    ENDDO !Fsr Sweep
  ENDDO !End of Fxr Sweep
!  CALL ReturnXsMacDat(XsMac)   !Memory leak problem !! 17/01/09   big memory but stable
ENDDO !End of Pin
!$OMP END DO
!
!IF(.NOT. lxsLib) DEALLOCATE(XsMacP1sm) ! modified because of crash! in benchmark XS
!IF(.NOT. lxsLib .AND. ScatOd .GE. 2) DEALLOCATE(XsMacP2sm)
!IF(.NOT. lxsLib .AND. ScatOd .EQ. 3) DEALLOCATE(XsMacP3sm)
!$OMP END PARALLEL
!Axail Source Contribution
IF(ScatOd .EQ. 1) THEN
  DO ipin = 1, nxy
    FsrIdxSt = Pin(ipin)%FsrIdxSt; icel = Pin(ipin)%Cell(iz);
    DO j = 1, CellInfo(icel)%nFsr
      ifsr = FsrIdxSt + j - 1
      srcm(1:2, ifsr) = 3._8 * srcm(1:2, ifsr) / xstr1g(ifsr)
    ENDDO
  ENDDO
ELSEIF(ScatOd .EQ. 2) THEN
  DO ipin = 1, nxy
    FsrIdxSt = Pin(ipin)%FsrIdxSt; icel = Pin(ipin)%Cell(iz);
    DO j = 1, CellInfo(icel)%nFsr
      ifsr = FsrIdxSt + j - 1
      srcm(1:2, ifsr) = 3._8 * srcm(1:2, ifsr) / xstr1g(ifsr)
      srcm(3:5, ifsr) = 5._8 * srcm(3:5, ifsr) / xstr1g(ifsr)
    ENDDO
  ENDDO
ELSEIF(ScatOd .EQ. 3) THEN
  DO ipin = 1, nxy
    FsrIdxSt = Pin(ipin)%FsrIdxSt; icel = Pin(ipin)%Cell(iz);
    DO j = 1, CellInfo(icel)%nFsr
      ifsr = FsrIdxSt + j - 1
      srcm(1:2, ifsr) = 3._8 * srcm(1:2, ifsr) / xstr1g(ifsr)
      srcm(3:5, ifsr) = 5._8 * srcm(3:5, ifsr) / xstr1g(ifsr)
      srcm(6:9, ifsr) = 7._8 * srcm(6:9, ifsr) / xstr1g(ifsr)
    ENDDO
  ENDDO
ENDIF
END SUBROUTINE

SUBROUTINE GamPseudoAbsorption(Core, Fxr, AxPXS, xstr1g, iz, l3dim)
USE PARAM
USE TYPEDEF,      ONLY : coreinfo_type,          Fxrinfo_type,          Cell_Type,      &
                         pin_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type) :: Fxr(:)
REAL, POINTER :: AxPXS(:), xstr1g(:)
INTEGER :: iz
LOGICAL :: l3dim

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
INTEGER :: nxy
INTEGER :: FsrIdxSt, FxrIdxSt, nLocalFxr, nFsrInFxr
INTEGER :: i, j, ipin, icel, ifsr, ifxr

IF(.NOT. l3dim) RETURN

Pin => Core%Pin
CellInfo => Core%CellInfo
nxy = Core%nxy

DO ipin = 1, nxy
  FsrIdxSt = Pin(ipin)%FsrIdxSt; icel = Pin(ipin)%Cell(iz);
  FxrIdxSt = Pin(ipin)%FxrIdxSt; nLocalFxr = CellInfo(icel)%nFxr
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j - 1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      IF(.NOT. Fxr(ifxr)%lVoid) xstr1g(ifsr) = xstr1g(ifsr)+AxPXS(ipin)
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE GamPseudoAbsorptionNM(Core, Fxr, AxPXS, xstnm, iz, ng, l3dim)
USE PARAM
USE TYPEDEF,      ONLY : coreinfo_type,          Fxrinfo_type,          Cell_Type,      &
                         pin_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type) :: Fxr(:)
REAL, POINTER :: AxPXS(:, :, :), xstnm(:, :)
INTEGER :: iz, ng
LOGICAL :: l3dim

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
INTEGER :: nxy
INTEGER :: FsrIdxSt, FxrIdxSt, nLocalFxr, nFsrInFxr
INTEGER :: i, j, ipin, icel, ifsr, ifxr, ig

IF(.NOT. l3dim) RETURN

Pin => Core%Pin
CellInfo => Core%CellInfo
nxy = Core%nxy

DO ipin = 1, nxy
  FsrIdxSt = Pin(ipin)%FsrIdxSt; icel = Pin(ipin)%Cell(iz);
  FxrIdxSt = Pin(ipin)%FxrIdxSt; nLocalFxr = CellInfo(icel)%nFxr
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j - 1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      DO ig = 1, ng
        IF(.NOT. Fxr(ifxr)%lVoid) xstnm(ig, ifsr) = xstnm(ig, ifsr) + AxPXS(ipin, ig, iz)
      ENDDO
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE
#endif
