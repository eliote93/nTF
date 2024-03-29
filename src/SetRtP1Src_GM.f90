SUBROUTINE SetRtP1SrcGM(Core, Fxr, srcm1g, phimNg, xst1g, iz, ig, ng, GroupInfo, l3dim, lxslib, lscat1, ScatOd, PE)

USE OMP_LIB
USE PARAM,        ONLY : ZERO, TRUE, FALSE, ONE, FORWARD, BACKWARD
USE TYPEDEF,      ONLY : coreinfo_type, Fxrinfo_type, Cell_Type, pin_Type, GroupInfo_Type, XsMac_Type, PE_Type
USE BenchXs,      ONLY : GetChiBen, xssben, xssm1ben, xssm2ben, xssm3ben, xssm1DynBen, xssm2DynBen, xssm3DynBen
USE MacXsLib_mod, ONLY : MacP1XsScatMatrix, MacP2XsScatMatrix, MacP3XsScatMatrix
USE XsUtil_mod,   ONLY : GetXsMacDat, ReturnXsMacDat
USE TRAN_MOD,     ONLY : TranInfo, TranCntl
USE MOC_MOD,      ONLY : mwt

IMPLICIT NONE

TYPE (CoreInfo_Type)  :: Core
TYPE (GroupInfo_Type) :: GroupInfo
TYPE (PE_Type)        :: PE

TYPE(FxrInfo_Type), DIMENSION(:) :: Fxr

REAL, POINTER, DIMENSION(:)     :: xst1g
REAL, POINTER, DIMENSION(:,:)   :: srcm1g
REAL, POINTER, DIMENSION(:,:,:) :: phimNg

INTEGER :: ig, ng, iz, ScatOd
LOGICAL :: lxslib, lscat1, l3dim
! ----------------------------------------------------
TYPE (XsMac_Type), POINTER :: XsMac

TYPE (Pin_Type),  POINTER, DIMENSION(:) :: Pin
TYPE (Cell_Type), POINTER, DIMENSION(:) :: CellInfo

INTEGER :: nxy, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr
INTEGER :: ifsr, jfsr, ifxr, jfxr, ipin, icel, itype, jg, tid

REAL :: xstinv

REAL, POINTER, DIMENSION(:,:) :: XsMacP1sm, XsMacP2sm, XsMacP3sm
! ----------------------------------------------------

Pin      => Core%Pin
CellInfo => Core%CellInfo
nxy       = Core%nxy
! ----------------------------------------------------
IF (.NOT. lxsLib) THEN
  ALLOCATE (XsMacP1sm (ng, ng))
  
  IF (ScatOd .GE. 2) ALLOCATE (XsMacP2sm (ng, ng))
  IF (ScatOd .EQ. 3) ALLOCATE (XsMacP3sm (ng, ng))
  
  XsMacP1sm = ZERO
  
  IF (ScatOd .GE. 2) XsMacP2sm = ZERO
  IF (ScatOd .EQ. 3) XsMacP3sm = ZERO
END IF

srcm1g = ZERO
tid  = 1

!$ call omp_set_dynamic(FALSE)
!$ call omp_set_num_threads(PE%nThread)
!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(XsMac, ifsr, ifxr, jfsr, jfxr, ipin, icel, itype, jg, tid, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr, XsMacP1Sm, XsMacP2Sm, XsMacP3Sm)
!$  tid = omp_get_thread_num() + 1
!$OMP DO
DO ipin = 1, nxy
  CALL GetXsMacDat(XsMac, ng, TRUE)
  
  FsrIdxSt  = Pin(ipin)%FsrIdxSt
  FxrIdxSt  = Pin(ipin)%FxrIdxSt
  icel      = Pin(ipin)%Cell(iz)
  nlocalFxr = CellInfo(icel)%nFxr
  
  ! SET : XS
  DO ifxr = 1, nLocalFxr
    jfxr      = FxrIdxSt + ifxr -1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(ifxr)
    
    IF (lXsLib) Then
      CALL MacP1XsScatMatrix(XsMac, Fxr(jfxr), ig, ig, ng, GroupInfo)
      
      IF (ScatOd .GE. 2) CALL MacP2XsScatMatrix(XsMac, Fxr(jfxr), ig, ig, ng, GroupInfo)
      IF (ScatOd .EQ. 3) CALL MacP3XsScatMatrix(XsMac, Fxr(jfxr), ig, ig, ng, GroupInfo)
      
      XsMacP1Sm => XsMac%XsMacP1Sm
      XsMacP2Sm => XsMac%XsMacP2Sm
      XsMacP3Sm => XsMac%XsMacP3Sm
    ELSE
      itype = Fxr(jfxr)%imix
      
      IF (TranCntl%lDynamicBen) THEN
        CALL XsSm1DynBen(itype, TranInfo%fuelTemp(ipin, iz), 1, ng, 1, ng, XsMacP1Sm)
        
        IF (ScatOd .GE. 2) CALL XsSm2DynBen(itype, TranInfo%fuelTemp(ipin, iz), 1, ng, 1, ng, XsMacP2Sm) ! OPTIMZIE
        IF (ScatOd .EQ. 3) CALL XsSm3DynBen(itype, TranInfo%fuelTemp(ipin, iz), 1, ng, 1, ng, XsMacP3Sm)
      ELSE
        CALL XsSm1Ben(itype, 1, ng, 1, ng, XsMacP1Sm)
        
        IF (ScatOd .GE. 2) CALL XsSm2Ben(itype, 1, ng, 1, ng, XsMacP2Sm) ! OPTIMZIE
        IF (ScatOd .EQ. 3) CALL XsSm3Ben(itype, 1, ng, 1, ng, XsMacP3Sm)
      END IF
    END IF
    
    ! CAL : src m
    DO ifsr = 1, nFsrInFxr
      jfsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(ifsr, ifxr) - 1
      
      SELECT CASE (ScatOd)
      CASE (1)
        DO jg = 1, ng
          srcm1g(1:2, jfsr) = srcm1g(1:2, jfsr) + XsMacP1Sm(jg, ig) * phimNg(1:2, jfsr, jg)
        END DO
      CASE (2)
        DO jg = 1, ng
          srcm1g(1:2, jfsr) = srcm1g(1:2, jfsr) + XsMacP1Sm(jg, ig) * phimNg(1:2, jfsr, jg)
          srcm1g(3:5, jfsr) = srcm1g(3:5, jfsr) + XsMacP2Sm(jg, ig) * phimNg(3:5, jfsr, jg)
        END DO
      CASE (3)
        DO jg = 1, ng
          srcm1g(1:2, jfsr) = srcm1g(1:2, jfsr) + XsMacP1Sm(jg, ig) * phimNg(1:2, jfsr, jg)
          srcm1g(3:5, jfsr) = srcm1g(3:5, jfsr) + XsMacP2Sm(jg, ig) * phimNg(3:5, jfsr, jg)
          srcm1g(6:9, jfsr) = srcm1g(6:9, jfsr) + XsMacP3Sm(jg, ig) * phimNg(6:9, jfsr, jg)
        END DO
      END SELECT
    END DO
  END DO
  
  CALL ReturnXsMacDat(XsMac)
END DO
!$OMP END DO
!$OMP END PARALLEL

IF (.NOT. lxsLib) THEN
  DEALLOCATE (XsMacP1sm)
  
  IF (ScatOd .GE. 2) DEALLOCATE (XsMacP2sm)
  IF (ScatOd .EQ. 3) DEALLOCATE (XsMacP3sm)
ELSE
  NULLIFY (XsMacP1Sm)
  
  IF (ScatOd .GE. 2) NULLIFY (XsMacP2sm)
  IF (ScatOd .EQ. 3) NULLIFY (XsMacP3sm)
END IF
! ----------------------------------------------------
!$OMP PARALLEL PRIVATE(ipin, FsrIdxSt, icel, ifxr, jfsr, xstinv)
!$OMP DO
DO ipin = 1, nxy
  FsrIdxSt = Pin(ipin)%FsrIdxSt
  icel     = Pin(ipin)%Cell(iz)
  
  DO ifxr = 1, CellInfo(icel)%nFsr
    jfsr   = FsrIdxSt + ifxr - 1
    xstinv = ONE / xst1g(jfsr)
    
    srcm1g(1:2, jfsr) = srcm1g(1:2, jfsr) * xstinv
    
    IF (ScatOd .GE. 2) srcm1g(3:5, jfsr) = srcm1g(3:5, jfsr) * xstinv
    IF (ScatOd .EQ. 3) srcm1g(6:9, jfsr) = srcm1g(6:9, jfsr) * xstinv
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

NULLIFY (Pin)
NULLIFY (CellInfo)
! ----------------------------------------------------

END SUBROUTINE SetRtP1SrcGM