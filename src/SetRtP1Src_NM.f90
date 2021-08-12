SUBROUTINE SetRtP1srcNM(Core, Fxr, srcmNg, phimNg, xstNg, iz, gb, ge, ng, GroupInfo, lxslib, ScatOd, PE, Offset)

USE OMP_LIB
USE PARAM,        ONLY : ZERO, ONE, TRUE, FALSE
USE TYPEDEF,      ONLY : coreinfo_type, Fxrinfo_type, Cell_Type, pin_Type, GroupInfo_Type, XsMac_Type, PE_Type
USE BenchXs,      ONLY : GetChiBen, xssben, xssm1ben, xssm2ben, xssm3ben, xssm1DynBen, xssm2DynBen, xssm3DynBen
USE MacXsLib_mod, ONLY : MacP1XsScatMatrix, MacP2XsScatMatrix, MacP3XsScatMatrix
USE XsUtil_mod,   ONLY : GetXsMacDat, ReturnXsMacDat
USE TRAN_MOD,     ONLY : TranInfo, TranCntl

IMPLICIT NONE

TYPE (CoreInfo_Type)  :: Core
TYPE (GroupInfo_Type) :: GroupInfo
TYPE (PE_Type)        :: PE

TYPE (FxrInfo_Type), DIMENSION(:) :: Fxr

REAL, POINTER, DIMENSION(:,:)   :: xstNg
REAL, POINTER, DIMENSION(:,:,:) :: srcmNg, phimNg

INTEGER :: gb, ge, ng, iz, ScatOd
LOGICAL :: lxslib, lscat1, l3dim
INTEGER, OPTIONAL :: Offset
! ----------------------------------------------------
TYPE(XsMac_Type), POINTER :: XsMac

TYPE(Pin_Type),  POINTER, DIMENSION(:) :: Pin
TYPE(Cell_Type), POINTER, DIMENSION(:) :: CellInfo

INTEGER :: nxy, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr, xyb, xye, ioff
INTEGER :: ifsr, jfsr, ifxr, jfxr, ipin, icel, itype, ig, jg

REAL :: xstinv

REAL, DIMENSION(9)  :: scatsrc
REAL, POINTER, DIMENSION(:,:) :: xsmacP1sm, xsmacP2sm, xsmacP3sm
! ----------------------------------------------------

Pin      => Core%Pin
CellInfo => Core%CellInfo
nxy       = Core%nxy

xyb = PE%myPinBeg ! Domain Decomposition + MPI
xye = PE%myPinEnd

ioff = 0
IF (PRESENT(Offset)) ioff = Offset
! ----------------------------------------------------
IF (.NOT. lxsLib) THEN
  ALLOCATE (XsMacP1sm (ng, ng))
  
  IF (ScatOd .GE. 2) ALLOCATE (XsMacP2sm (ng, ng))
  IF (ScatOd .EQ. 3) ALLOCATE (XsMacP3sm (ng, ng))
  
  XsMacP1sm = ZERO
  
  IF (ScatOd .GE. 2) XsMacP2sm = ZERO
  IF (ScatOd .EQ. 3) XsMacP3sm = ZERO
END IF

srcmNg = ZERO

!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(XsMac, ifsr, ifxr, jfsr, jfxr, icel, itype, scatSrc, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr, XsMacP1Sm, XsMacP2Sm, XsMacP3Sm)
!$OMP DO
DO ipin = xyb, xye
  CALL GetXsMacDat(XsMac, ng, TRUE)
  
  FsrIdxSt  = Pin(ipin)%FsrIdxSt
  FxrIdxSt  = Pin(ipin)%FxrIdxSt
  icel      = Pin(ipin)%Cell(iz)
  nlocalFxr = CellInfo(icel)%nFxr
  
  DO ifxr = 1, nLocalFxr
    jfxr      = FxrIdxSt + ifxr - 1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(ifxr)
    
    IF (lXsLib) Then
      CALL MacP1XsScatMatrix(XsMac, Fxr(jfxr), gb, ge, ng, GroupInfo)
      
      IF (ScatOd .GE. 2) CALL MacP2XsScatMatrix(XsMac, Fxr(jfxr), gb, ge, ng, GroupInfo)
      IF (ScatOd .EQ. 3) CALL MacP3XsScatMatrix(XsMac, Fxr(jfxr), gb, ge, ng, GroupInfo)
      
      XsMacP1Sm => XsMac%XsMacP1Sm
      XsMacP2Sm => XsMac%XsMacP2Sm
      XsMacP3Sm => XsMac%XsMacP3Sm
    ELSE
      itype = Fxr(jfxr)%imix
      
      IF(TranCntl%lDynamicBen) THEN
        CALL XsSm1DynBen(itype, TranInfo%fuelTemp(ipin, iz), 1, ng, 1, ng, XsMacP1Sm)
        
        IF (ScatOd .GE. 2) CALL XsSm2DynBen(itype, TranInfo%fuelTemp(ipin, iz), 1, ng, 1, ng, XsMacP2Sm) ! OPTIMZIE
        IF (ScatOd .EQ. 3) CALL XsSm3DynBen(itype, TranInfo%fuelTemp(ipin, iz), 1, ng, 1, ng, XsMacP3Sm)
      ELSE
        CALL XsSm1Ben(itype, 1, ng, 1, ng, XsMacP1Sm)
        
        IF (ScatOd .GE. 2) CALL XsSm2Ben(itype, 1, ng, 1, ng, XsMacP2Sm) ! OPTIMZIE
        IF (ScatOd .EQ. 3) CALL XsSm3Ben(itype, 1, ng, 1, ng, XsMacP3Sm)
      END IF
    END IF
    
    DO ifsr = 1, nFsrInFxr
      jfsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(ifsr, ifxr) - 1
      
      DO ig = gb, ge
        scatSrc = ZERO
        
        SELECT CASE (ScatOd)
        CASE (1)
          DO jg = 1, ng
            scatSrc(1:2) = scatSrc(1:2) + XsMacP1Sm(jg, ig) * phimNg(1:2, jg, jfsr)
          END DO
          
          srcmNg(1:2, ig - ioff, jfsr) = scatSrc(1:2)
        CASE (2)
          DO jg = 1, ng
            scatSrc(1:2) = scatSrc(1:2) + XsMacP1Sm(jg, ig) * phimNg(1:2, jg, jfsr)
            scatSrc(3:5) = scatSrc(3:5) + XsMacP2Sm(jg, ig) * phimNg(3:5, jg, jfsr)
          ENDDO
          
          srcmNg(1:5, ig - ioff, jfsr) = scatSrc(1:5)
        CASE (3)
          DO jg = 1, ng
            scatSrc(1:2) = scatSrc(1:2) + XsMacP1Sm(jg, ig) * phimNg(1:2, jg, jfsr)
            scatSrc(3:5) = scatSrc(3:5) + XsMacP2Sm(jg, ig) * phimNg(3:5, jg, jfsr)
            scatSrc(6:9) = scatSrc(6:9) + XsMacP3Sm(jg, ig) * phimNg(6:9, jg, jfsr)
          ENDDO
          
          srcmNg(1:9, ig - ioff, jfsr) = scatSrc(1:9)
        END SELECT
      END DO
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
  NULLIFY (XsMacP1sm)
  
  IF (ScatOd .GE. 2) NULLIFY (XsMacP2sm)
  IF (ScatOd .EQ. 3) NULLIFY (XsMacP3sm)
END IF
! ----------------------------------------------------
!$OMP PARALLEL PRIVATE(ipin, FsrIdxSt, icel, ifxr, jfsr, ig, xstinv)
!$OMP DO
DO ipin = xyb, xye
  FsrIdxSt = Pin(ipin)%FsrIdxSt
  icel     = Pin(ipin)%Cell(iz)
  
  DO ifxr = 1, CellInfo(icel)%nFsr
    jfsr = FsrIdxSt + ifxr - 1
    
    DO ig = gb, ge
      xstinv = ONE / xstNg(ig, jfsr)
      
      srcmNg(1:2, ig - ioff, jfsr) = 3._8 * srcmNg(1:2, ig - ioff, jfsr) * xstinv
      
      IF (ScatOd .GE. 2) srcmNg(3:5, ig - ioff, jfsr) = 5._8 * srcmNg(3:5, ig - ioff, jfsr) * xstinv
      IF (ScatOd .EQ. 3) srcmNg(6:9, ig - ioff, jfsr) = 7._8 * srcmNg(6:9, ig - ioff, jfsr) * xstinv
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

NULLIFY (Pin)
NULLIFY (CellInfo)
! ----------------------------------------------------

END SUBROUTINE SetRtP1srcNM