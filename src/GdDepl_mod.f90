MODULE GdDepl_mod
USE PARAM
IMPLICIT NONE
INTEGER, SAVE :: IsoList(7), IsoLoc(7)
DATA IsoList /64152, 64154, 64155, 64156, 64157, 64158, 64160/

INTERFACE

SUBROUTINE GdFxrBurnUp(DeplVars, DeplLib, DeplCntl)
USE PARAM
USE DeplType,    ONLY : DeplVars_Type,      DeplLib_Type,        DeplCntl_Type
IMPLICIT NONE

TYPE(DeplVars_Type) :: DeplVars
TYPE(DeplLib_Type) :: DeplLib
TYPE(DeplCntl_Type) :: DeplCntl

END SUBROUTINE

SUBROUTINE UpdateDeplGdFxrInfo(Fxr, DeplVars, GroupInfo, lCorrectStep, lXeDyn)
USE PARAM
USE TYPEDEF,            ONLY : FxrInfo_Type,       GroupInfo_Type
USE DeplType,           ONLY : DeplVars_Type
IMPLICIT NONE
TYPE(FxrInfo_Type) :: Fxr
TYPE(DeplVars_Type) :: DeplVars
TYPE(GroupInfo_Type) :: GroupInfo
LOGICAL :: lCorrectStep, lXeDyn

END SUBROUTINE

SUBROUTINE GdXsFunction(DeplVars, DeplXs, Fxr, lCorrectStep)
USE PARAM
USE TypeDef,     ONLY : FxrInfo_Type
USE DeplType,    ONLY : DeplXs_Type,      DeplVars_Type
IMPLICIT NONE

TYPE(DeplVars_Type) :: DeplVars
TYPE(DeplXs_Type) :: DeplXS
TYPE(FxrInfo_Type) :: Fxr
LOGICAL :: lCorrectStep
END SUBROUTINE

SUBROUTINE UpdateGdXs(DeplVars, N155)
USE PARAM
USE DeplType,    ONLY : DeplVars_Type
IMPLICIT NONE
TYPE(DeplVars_Type) :: DeplVars
REAL :: N155

END SUBROUTINE

END INTERFACE

CONTAINS
SUBROUTINE Init_GdDepl(DeplVars)
USE DeplType,      ONLY : DeplVars_Type
USE nuclidmap_mod, ONLY : iposiso
IMPLICIT NONE
TYPE(DeplVars_Type) :: DeplVars
INTEGER :: i, id, id_lib
DO i = 1, 7
 id = IsoList(i); id_lib = iposiso(id) 
 IsoLoc(i) = DeplVars%MapXs2Dep(id_lib)
ENDDO
END SUBROUTINE

END MODULE

!SUBROUTINE GdDepl_mod(Fxr, DeplVars, DeplXS, DeplCntl)
!USE PARAM
!USE TypeDef,     ONLY : FxrInfo_Type
!USE DeplType,    ONLY : DeplCntl_Type,    DeplXs_Type,      DeplVars_Type
!USE nuclidmap_mod,ONLY : iposiso
!TYPE(FxrInfo_Type) :: Fxr
!TYPE(DeplVars_Type) :: DeplVars
!TYPE(DeplXs_Type) :: DeplXS
!TYPE(DeplCntl_Type) :: DeplCntl
!
!REAL :: DelT0, DelT
!REAL :: xs1g(64152:64160), pnum0(64152:64160), pnum(64152:64160, 0:1), RR(64152:64160)
!
!INTEGER :: niso
!INTEGER :: i, istep, id_lib
!
!xs1g = 0
!niso = Fxr%niso_depl
!DO i = 1, niso
!  id = Fxr%idiso(i)
!  IF(id .LT. 64152 .OR. id .GT. 64160) CYCLE
!  IF(id .EQ. 64153) CYCLE
!  IF(id .EQ. 64159) CYCLE  
!  xs1g(id) = DeplXs%xsa(i)
!ENDDO
!
!pnum0 = 0
!niso = Fxr%niso_past
!DO i = 1, niso
!  id = Fxr%idiso_past(i)
!  IF(id .LT. 64152 .OR. id .GT. 64160) CYCLE
!  IF(id .EQ. 64153) CYCLE
!  IF(id .EQ. 64159) CYCLE  
!  pnum0(id) = Fxr%pnum_past(i)
!ENDDO
!DelT0 = DeplCntl%tsec
!DelT = DelT0/10000
!RR = 0
!DO id = 64152, 64160
!  IF(id .EQ. 64153) CYCLE
!  IF(id .EQ. 64159) CYCLE
!  RR(id) = xs1g(id) * DeplXs%phi1g * 1.0E-24_8
!ENDDO
!pnum(:, 0) = pnum0(:)
!DO istep = 1, 10000
!  id = 64152
!  pnum(id, 1) = pnum(id, 0) - DelT * pnum(id, 0) * RR(id)
!  id = 64160
!  pnum(id, 1) = pnum(id, 0) - DelT * pnum(id, 0) * RR(id)
!  pnum(64153, 0) = 0
!  RR(64153) = 0
!  DO id = 64154, 64158
!    pnum(id, 1) = pnum(id,0) + DelT * (pnum(id-1, 0) * RR(id-1) - pnum(id, 0) * RR(id))
!  ENDDO 
!  pnum(:, 0) = pnum(:, 1)
!ENDDO
!
!DO id = 64152, 64160
!  IF(id .EQ. 64153) CYCLE
!  IF(id .EQ. 64159) CYCLE
!  id_lib = iPosIso(id)
!  i = DeplVars%MapXs2Dep(id_lib)
!  
!  DeplVars%IsoNum(i) = pnum(id, 1)
!ENDDO
!END SUBROUTINE