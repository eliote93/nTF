MODULE Cooling_mod
!Module for the fuel cooling treatment at the SNF pool  
USE PARAM
USE TYPEDEF,        ONLY : FXRInfo_Type,     GroupInfo_TYPE,           PE_TYPE,    &
                           CoreInfo_Type,    FmInfo_Type
USE DeplType,       ONLY : DeplVars_Type,      DeplLib_Type,        DeplCntl_Type,          &
                           MatExp_Type
USE DeplLib_Mod,    ONLY : DeplLib 
USE Depl_Mod,       ONLY : DeplVars,          DeplCntl,                            &
                           Init_Depl,         AllocDeplFxrMem
USE BasicOperation, ONLY : CP_CA
USE ALLOCS
IMPLICIT NONE

REAL, PRIVATE, POINTER :: IsoNum(:, :)
INTEGER, PRIVATE, POINTER :: MapXs2Dep(:), MapDep2Xs(:)
INTEGER, PRIVATE :: nIsoDep, nIsoLib


CONTAINS
SUBROUTINE Init_Cooling(Core, FmInfo, GroupInfo, PE)
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(PE_Type) :: PE

IF(.NOT. DeplCntl%lInitDepl) THEN
  DeplCntl%lInitDepl = .TRUE.
  CALL Init_Depl(DeplLib, DeplVars, GroupInfo, PE)
  CALL AllocDeplFxrMem(Core, FmInfo%Fxr, GroupInfo, PE)
ENDIF

nIsoLib = GroupInfo%ntiso; nIsoDep = DeplLib(1)%nIsoDep
ALLOCATE(IsoNum(NIsoDep, Pe%nDeplThread))
MapXs2Dep => DeplVars(1)%MapXs2Dep
MapDep2Xs => DeplVars(1)%MapDep2Xs
END SUBROUTINE

SUBROUTINE FxrCooling(Fxr, CoolingTime, tid, lXeDyn)
USE Depl_mod,      ONLY : MatExp,             AllocMatExpVec,      InitDeplXS,             &
                         CopyIsoNumVector,    MakeDeplMat
USE nuclidmap_mod, ONLY : iposiso,            PndCrit,             nuclide
USE MatExp_mod,    ONLY : MatExpSolver
USE XsUtil_mod,    ONLY : SetXeDynEnv
IMPLICIT NONE
TYPE(FxrInfo_Type) :: Fxr
REAL :: CoolingTime
INTEGER :: tid
LOGICAL :: lXeDyn

REAL, POINTER :: pnum(:)
REAL :: phi1g
INTEGER, POINTER :: idiso(:)
INTEGER :: niso
INTEGER :: i, j, id_xs, id_lib, id_depl


niso = Fxr%niso_depl
pnum => Fxr%Pnum; idiso => Fxr%idiso
MatExp(Tid)%Mat => DeplVars(tid)%Dmat; MatExp(Tid)%nIsoDepl = nIsoDep
!COPY Isotope Information
CALL CP_CA(IsoNum(1:nIsoDep, tid), 0._8, nIsoDep)
CALL SetXeDynEnv(IdIso, pnum, Fxr%niso, Fxr%niso_depl, nIsoLib)
DO i = 1, niso
  j = IdIso(i); id_lib = iPosIso(j)
  id_depl = MapXs2Dep(id_lib)
  IsoNum(id_depl, tid) = pnum(i)
!  WRITE(98, *) IdIso(i), pnum(i)
ENDDO
!
!RETURN
!Allocate Depletion Matrix
IF(.NOT. MatExp(Tid)%lAllocVec) THEN
  CALL AllocMatExpVec(MatExp(tid), nIsoDep)
ENDIF
CALL InitDeplXS(DeplLib(tid))
CALL CopyIsoNumVector(MatExp(tid)%Viso0, IsoNum(1:nIsoDep, tid), nIsoDep, epsm30)
phi1g = 0
!DO i = 1, nIsoDep
!  WRITE(97, *)  MatExp(tid)%Viso0
!ENDDO
!Depletion Calculation
CALL MakeDeplMat(MatExp(tid)%Mat, DeplLib(tid), Phi1g, CoolingTime)
CALL MatExpSolver(MatExp(tid))
CALL CopyIsoNumVector(IsoNum(1:nIsoDep, tid), MatExp(tid)%VIso, nIsoDep, epsm30)
!
!Assign the result
CALL CP_CA(IdIso, 0, nIsoLib)
CALL CP_CA(pnum, 0._8, nIsoLib)
niso = 0
!DO i = 1, nIsoDep
!  WRITE(98, *)  MatExp(tid)%Viso0(i), MatExp(tid)%VIso(i), MatExp(tid)%VIso(i)/MatExp(tid)%Viso0(i)
!ENDDO
!

DO i = 1, nIsoDep
  j = MapDep2Xs(i)
  IF(j .EQ. 0) CYCLE
  IF(abs(IsoNum(i, tid)) .LT. epsm20) CYCLE
  IF(IsoNum(i, tid) .LT. pndcrit(j)) CYCLE
  niso = niso + 1
  Id_Xs = nuclide(j)
  pnum(niso) = IsoNum(i, tid)
  idiso(niso) = Id_Xs
ENDDO
Fxr%niso = niso

DO i = 1, nIsoDep
  j = MapDep2Xs(i)
  IF(j .EQ. 0) CYCLE
  IF(abs(IsoNum(i, tid)) .LT. epsm20) CYCLE
  IF(IsoNum(i, tid) .GT. pndcrit(j)) CYCLE
  niso = niso + 1
  Id_Xs = nuclide(j)
  pnum(niso) = IsoNum(i, tid)
  idiso(niso) = Id_Xs
ENDDO
Fxr%niso_depl = niso

IF(lXeDyn) CALL SetXeDynEnv(IdIso, pnum, Fxr%niso, Fxr%niso_depl, nIsoLib)

NULLIFY(MatExp(tid)%Mat)
NULLIFY(pnum, idiso)
END SUBROUTINE

SUBROUTINE Fin_Cooling
  IMPLICIT NONE
  DEALLOCATE(IsoNum)
END SUBROUTINE Fin_Cooling

END MODULE

