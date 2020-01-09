MODULE DeplLib_MOD
USE PARAM
USE DeplType
USE MatExp_Mod, ONLY : Mat_Type

!TYPE(FisYield_TYPE), POINTER :: FisYeild(:)
!TYPE(ATOMKIND), POINTER :: AtomLib0(:), AtomLib1(:)


!INTEGER, POINTER :: MapMatId2ANum(:)
!INTEGER, POINTER :: MapMatId2IsoWt(:)
!INTEGER, POINTER :: MapMatId2State(:)

TYPE(DeplLib_Type) :: DeplLib(nThreadMax)
!TYPE(DeplCntl_Type) :: DeplCntl
!TYPE(DeplVars_Type) :: DeplVars

INTEGER :: NISODEP, NATOMDEP
!INTEGER, POINTER :: IFISYD(:, :)
!INTEGER, POINTER :: IM2Z(:),IM2A(:),IM2ST(:)

CONTAINS
FUNCTION GetnDeplIso()
INTEGER :: GetnDeplIso
GetnDeplIso = DeplLib(1)%NISODEP
END FUNCTION

FUNCTION GetDeplIsoList(NISO)
INTEGER :: GetDeplIsoList(NISO)
TYPE(ATOMKIND), POINTER :: Lib(:)
INTEGER :: I, J, K
INTEGER :: IMAT

Lib => DeplLib(1)%AtomLib0
DO I = 1, DeplLib(1)%nAtomDep
  DO J = Lib(I)%IB, Lib(I)%IE
    DO K = 0, Lib(I)%A(J)%NSTATE
      IMAT = Lib(I)%A(J)%STAT(K)%IMAT
      IF (IMAT.EQ.0) CYCLE
      GetDeplIsoList(IMAT) = 1000 * I + J + 100 * K      
    ENDDO
  ENDDO
ENDDO
NULLIFY(Lib)
END FUNCTION
END MODULE