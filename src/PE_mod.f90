MODULE PE_Mod
USE PARAM
USE TYPEDEF, ONLY : PE_TYPE
IMPLICIT NONE
TYPE(PE_TYPE) :: PE
TYPE(PE_TYPE) :: DcplPE(100)
INTEGER :: myzb, myze
!LOGICAL :: master = .TRUE.
!LOGICAL :: slave = .FALSE.
END MODULE