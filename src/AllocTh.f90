SUBROUTINE AllocTH()
USE PARAM
USE TYPEDEF,      ONLY : THInfo_Type, PE_TYPE, CoreInfo_Type
USE Geom,         ONLY : Core
USE CORE_mod,     ONLY : THInfo
USE PE_Mod,           ONLY : PE
USE Allocs
IMPLICIT NONE

INTEGER :: myzb, myze, nxy, nfxr, nz

myzb = PE%myzb; myze = PE%myze; nz = Core%nz
nxy = Core%nxy; nfxr = Core%nCoreFxr
CALL dmalloc0(THInfo%RefFuelTemp, 0, nz)
CALL dmalloc0(THinfo%FxrTemp, 1, nfxr, myzb, myze)
END SUBROUTINE