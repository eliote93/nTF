MODULE Material_Mod
USE PARAM
USE TYPEDEF, ONLY : Mixture_Type

TYPE(Mixture_Type), POINTER :: Mixture(:)
INTEGER :: nMixType = 0
ENd MODULE