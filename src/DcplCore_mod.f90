MODULE DcplCore_Mod
USE PARAM
USE TYPEDEF,  ONLY : CoreInfo_Type     ,DcplInfo_Type     ,DcplFxrInfo_Type        &
                    ,GroupInfo_Type    ,FMInfo_Type       ,CMInfo_Type             &
                    ,CMFDLS_TYPE       ,DcplCmfdLs_Type   ,ThInfo_Type             &
                    ,MicXsFtn_Type
IMPLICIT NONE
SAVE
TYPE(DcplInfo_Type) :: DcplInfo
TYPE(CoreInfo_Type), POINTER :: DcplCore(:)
TYPE(DcplFxrInfo_Type), POINTER :: DcplFxr(:)
TYPE(FmInfo_Type), POINTER :: DcplFmInfo(:, :)
TYPE(CmInfo_Type), POINTER :: DcplCmInfo(:, :)
TYPE(DcplCmfdLs_Type), POINTER :: DcplCmfdLs(:, :)
TYPE(ThInfo_Type), POINTER :: DcplThInfo(:)
TYPE(MicXsFtn_Type), POINTER :: MicXsFtn(:, :)
END MODULE