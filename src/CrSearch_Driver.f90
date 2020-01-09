SUBROUTINE CrSearch()
USE PARAM
USE GEOM,             ONLY : Core, ng
USE Core_mod,         ONLY : CMInfo, FmInfo, THInfo, GroupInfo, GcGroupInfo,   &
                             eigv, peigv, xkconv, xkcrit  
USE PE_MOD,           ONLY : PE
USE CNTL,             ONLY : nTracerCntl
USE CNTLROD_mod,      ONLY : UpdtCrSearchPosition,   SetCrBankPosition
IMPLICIT NONE
INTEGER :: CrIter
CALL SSEIG()
DO CrIter = 1, 5
  CALL UpdtCrSearchPosition(nTracerCntl%Target_eigv, eigv, .FALSE., PE%MASTER)
  CALL SetCrBankPosition(Core, FmInfo, CmInfo, GroupInfo, nTracerCntl, PE)
  CALL SSEIG()
ENDDO
END SUBROUTINE
