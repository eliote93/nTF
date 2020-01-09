MODULE FXRVAR_MOD
    IMPLICIT NONE
    LOGICAL :: lBranchRun = .FALSE.
    LOGICAL :: ltmod   = .FALSE.
    LOGICAL :: ltfuel  = .FALSE.
    LOGICAL :: lboron  = .FALSE.
    LOGICAL :: lrho    = .FALSE.
    LOGICAL :: lfxrvar = .FALSE.
    REAL :: tmod(0:10)=0.0_8
    REAL :: tfuel(0:10)=0.0_8
    REAL :: rho(0:10)=0.0_8
    REAL :: bboron(0:10)=0.0_8 ! Branch Boron
    INTEGER :: ntmod, ntfuel, nboron, nrho
    INTEGER :: CurrentVarId, CurrentCaseId
    INTEGER, PARAMETER :: branch_base  = 0
    INTEGER, PARAMETER :: branch_tmod  = 1
    INTEGER, PARAMETER :: branch_tfuel = 2
    INTEGER, PARAMETER :: branch_boron = 3
    INTEGER, PARAMETER :: branch_rho   = 4
ENDMODULE
    
    