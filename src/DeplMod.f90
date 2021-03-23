!*********************************************************************************************.
!                                                                                             !
!                   *  ------------ Made by LHG, Sep. 2018 ---------------  *                 !
!                   *  Purpose : Atttachable Depletion Module to Any Codes  *                 !
!                   *  ---------------------------------------------------- *                 !
!                   *  Notification : You should initialize the structures  *                 !
!                   *   , DeplFxr, GdFxr and DeplFxrBundles.                *                 !
!                   *   Also, this module requires MatExpMod. If possible,  *                 !
!                   *   you can replace them with your own matrix exponen-  *                 !
!                   *   tial treating module. But you need to revise Depl-  *                 !
!                   *   SolvePostProc.f90 file.                             *                 !
!                   *  ---------------------------------------------------- *                 !
!                   *  Methods : Semi-PC, Utilize SSEIG xs, Whole Iso.(No-  *                 !
!                   *   Split Short/Long Lived), GdQD, XeEq/Tr              *                 !
!                   *  ---------------------------------------------------- *                 !
!                                                                                             !
!                                                                                             !
!*********************************************************************************************'
MODULE HPDeplMod
  USE HPDeplType
  USE CSRMATRIX
  IMPLICIT NONE
  !TYPE(DeplFxrBundle_Type) :: DeplFxrBundle
  !TYPE(DeplLib_Type) :: DeplLib
#if defined(__INTEL_MKL) || defined(__PGI)
  INTERFACE
  ! -------------------------------------- DeplInit.f90 --------------------------------------.
    SUBROUTINE DeplLibInit(libIO, Filename, DeplLib)                                          !
    USE HPDeplType                                                                              !
    INTEGER :: libIO                                                                          !
    CHARACTER(*) :: Filename                                                                  !
    TYPE(DeplLib_Type) :: DeplLib                                                             !
    END SUBROUTINE DeplLibInit                                                                    !

    INTEGER FUNCTION IsoidSrch(IsoArr, NofIso, Idiso)                                                                   !
    INTEGER, INTENT(IN) :: NofIso, Idiso                                                      !
    INTEGER, INTENT(IN) :: IsoArr(*)                                                          !
    END FUNCTION                                                                              !
  ! ------------------------------------- DeplSetSys.f90 -------------------------------------!
    SUBROUTINE SetHmkg0(DeplFxrBundle, DeplLib)                                               !
    USE HPDeplType                                                                              !
    TYPE(DeplFxrBundle_Type) :: DeplFxrBundle                                                 !
    TYPE(DeplLib_Type) :: DeplLib                                                             !
    END SUBROUTINE                                                                            !
                                                                                              !
    SUBROUTINE GetDelT(DeplFxrBundle, BurnUP, lDelBU)                                         !
    USE HPDeplType                                                                              !
    TYPE(DeplFxrBundle_Type) :: DeplFxrBundle                                                 !
    REAL(8) :: BurnUP  ! lDelBU = .TRUE. -> BurnUP(MWD/kg) / lDelBU = .FALSE. -> BurnUP(days) !
    LOGICAL :: lDelBU                                                                         !
    END SUBROUTINE                                                                            !
                                                                                              !
    INTEGER FUNCTION CalSizeSysBun(DeplLib, ByteSys, Scale, NsysMax, IsCRAM)
    USE HPDeplType
    IMPLICIT NONE
    TYPE(DeplLib_Type) :: DeplLib
    INTEGER :: ByteSys, Scale, NsysMax
    LOGICAL :: IsCRAM
    END FUNCTION

    SUBROUTINE SetDeplSys(DeplLib, DeplFxrBundle, Nsys, ifxrbeg, DeplSysBundle, lGd, Nths, nSubStp)   !
    USE HPDeplType                                                                              !
    TYPE(DeplLib_Type) :: DeplLib                                                             !
    TYPE(DeplFxrBundle_Type) :: DeplFxrBundle                                                 !
    INTEGER :: Nsys, ifxrbeg                                                                  !
    TYPE(DeplSysBundle_Type) :: DeplSysBundle                                                 !
    LOGICAL :: lGd
    INTEGER :: Nths
    INTEGER, OPTIONAL :: nSubStp                                                              !
    END SUBROUTINE                                                                            !
                                                                                              !
    SUBROUTINE DestroySys(DeplSysBundle)                                                      !
    USE HPDeplType                                                                            !
    IMPLICIT NONE                                                                             !
    TYPE(DeplSysBundle_Type) :: DeplSysBundle                                                 !
    END SUBROUTINE                                                                            !
                                                                                              !
#ifdef __PGI                                                                                  !
    SUBROUTINE DestroySys_wPoint(DeplSysBundle, Nvec, Solvec, lcsrT)                                 !
    USE HPDeplType                                                                            !
    TYPE(DeplSysBundle_Type) :: DeplSysBundle                                                 !
    REAL(8), POINTER :: Nvec(:), Solvec(:)                                                    !
    LOGICAL :: lcsrT
    END SUBROUTINE                                                                            !

    SUBROUTINE DestroySysnVec(DeplSysBundle, lcsrT)
    USE CSRMATRIX
    USE HPDeplType
    IMPLICIT NONE
    TYPE(DeplSysBundle_Type) :: DeplSysBundle
    LOGICAL :: lcsrT
    END SUBROUTINE

    SUBROUTINE SetDeplSys_woCSRT(DeplLib, DeplFxrBundle, Nsys, ifxrbeg, DeplSysBundle, lGd, Nths, nSubStp)
    USE HPDeplType
    USE CSRMATRIX
    IMPLICIT NONE
    TYPE(DeplLib_Type) :: DeplLib
    TYPE(DeplFxrBundle_Type) :: DeplFxrBundle
    INTEGER :: Nsys, ifxrbeg
    TYPE(DeplSysBundle_Type) :: DeplSysBundle
    LOGICAL :: lGd
    INTEGER :: Nths
    INTEGER, OPTIONAL :: nSubStp
    END SUBROUTINE
                                                                                              !
    SUBROUTINE CopySolVec(DeplSysBundle, NofIso, Solvec, lFM, stream)                              !
    USE HPDeplType                                                                            !
    USE CUDAFOR                                                                               !
    TYPE(DeplSysBundle_Type) :: DeplSysBundle                                                 !
    INTEGER :: NofIso                                                                         !
    REAL(8) :: Solvec(:)                                                             !
    LOGICAL :: lFM
    INTEGER(KIND = cuda_stream_kind) :: stream                                                !
    END SUBROUTINE                                                                            !
                                                                                              !
    SUBROUTINE DestroyVecs_wCopy(DeplFxrBundle, ifxrbeg, Nsys, NofIso, Nvec,&                 !
      Solvec, lCorrector, lGd, lFM, stream)                                                        !
    USE HPDeplType                                                                            !
    USE CUDAFOR                                                                               !
    TYPE(DeplFxrBundle_Type) :: DeplFxrBundle                                                 !
    INTEGER :: ifxrbeg, Nsys, NofIso                                                          !
    REAL(8), POINTER :: Nvec(:), Solvec(:)                                                    !
    LOGICAL :: lCorrector, lGd, lFM                                                                !
    INTEGER(KIND = cuda_stream_kind) :: stream                                                !
    END SUBROUTINE                                                                            !
#endif                                                                                        !
  ! ------------------------------------ DeplSetVars.f90 -------------------------------------!
    SUBROUTINE QuadCoefGen(xs,ys,cs)                                                          !
    REAL(8) :: xs(1:), ys(1:), cs(1:)                                                         !
    END SUBROUTINE QuadCoefGen                                                                !
                                                                                              !
    REAL(8) FUNCTION QuadFunc(cs,x)                                                           !
    REAL(8) :: cs(:), x                                                                       !
    END FUNCTION                                                                              !
                                                                                              !
    SUBROUTINE PostCorrection(DeplFxrBundle, lCorrector, lfirst, Nths, pndcrit)                        !
    USE HPDeplType                                                                            !
    USE OMP_LIB                                                                               !
    TYPE(DeplFxrBundle_Type) :: DeplFxrBundle                                                 !
    LOGICAL :: lCorrector, lfirst                                                             !
    INTEGER :: Nths                                                                           !
    REAL(8), OPTIONAL :: pndcrit(:)
    END SUBROUTINE                                                                            !
                                                                                              !
    SUBROUTINE SetGdVars(DeplFxrBundle, lCorrector, Nths, pndcrit)                                     !
    USE HPDeplType                                                                            !
    USE OMP_LIB                                                                               !
    TYPE(DeplFxrBundle_Type) :: DeplFxrBundle                                                 !
    LOGICAL :: lCorrector                                                                     !
    INTEGER :: Nths                                                                           !
    REAL(8), OPTIONAL :: pndcrit(:)
    END SUBROUTINE                                                                            !
  ! -------------------------------- DeplSolvePostproc.f90 -----------------------------------!
#ifdef __INTEL_MKL
    SUBROUTINE SolveDeplSys(DeplSysBundle, DeplFxrBundle, lCorrector, lGd, SolverTyp, Nths)                         !
    USE HPDeplType                                                                              !
    USE CSRMATRIX                                                                             !
    TYPE(DeplSySBundle_Type) :: DeplSysBundle                                                 !
    TYPE(DeplFxrBundle_Type) :: DeplFxrBundle                                                 !
    LOGICAL :: lCorrector, lGd                                                                !
    INTEGER :: SolverTyp, Nths
    END SUBROUTINE                                                                            !
#endif
                                                                                              !
    SUBROUTINE UpdatePnumDepl(DeplFxrBundle, lCorrector, lSavePre)                                      !
    USE HPDeplType                                                                              !
    TYPE(DeplFxrBundle_Type) :: DeplFxrBundle                                                 !
    LOGICAL :: lSavePre
    LOGICAL :: lCorrector                                                                     !
    END SUBROUTINE                                                                            !

    SUBROUTINE UpdatePnumDeplGd(DeplFxrBundle, lCorrector)
    USE HPDeplType
    IMPLICIT NONE
    TYPE(DeplFxrBundle_Type) :: DeplFxrBundle
    LOGICAL :: lCorrector
    END SUBROUTINE UpdatePnumDeplGd!

    SUBROUTINE UpdatePnumSS(DeplFxrBundle, lCorrector)                                        !
    USE HPDeplType                                                                              !
    TYPE(DeplFxrBundle_Type) :: DeplFxrBundle                                                 !
    LOGICAL :: lCorrector                                                                     !
    END SUBROUTINE                                                                            !
                                                                                              !
    SUBROUTINE UpdateBU(DeplFxrBundle, lCorrector)                                                        !
    USE HPDeplType                                                                              !
    TYPE(DeplFxrBundle_Type) :: DeplFxrBundle                                                 !
    LOGICAL :: lCorrector
    END SUBROUTINE                                                                            !
                                                                                              !
    SUBROUTINE UpdateEqXe(DeplLib, DeplFxrBundle)                                             !
    USE HPDeplType                                                                              !
    TYPE(DeplLib_Type) :: DeplLib                                                             !
    TYPE(DeplFxrBundle_Type) :: DeplFxrBUndle                                                 !
    END SUBROUTINE                                                                            !
#ifdef __PGI
    SUBROUTINE cuSolveDeplSys(DeplSysBundle, DeplFxrBundle, lCorrector, lGd, SolverTyp, Nths)
    USE HPDeplType
    TYPE(DeplSysBundle_Type) :: DeplSysBundle
    TYPE(DeplFxrBundle_Type) :: DeplFxrBundle
    LOGICAL :: lCorrector, lGd
    INTEGER :: Nths, SolverTyp
    END SUBROUTINE
#endif
  ! ------------------------------------------------------------------------------------------'
  END INTERFACE
#endif
END MODULE
