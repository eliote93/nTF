#include <defines.h>
SUBROUTINE Init_Depl(DeplLib, DeplVars, GroupInfo, PE)
USE PARAM
USE TYPEDEF,        ONLY : GroupInfo_Type,    PE_TYPE
USE DeplType,       ONLY : DeplLib_Type,      DeplVars_Type,               &
                           DeplCntl_Type,                                  &
                           FisYield_Type,     ATOMKIND
!USE Depl_Mod,       ONLY : SigFold
USE GdDepl_mod,     ONLY : Init_GdDepl
USE FILES,          ONLY : FILENAME,          DeplFileIdx,                 &
                           io_quick
USE MatExp_Mod,     ONLY : Mat_Type
USE CNTL,           ONLY : nTracerCntl_Type
USE Allocs
IMPLICIT NONE

TYPE(DeplLib_Type) :: DeplLib(nThreadMax)
TYPE(DeplVars_Type) :: DeplVars(nThreadMax)
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(PE_TYPE) :: PE

TYPE(FisYield_Type), POINTER :: FisYield(:)
TYPE(ATOMKIND), POINTER :: AtomLib0(:), AtomLib1(:)

INTEGER :: i,j

!CALL ReadDeplFile(io_quick, FileName(DeplFileIdx), DeplLib(1))

!CALL Dmalloc(SigNfOld, GroupInfo%ntiso, GroupInfo%ng, nThreadMax)
!CALL Dmalloc(SigFOld, GroupInfo%ntiso, GroupInfo%ng, nThreadMax)
FisYield => DeplLib(1)%FisYield
AtomLib0 => DeplLib(1)%AtomLib0

CALL SetDepMappingVec(DeplLib(1))
CALL SetChildrenNFisProduct(DeplLib(1))
!CALL setXeDynVar(DeplLib(1))
!CALL setKappa(DeplLib(1))

DO i = 2, PE%nDeplThread, 1
  CALL ReadDeplFile(io_quick, FileName(DeplFileIdx), DeplLib(i),j)
  CALL SetDepMappingVec(DeplLib(i))
  CALL SetChildrenNFisProduct(DeplLib(i))
!  DeplLib(i)%AtomLib0 => DeplLib(1)%AtomLib0
!  DeplLib(i)%FisYield => DeplLib(1)%FisYield
!  DeplLib(i)%NISODEP = DeplLib(1)%NISODEP
!  DeplLib(i)%NFISYLD = DeplLib(1)%NFISYLD
!  DeplLib(i)%NATOMDEP = DeplLib(1)%NATOMDEP
!  DeplLib(i)%MapMatId2ANum => DeplLib(1)%MapMatId2ANum
!  DeplLib(i)%MapMatId2IsoWt => DeplLib(1)%MapMatId2IsoWt
!  DeplLib(i)%MapMatId2State => DeplLib(1)%MapMatId2State
ENDDO
DO i = 1, PE%nDeplThread
  ALLOCATE(DeplVars(i)%DMAT); DeplVars(i)%nIsoDepl = DeplLib(i)%NISODEP

  CALL DMALLOC(DeplVars(i)%IsoNum, DeplVars(i)%nIsoDepl)
  CALL DMALLOC(DeplVars(i)%BurnUpXs, 4, DeplVars(i)%nIsoDepl)

  CALL SetDeplMatDim(DeplVars(i)%DMAT, DeplLib(i))
  CALL MakeAtomLib1(DeplLib(i))
  CALL MakeXsDeplLibMap(DeplVars(i), GroupInfo)
  CALL Init_GdDepl(DeplVars(i))
  DeplVars(i)%tid = i
ENDDO

END SUBROUTINE

SUBROUTINE InitCoreFollow(DeplCntl, nTracerCntl)
USE PARAM
USE TYPEDEF,        ONLY :
USE DEPLTYPE,       ONLY : DeplCntl_Type,         CoreState_Type
USE CNTL,           ONLY : nTracerCntl_Type
IMPLICIT NONE
TYPE(DeplCntl_Type) :: DeplCntl
TYPE(nTracerCntl_Type) :: nTracerCntl

LOGICAL :: lchg
INTEGER :: nburnup
INTEGER :: i

DeplCntl%CoreState%RelPow(0) = nTracerCntl%PowerLevel
DeplCntl%CoreState%FlowRate(0) = 1.00
DeplCntl%CoreState%BoronPPM(0) = nTracerCntl%BoronPPM
DeplCntl%CoreState%T_efpd(0) = 0
DeplCntl%CoreState%T_mwdkghm(0) = 0
nTracerCntl%lBoronSearch = .FALSE.

nburnup = DeplCntl%nBurnupStep
!Copy Depeltion Time
DeplCntl%CoreState%T_efpd(1:nburnup) = DeplCntl%T_efpd(1:nburnup)
DeplCntl%CoreState%T_mwdkghm(1:nburnup) = DeplCntl%T_mwdkghm(1:nburnup)
DeplCntl%CoreState%lStateChg(0:nburnup) = .FALSE.

DO i = 1, nburnup
  lChg = .FALSE.
  IF(abs(DeplCntl%CoreState%RelPow(i-1) - DeplCntl%CoreState%RelPow(i)) .GT. epsm3) lChg = .TRUE.
  IF(abs(DeplCntl%CoreState%FlowRate(i-1) - DeplCntl%CoreState%FlowRate(i)) .GT. epsm3) lChg = .TRUE.
  DeplCntl%CoreState%lStateChg(i) = lChg
ENDDO

END SUBROUTINE

SUBROUTINE SetDepMappingVec(DeplLib)
USE PARAM
USE DeplType,       ONLY : DeplLib_Type,  ATOMKIND
USE ALLOCS
IMPLICIT NONE
TYPE(DeplLib_Type) :: DeplLib
TYPE(ATOMKIND), POINTER :: AtomLib0(:)
INTEGER :: NISODEP, NATOMDEP
INTEGER :: I, J, K, IMAT
AtomLib0 => DeplLib%AtomLib0
NISODEP = DeplLib%NISODEP
NATOMDEP = DeplLib%NATOMDEP
CALL Dmalloc(DeplLib%MapMatId2ANum, NISODEP)
CALL Dmalloc(DeplLib%MapMatId2IsoWt, NISODEP)
CALL Dmalloc(DeplLib%MapMatId2State, NISODEP)

DO I = 1, NATOMDEP
  DO J = AtomLib0(i)%IB, AtomLib0(i)%IE
    DO K = 0, AtomLib0(i)%A(J)%NSTATE
      IF(AtomLib0(i)%A(J)%STAT(K)%IMAT .EQ. 0) CYCLE
      IMAT = AtomLib0(i)%A(J)%STAT(K)%IMAT
      DeplLib%MapMatId2ANum(IMAT) = I
      DeplLib%MapMatId2IsoWt(IMAT) = J
      DeplLib%MapMatId2State(IMAT) = K
    ENDDO
  ENDDO
ENDDO
NULLIFY(AtomLib0)
END SUBROUTINE

SUBROUTINE SetChildrenNFisProduct(DeplLib)
USE PARAM
USE DeplType,       ONLY : DeplLib_Type,  ATOMKIND
USE ALLOCS
IMPLICIT NONE
TYPE(DeplLib_Type) :: DeplLib
TYPE(ATOMKIND), POINTER :: AtomLib0(:)

REAL, PARAMETER :: ERR = 1.E-25_8

INTEGER :: NISODEP, NATOMDEP, N
INTEGER :: I, J, K
REAL :: RAMBDA
NISODEP = DeplLib%NISODEP; NATOMDEP = DeplLib%NATOMDEP

AtomLib0 => DeplLib%AtomLib0

DO I = 1, NATOMDEP
  DO J = AtomLib0(i)%IB, AtomLib0(i)%IE
    DO K = 0, AtomLib0(i)%A(J)%NSTATE
      IF(AtomLib0(i)%A(J)%STAT(K)%IMAT .EQ. 0) CYCLE
      N = 0
      IF(AtomLib0(I)%A(J)%STAT(K)%RAMBDA .GT. ERR) THEN
        CALL SetChildrenFromDecay(DeplLib, I, J, K)
      ENDIF
      IF(AtomLib0(I)%A(J)%STAT(K)%XS(0).GT.ERR) THEN
        CALL SetChildrenFromReaction(DeplLib, I, J, K)
      ENDIF
      IF(AtomLib0(I)%A(J)%STAT(K)%IGRP .EQ. 3) THEN
        CALL SetFissionReaction(DeplLib, I, J, K)
      ENDIF
    ENDDO
  ENDDO
ENDDO
CONTINUE
END SUBROUTINE

SUBROUTINE SetChildrenFromDecay(DeplLib, IZ, IA, IST)
USE PARAM
USE DeplType,       ONLY : DeplLib_Type,  ATOMKIND,   STATEKIND
USE ALLOCS
IMPLICIT NONE
TYPE(DeplLib_Type) :: DeplLib
TYPE(ATOMKIND), POINTER :: AtomLib0(:)
TYPE(STATEKIND), POINTER :: mySTAT

REAL, PARAMETER :: ERR = 1.E-25_8

INTEGER ::  IDECAYTO(3,8)
!DATA IDECAYTO/1,0,1, -1,0,-1, -1,0,1, -2,-4,-1, 0,0,0, 1,-1,-1, 0,0,0, 1,0,0/
DATA IDECAYTO/1,0,1, -1,0,0, -1,0,1, -2,-4,0, 0,0,0, 0,0,0, 1,-1,0, 1,0,0/

INTEGER :: NISODEP, NATOMDEP, N, M
INTEGER :: IZ, IA, IST
INTEGER :: IZ2, IA2, IST2
INTEGER :: IA2IND
INTEGER :: I, J, K
REAL :: RAMBDA

NISODEP = DeplLib%NISODEP; NATOMDEP = DeplLib%NATOMDEP
AtomLib0 => DeplLib%AtomLib0
mySTAT => AtomLib0(IZ)%A(IA)%STAT(IST)
N=0
RAMBDA = mySTAT%RAMBDA
DO M = 1, 8
!  IF(M .EQ. 7 .OR. mySTAT%FRAC(M) .LT. ERR) CYCLE
  IF(M .EQ. 6 .OR. mySTAT%FRAC(M) .LT. ERR) CYCLE  
  mySTAT%FRAC(M)=RAMBDA*mySTAT%FRAC(M)
  N=N+1
ENDDO
mySTAT%NTO1 = N
CALL DMALLOC(mySTAT%ITO1, 3, N)
N = 0
DO M = 1, 8
  !IF(M .EQ. 7 .OR. mySTAT%FRAC(M) .LT. ERR) CYCLE
  IF(M .EQ. 6 .OR. mySTAT%FRAC(M) .LT. ERR) CYCLE
  IZ2 = IZ + IDECAYTO(1, M)           !Daugther Atom Number(Z)
  IA2 = IA + IDECAYTO(2, M)           !Daugther Atomic Weight(A)
  !IST2 = IST
  !IF(IDECAYTO(3,M).NE.-1) IST2=IDECAYTO(3,M)
  IST2 = IDECAYTO(3, M)
  IF(IZ2 .GT. NATOMDEP) CYCLE
  IA2IND=(IA2-AtomLib0(IZ2)%IB)*(IA2-AtomLib0(IZ2)%IE)
  IF(AtomLib0(IZ2)%IE .EQ. 0 .OR. IA2IND .GT. 0) CYCLE
  IF(IST2 .GT. AtomLib0(IZ2)%A(IA2)%NSTATE) CYCLE
  IF(AtomLib0(IZ2)%A(IA2)%STAT(IST2)%IMAT.EQ.0) CYCLE
  N = N + 1
  mySTAT%ITO1(1, N) = M
  mySTAT%ITO1(2, N) = AtomLib0(IZ2)%A(IA2)%STAT(IST2)%IMAT
ENDDO
mySTAT%NTO1=N
END SUBROUTINE

SUBROUTINE SetChildrenFromReaction(DeplLib, IZ, IA, IST)
USE PARAM
USE DeplType,       ONLY : DeplLib_Type,  ATOMKIND,   STATEKIND
USE ALLOCS
IMPLICIT NONE
TYPE(DeplLib_Type) :: DeplLib
TYPE(ATOMKIND), POINTER :: AtomLib0(:)
TYPE(STATEKIND), POINTER :: mySTAT

REAL, PARAMETER :: ERR = 1.E-25_8

INTEGER :: IREACTTO(3,7)
DATA IREACTTO/0,1,0, 0,-1,0, -2,-3,0, -1,0,0, 0,1,1, 0,-1,1, 0,-2,0/

INTEGER :: NISODEP, NATOMDEP, N, M, M2
INTEGER :: IZ, IA, IST
INTEGER :: IZ2, IA2, IST2
INTEGER :: IA2IND
INTEGER :: I, J, K
!REAL :: RAMBDA

! REACTION TYPE
!  1: (N,R)  (Z,A,0)->(Z,A+1,0)
!  2: (N,2N) (Z,A,0)->(Z,A-1,0)
!  3: (N,A)  (Z,A,0)->(Z-2,A-3,0) FOR ACTIVATION AND F.P.
!  4: (N,P)  (Z,A,0)->(Z-1,A,0) FOR ACTIVATION AND F.P.
!  5: (N,R') (Z,A,0)->(Z,A+1,1)
!  6: (N,2N') (Z,A,0)->(Z,A-1,1)
!  3: (N,3N) (Z,A,0)->(Z,A-2,0) FOR ACTINIDE
!  4: FISSION FOR ACTINIDE

NISODEP = DeplLib%NISODEP; NATOMDEP = DeplLib%NATOMDEP
AtomLib0 => DeplLib%AtomLib0
mySTAT => AtomLib0(IZ)%A(IA)%STAT(IST)

N=0
DO M = 1, 6
  IF(mySTAT%XS(M) .LT. ERR) CYCLE
  IF(mySTAT%IGRP .EQ. 2 .AND. M .EQ. 4) CYCLE !Actinide => No (N, P) Reaction
  N = N + 1
  IF(M .EQ. 4) N = N + 1 ! PROTON TO HYDROGEN
ENDDO
mySTAT%NTO2 = N
CALL Dmalloc(mySTAT%ITO2, 3, N)
N = 0
DO M = 1, 6
  IF(mySTAT%XS(M) .LT. ERR) CYCLE
  IF(mySTAT%IGRP .EQ. 2 .AND. M .EQ. 4) CYCLE
  M2 = M
  IF(mySTAT%IGRP .EQ. 2 .AND. M .EQ. 3) M2 = 7
  IF(M .EQ. 4 .AND. AtomLib0(IZ)%IE .GT. 1) THEN
    N = N + 1
    mySTAT%ITO2(1, N) = M
    mySTAT%ITO2(2, N) = AtomLib0(1)%A(1)%STAT(0)%IMAT
  ENDIF
  IZ2 = IZ + IREACTTO(1, M2)
  IA2 = IA + IREACTTO(2, M2)
  IST2 = IREACTTO(3, M2)
  IF(IZ2 .GT. NATOMDEP) CYCLE
  IA2IND=(IA2-AtomLib0(IZ2)%IB)*(IA2-AtomLib0(IZ2)%IE)
  IF(AtomLib0(IZ2)%IE .EQ. 0 .OR. IA2IND .GT. 0) CYCLE
  IF(IST2 .GT. ATOMLIB0(IZ2)%A(IA2)%NSTATE)  CYCLE
  IF(AtomLib0(IZ2)%A(IA2)%STAT(IST2)%IMAT .EQ. 0) CYCLE
  N = N + 1
  mySTAT%ITO2(1,N) = M
  mySTAT%ITO2(2,N) = AtomLib0(IZ2)%A(IA2)%STAT(IST2)%IMAT
ENDDO
mySTAT%NTO2 = N

END SUBROUTINE

SUBROUTINE SetFissionReaction(DeplLib, IZ, IA, IST)
USE PARAM
USE DeplType,       ONLY : DeplLib_Type,  ATOMKIND,   STATEKIND
USE ALLOCS
IMPLICIT NONE
TYPE(DeplLib_Type) :: DeplLib
TYPE(ATOMKIND), POINTER :: AtomLib0(:)
TYPE(STATEKIND), POINTER :: mySTAT

REAL, PARAMETER :: ERR = 1.E-25_8

INTEGER :: NISODEP, NATOMDEP, NFISYLD
INTEGER :: N, M
INTEGER :: IZ, IA, IST
INTEGER :: IZ2, IA2, IST2
INTEGER :: IA2IND
INTEGER :: I, J, K

NISODEP = DeplLib%NISODEP; NATOMDEP = DeplLib%NATOMDEP
NFISYLD = DeplLib%NFISYLD
AtomLib0 => DeplLib%AtomLib0
mySTAT => AtomLib0(IZ)%A(IA)%STAT(IST)

N = 0
DO M = 1, NFISYLD
  IF(mySTAT%Y(M) .LT. ERR) CYCLE
  N = N + 1
ENDDO
mySTAT%NFR3 = N
CALL Dmalloc(mySTAT%IFR3, 2, N)

N = 0
DO M = 1, NFISYLD
  IF(mySTAT%Y(M) .LT. ERR) CYCLE
  N = N +1
  mySTAT%IFR3(1, N) = M
ENDDO

END SUBROUTINE

SUBROUTINE SetDeplMatDim(DMAT, DeplLib)
USE PARAM
USE MatExp_Mod,     ONLY : Mat_Type
USE DeplType,       ONLY : DeplLib_Type,      ATOMKIND,       STATEKIND,    &
                           FisYield_Type
USE DEPL_MOD,       ONLY : ElmtLocFindnUpdt
USE ALLOCS
IMPLICIT NONE
TYPE(Mat_Type) :: DMAT
TYPE(DeplLib_Type) :: DeplLib

TYPE(ATOMKIND), POINTER :: AtomLib0(:)
TYPE(STATEKIND), POINTER :: mySTAT
TYPE(FisYield_Type), POINTER :: FisYield(:)

INTEGER :: NISODEP, NATOMDEP, NFISYLD
INTEGER :: nmaxoffdiag
INTEGER :: I, J, K, M, N
INTEGER :: ITO, IMAT, IYD, IFR

NISODEP = DeplLib%NISODEP; NATOMDEP = DeplLib%NATOMDEP
NFISYLD = DeplLib%NFISYLD
AtomLib0 => DeplLib%AtomLib0
FisYield => DeplLib%FisYield

CALL DMALLOC(DMAT%nlmnt, NISODEP)
CALL DMALLOC(DMAT%Diag, NISODEP)

DO I = 1, NATOMDEP
  DO J = AtomLib0(i)%IB, AtomLib0(i)%IE
    DO K = 0, AtomLib0(i)%A(J)%NSTATE
      IF(AtomLib0(i)%A(J)%STAT(K)%IMAT .EQ. 0) CYCLE
      mySTAT => AtomLib0(i)%A(J)%STAT(K)
      DO M = 1, mySTAT%NTO1
        ITO = mySTAT%ITO1(2, M)
        DMAT%nlmnt(ITO) = DMAT%nlmnt(ITO) + 1
      ENDDO

      DO M = 1, mySTAT%NTO2
        ITO = mySTAT%ITO2(2, M)
        DMAT%nlmnt(ITO) = DMAT%nlmnt(ITO) + 1
      ENDDO
      IF(mySTAT%IGRP .NE. 3) CYCLE
      IMAT = mySTAT%IMAT
      DMAT%nlmnt(IMAT) = DMAT%nlmnt(IMAT) + mySTAT%NFR3
      NULLIFY(mySTAT)
    ENDDO
  ENDDO
ENDDO

nmaxoffdiag = 0
DO I = 1, NISODEP
  nmaxoffdiag = MAX(DMAT%nlmnt(I), nmaxoffdiag)
  DMAT%nlmnt(I) = 0
ENDDO
DMAT%nmaxoffdiag = nmaxoffdiag
CALL DMALLOC(DMAT%OffDiag, nmaxoffdiag, NISODEP)
CALL DMALLOC(DMAT%lmntIdx, nmaxoffdiag, NISODEP)

CALL DMALLOC(DMAT%Diag, NISODEP)
CALL DMALLOC(DMAT%OffDiag, nmaxoffdiag, NISODEP)

DO I = 1, NATOMDEP
  DO J = AtomLib0(i)%IB, AtomLib0(i)%IE
    DO K = 0, AtomLib0(i)%A(J)%NSTATE
      IF(AtomLib0(i)%A(J)%STAT(K)%IMAT .EQ. 0) CYCLE
      mySTAT =>  AtomLib0(i)%A(J)%STAT(K)
      IMAT = myStat%IMAT
      DO M = 1, mySTAT%NTO1   !  FROM DECAY CHAIN
        ITO = mySTAT%ITO1(2,M)
        mySTAT%ITO1(3, M) = ElmtLocFindnUpdt(DMAT, ITO, IMAT)
      ENDDO

      DO M = 1, mySTAT%NTO2   !  FROM REACTION TYPE
        ITO = mySTAT%ITO2(2,M)
        mySTAT%ITO2(3, M) = ElmtLocFindnUpdt(DMAT, ITO, IMAT)
      ENDDO

      IF(mySTAT%IGRP.NE.3) CYCLE

      DO M = 1, mySTAT%NFR3   !  FROM FISSION
        IYD = mySTAT%IFR3(1,M)
        IFR = FisYield(IYD)%MATID
        mySTAT%IFR3(2, M) = ElmtLocFindnUpdt(DMAT, IMAT, IFR)
      ENDDO

      NULLIFY(mySTAT)
    ENDDO
  ENDDO
ENDDO
NULLIFY(AtomLib0)
END SUBROUTINE

SUBROUTINE MakeAtomLib1(DeplLib)
USE PARAM
USE MatExp_Mod,     ONLY : Mat_Type
USE DeplType,       ONLY : DeplLib_Type,      ATOMKIND,       STATEKIND,    &
                           FisYield_Type
!USE DEPL_MOD,       ONLY : ElmtLocFindnUpdt
USE ALLOCS
IMPLICIT NONE

TYPE(DeplLib_Type) :: DeplLib
TYPE(ATOMKIND), POINTER :: AtomLib0(:), AtomLib1(:)
TYPE(STATEKIND), POINTER :: mySTAT1, mySTAT0
TYPE(FisYield_Type), POINTER :: FisYield(:)

INTEGER :: NISODEP, NATOMDEP, NFISYLD
INTEGER :: NSTATE
INTEGER :: I, J, K
NISODEP = DeplLib%NISODEP; NATOMDEP = DeplLib%NATOMDEP
NFISYLD = DeplLib%NFISYLD
AtomLib0 => DeplLib%AtomLib0
FisYield => DeplLib%FisYield

ALLOCATE(DeplLib%AtomLib1(NATOMDEP))
AtomLib1 => DeplLib%AtomLib1

DO I = 1, NATOMDEP
  AtomLib1(I)%IB = AtomLib0(I)%IB
  AtomLib1(I)%IE = AtomLib0(I)%IE
  IF(AtomLib0(I)%IE .EQ. 0) CYCLE
  ALLOCATE(AtomLib1(I)%A(AtomLib1(I)%IB:AtomLib1(I)%IE))
  DO J=AtomLib1(I)%IB, AtomLib1(I)%IE
    AtomLib1(I)%A(J)%NSTATE = AtomLib0(I)%A(J)%NSTATE
    IF(AtomLib1(I)%A(J)%NSTATE .EQ. -1) CYCLE
    NSTATE = AtomLib1(I)%A(J)%NSTATE
    ALLOCATE(AtomLib1(I)%A(J)%STAT(0:NSTATE))
    DO K = 0, NSTATE
      mySTAT0 => AtomLib0(I)%A(J)%STAT(K)
      mySTAT1 => AtomLib1(I)%A(J)%STAT(K)
      mySTAT1%IGRP=mySTAT0%IGRP; mySTAT1%NTO1=mySTAT0%NTO1
      mySTAT1%NTO2=mySTAT0%NTO2; mySTAT1%IMAT=mySTAT0%IMAT
      mySTAT1%RAMBDA=mySTAT0%RAMBDA; mySTAT1%FRAC=mySTAT0%FRAC
      mySTAT1%XS=mySTAT0%XS
      CALL DMALLOC(mySTAT1%ITO1, 3, mySTAT1%NTO1)
      CALL DMALLOC(mySTAT1%ITO2, 3, mySTAT1%NTO2)
      IF(mySTAT1%NTO1.GT.0) mySTAT1%ITO1 = mySTAT0%ITO1(:, 1 : mySTAT1%NTO1)
      IF(mySTAT1%NTO2.GT.0) mySTAT1%ITO2 = mySTAT0%ITO2(:, 1 : mySTAT1%NTO2)
      mySTAT1%NFR3 = mySTAT0%NFR3
      IF(mySTAT1%IGRP .NE. 3) CYCLE
      CALL Dmalloc(mySTAT1%Y, NFISYLD)
      CALL Dmalloc(mySTAT1%IFR3, 2,  mySTAT1%NFR3)
      mySTAT1%IFR3 = mySTAT0%IFR3
      mySTAT1%Y = mySTAT0%Y
    ENDDO
  ENDDO
ENDDO
IF(ASSOCIATED(mySTAT0)) NULLIFY(mySTAT0)
IF(ASSOCIATED(mySTAT1)) NULLIFY(mySTAT1)
NULLIFY(AtomLib0, AtomLib1)
NULLIFY(FisYield)
END SUBROUTINE

SUBROUTINE MakeXsDeplLibMap(DeplVars, GroupInfo) !Mapping Between
USE PARAM
USE TYPEDEF,      ONLY : GroupInfo_TYPE
USE DeplType,     ONLY : DeplVars_Type
USE DeplLib_MOD,  ONLY : GetnDeplIso,      GetDeplIsoList
USE nuclidmap_mod,ONLY : iposiso
USE Allocs
IMPLICIT NONE
TYPE(DeplVars_TYPE) :: DeplVars
TYPE(GroupINfo_Type) :: GroupInfo

INTEGER, POINTER ::  MapXs2Dep(:), MapDep2XS(:), IsoList(:)
INTEGER :: nIsoDepl, nIsoXs
INTEGER :: ixs, i

nIsoXs = GroupInfo%ntiso      !Number of total XS isotope
nIsoDepl = GetnDeplIso()      !Number of Depletion Isotope

CALL Dmalloc(DeplVars%MapXs2Dep, nIsoXs)
CALL Dmalloc(DeplVars%MapDep2XS, nIsoDepl)
ALLOCATE(IsoList(nIsoDepl))
MapXs2Dep => DeplVars%MapXs2Dep
MapDep2XS => DeplVars%MapDep2XS
IsoList = GetDeplIsoList(nIsoDepl)
DO I = 1, nIsoDepl  !  MAPPING FILES
  ixs = iposiso(IsoList(i))
  IF(ixs .NE. 0) THEN
    MapDep2Xs(I) = IXS
    MapXs2Dep(IXS) = I
  ENDIF
ENDDO
DEALLOCATE(IsoList)
NULLIFY(MapXs2Dep, MapDep2Xs)
END SUBROUTINE

SUBROUTINE AllocDeplFXRMem(Core, Fxr, GroupInfo, PE)
USE PARAM
USE TYPEDEF,    ONLY : CoreInfo_Type,     FxrInfo_Type,         PE_TYPE,          &
                       GroupInfo_Type,    Pin_Type,             Cell_Type
USE ALLOCS
USE DeplLib_MOD,  ONLY : GetnDeplIso  ! 16/02/11 Depletion timestep bug fixed
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: FXR(:, :)
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(PE_Type) :: PE

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
TYPE(FxrInfo_Type), POINTER :: myFxr

INTEGER :: ntiso, nFxr, nxy, myzb, myze
INTEGER :: FxrIdxSt, nLocalFxr
INTEGER :: ixy, iz, ifxr, icell, i, j
INTEGER :: nIsoDepl   ! 16/02/11 Depletion timestep bug fixed
myzb = PE%myzb; myze = PE%myze
nFxr = Core%nCoreFxr; nxy = Core%nxy
ntiso = GroupInfo%ntiso
nIsoDepl = GetnDeplIso()      !Number of Depletion Isotope ! 16/02/11 Depletion timestep bug fixed
Pin => Core%Pin; Cell => Core%CellInfo

DO iz = myzb, myze
  DO ixy = 1, nxy
    icell = Pin(ixy)%Cell(iz)
    FxrIdxSt = Pin(ixy)%FxrIdxSt; nLocalFxr = Cell(icell)%nFxr
    DO i = 1, nLocalFxr
      iFxr = FxrIdxSt + i - 1
      IF(.NOT. Fxr(iFxr, iz)%lDepl) CYCLE
      myFxr => Fxr(iFxr, iz)

      CALL Dmalloc(myFxr%idiso_past, ntiso)
      CALL Dmalloc(myFxr%pnum_past, ntiso)
      !N.D. Vectors for saving all isotopes specified in Depl. Lib.
      CALL Dmalloc(myFxr%pnum_all, nisodepl)
      CALL Dmalloc(myFxr%pnum_past_all, nisodepl)
      !--- end  ! 16/02/11 Depletion timestep bug fixed
      myFxr%l_pnum_all = .FALSE.
      !myFxr%burnup = 0; myFxr%burnup_past = 0
#ifdef GdHO
      IF(.NOT. myFxr%lGD) CYCLE
#endif
      ALLOCATE(myFxr%DeplPCP)
      ALLOCATE(myFxr%DeplXS1g(-1:0))
      DO j = -1, 0
        CALL Dmalloc(myFxr%DeplXs1g(j)%idiso, ntiso)
        CALL Dmalloc(myFxr%DeplXs1g(j)%xsa, ntiso)
        CALL Dmalloc(myFxr%DeplXs1g(j)%xsf, ntiso)
        CALL Dmalloc(myFxr%DeplXs1g(j)%xsn2n, ntiso)
        CALL Dmalloc(myFxr%DeplXs1g(j)%xsn3n, ntiso)
      ENDDO
    ENDDO
  ENDDO
ENDDO

IF(Associated(myFXR)) NULLIFY(myFXR)

END SUBROUTINE

FUNCTION ElmtLocFindnUpdt(DMAT, I, IMAT)
USE MatExp_Mod, ONLY : Mat_Type
IMPLICIT NONE
TYPE(Mat_Type) :: DMAT
INTEGER :: I, IMAT
INTEGER :: ElmtLocFindnUpdt
INTEGER :: j, n

n = DMAT%nlmnt(I)
DO J = 1, N
  IF(DMAT%lmntIdx(J, I) .EQ. IMAT) EXIT
ENDDO
IF(J .GT. N) THEN
  DMAT%nlmnt(I) = DMAT%nlmnt(I) + 1
  DMAT%lmntIdx(DMAT%nLmnt(I), I) = IMAT
  ElmtLocFindnUpdt = DMAT%nlmnt(I)
ELSE
  ElmtLocFindnUpdt = J
ENDIF
END FUNCTION

SUBROUTINE setXeDynVar(DeplLib)
    USE Depl_Mod
    USE XsLib_mod, only : nelthel,ldiso    
    USE DeplType,  only : DeplLib_Type, ATOMKIND
    implicit none
    
    integer :: iz,ia,is,iz1,ia1,is1,iz2,ia2,is2,nid,i,ifis
    TYPE(DeplLib_Type) :: DeplLib
    TYPE(ATOMKIND), POINTER :: AtomLib0(:)
    
    AtomLib0 => DeplLib%AtomLib0
    
    iz = 54; ia = 635; is = 0
    decayXe135 = AtomLib0(iz)%A(ia)%Stat(is)%RAMBDA
    allocate(yieldXe135(nelthel))
    do i = 1, nelthel
        yieldXe135(i) = 0._8
        if (ldiso(i)%ifis.eq.0) cycle
        nid = ldiso(i)%nid
        iz1 = nid/1000
        ia1 = nid - 1000*iz1
        is1 = 0
        if (ia1.gt.300) then
            is1 = 1
            ia1 = ia1 - 100
        endif
        do ifis = 1, DeplLib%NFISYLD
            iz2 = DeplLib%FisYield(ifis)%AtomNum
            ia2 = DeplLib%FisYield(ifis)%AtomWeight
            is2 = DeplLib%FisYield(ifis)%Stat
            if ((iz1.eq.iz2).and.(ia1.eq.ia2).and.(is1.eq.is2)) exit
        enddo
        if (ifis.gt.DeplLib%NFISYLD) cycle
        yieldXe135(i) = AtomLib0(iz)%A(ia)%Stat(is)%Y(ifis)
    enddo
  
    iz = 53; ia = 635; is = 0
    decayI135 = AtomLib0(iz)%A(ia)%Stat(is)%RAMBDA
    allocate(yieldI135(nelthel))
    do i = 1, nelthel
        yieldI135(i) = 0._8
        if (ldiso(i)%ifis.eq.0) cycle
        nid = ldiso(i)%nid
        iz1 = nid/1000
        ia1 = nid - 1000*iz1
        is1 = 0
        if (ia1.gt.300) then
            is1 = 1
            ia1 = ia1 - 100
        endif
        do ifis = 1, DeplLib%NFISYLD
            iz2 = DeplLib%FisYield(ifis)%AtomNum
            ia2 = DeplLib%FisYield(ifis)%AtomWeight
            is2 = DeplLib%FisYield(ifis)%Stat
            if ((iz1.eq.iz2).and.(ia1.eq.ia2).and.(is1.eq.is2)) exit
        enddo
        if (ifis.gt.DeplLib%NFISYLD) cycle
        yieldI135(i) = AtomLib0(iz)%A(ia)%Stat(is)%Y(ifis)
    enddo

END SUBROUTINE

SUBROUTINE setKappa(DeplLib)
USE Depl_Mod
    USE XsLib_mod, only : nelthel,ldiso    
    USE DeplType,  only : DeplLib_Type, ATOMKIND
    implicit none
    integer :: i,nid,iz1,ia1,is1
    TYPE(DeplLib_Type) :: DeplLib
    TYPE(ATOMKIND), POINTER :: AtomLib0(:)
    
    AtomLib0 => DeplLib%AtomLib0

    do i = 1, nelthel
        ldiso(i)%kappa = 0._8
        if (ldiso(i)%ifis.eq.0) cycle
        nid = ldiso(i)%nid
        iz1 = nid/1000
        ia1 = nid - 1000*iz1
        is1 = 0
        if (ia1.gt.300) then
            is1 = 1
            ia1 = ia1 - 100
        endif
        if (iz1.gt.DeplLib%NATOMDEP) cycle
        if (ia1.gt.AtomLib0(iz1)%IE) cycle
        if (ia1.lt.AtomLib0(iz1)%IB) cycle
        if (is1.gt.AtomLib0(iz1)%A(ia1)%NSTATE) cycle
        if (AtomLib0(iz1)%A(ia1)%Stat(is1)%imat.eq.0) cycle
        ldiso(i)%kappa = AtomLib0(iz1)%A(ia1)%Stat(is1)%kappa
    enddo

END SUBROUTINE
