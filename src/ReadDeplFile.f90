SUBROUTINE ReadDeplFile(InDev, filename, DeplLib, NISODEP)
USE DeplType,    ONLY : STATEKIND,     ISOTOPEKIND,       ATOMKIND,        &
                        FisYield_TYPE, DeplLib_Type
USE MatExp_Mod,  ONLY : Mat_Type
!USE DEPL_MOD,    ONLY : NISODEP,       NATOMDEP,                          &
!                        AtomLib0,      AtomLib1,          FisYeild,       &
!                        MapMatId2ANum, MapMatId2IsoWt,    MapMatId2State, &
!                        VISO
!USE FILES,       ONLY : FILENAME, DeplFileIdx
USE IOUTIL,      ONLY : OpenFile
USE ALLOCS
IMPLICIT NONE
TYPE(DeplLib_Type) :: DeplLib
INTEGER :: InDev
character(256) :: filename

character(256) :: oneline

REAL :: TCONV(9)
DATA TCONV/1._8,            60._8,      3600._8,       86400._8,   3.15576E+07_8,   &
           1._8,    3.15576E+10_8, 3.15576E+13_8, 3.15576E+16_8/

TYPE(FisYield_TYPE), POINTER :: FisYield(:)
TYPE(ATOMKIND), POINTER :: AtomLib0(:), AtomLib1(:)

INTEGER :: NISO1, NISO2, NFISYLD, NISODEP
INTEGER :: NATOMDEP
INTEGER :: IZ, IA , IST, IMAT      !Z number, A number, # of States
INTEGER, POINTER :: NIDINP(:)
INTEGER :: I, J, K
INTEGER :: NID, IUNIT, N, JORG
REAL, PARAMETER :: ERR = 1.E-25_8
REAL :: Log2, RAMBDA

CALL OpenFile(InDev, .FALSE., .FALSE., .FALSE.,FileName)
READ(InDev, *)
READ(InDev, *) NISO1, NISO2, NFISYLD
NISODEP = NISO1 + NISO2
DeplLib%NISODEP = NISODEP; DeplLib%NFISYLD = NFISYLD
ALLOCATE(NIDINP(NISODEP))

!READ Fission Yield
ALLOCATE(DeplLib%FisYield(NFISYLD))
FisYield => DeplLib%FisYield
READ(InDev, *) (FisYield(i)%MatID, i = 1, NFISYLD)
READ(InDev, *)
NATOMDEP = 0
DO I = 1 , NISO1
  READ(InDev, *) NIDINP(I)
  READ(InDev, *)
  IZ = NIDINP(I)/10000
  NATOMDEP = MAX(IZ, NATOMDEP)
ENDDO
DeplLib%NATOMDEP = NATOMDEP

READ(InDev, *)
DO I = NISO1 + 1, NISODEP
  READ(InDev, *) NIDINP(I)
  IZ = NIDINP(I)/10000
  NATOMDEP = MAX(IZ, NATOMDEP)  
ENDDO

ALLOCATE(DeplLib%AtomLib0(NATOMDEP))
AtomLib0 => DeplLib%AtomLib0
DO I = 1, NATOMDEP
  AtomLib0%IB = 1000; AtomLib0%IE = 0
ENDDO

!DETERMINE ISOTOPE For Each Atom
DO I = 1, NISODEP
  IZ = NIDINP(I) / 10000
  IA = MOD(NIDINP(I), 10000)/10
  AtomLib0(IZ)%IB = MIN(AtomLib0(IZ)%IB, IA)
  AtomLib0(IZ)%IE = MAX(AtomLib0(IZ)%IE, IA)
ENDDO

DO I = 1, NATOMDEP
  IF(AtomLib0(I)%IE .EQ. 0) CYCLE
  ALLOCATE(AtomLib0(I)%A(AtomLib0(I)%IB:AtomLib0(I)%IE))
  DO J = AtomLib0(I)%IB, AtomLib0(I)%IE
    AtomLib0(I)%A(J)%NSTATE = -1
  ENDDO
ENDDO
!Determine # of states
DO I = 1, NISODEP
  IZ = NIDINP(I) / 10000
  IA = MOD(NIDINP(I), 10000)/10 
  IST = MOD(NIDINP(I), 10)
  AtomLib0(IZ)%A(IA)%NSTATE = MAX(AtomLib0(IZ)%A(IA)%NSTATE, IST)
ENDDO

!Initialize the 
DO I = 1, NATOMDEP
  DO J = AtomLib0(i)%IB, AtomLib0(i)%IE
    IF(AtomLib0(I)%A(J)%NSTATE .EQ. -1) CYCLE
    ALLOCATE(AtomLib0(I)%A(J)%STAT(0:AtomLib0(I)%A(J)%NSTATE))
    DO K = 0, AtomLib0(I)%A(J)%NSTATE
      AtomLib0(I)%A(J)%STAT(K)%IGRP = 1; AtomLib0(I)%A(J)%STAT(K)%NTO1 = 0
      AtomLib0(I)%A(J)%STAT(K)%NTO2 = 0; AtomLib0(I)%A(J)%STAT(K)%NFR3 = 0
      AtomLib0(I)%A(J)%STAT(K)%IMAT = 0; AtomLib0(I)%A(J)%STAT(K)%RAMBDA = 0
      AtomLib0(I)%A(J)%STAT(K)%FRAC = 0; AtomLib0(I)%A(J)%STAT(K)%XS = 0    
    ENDDO
  ENDDO
ENDDO

REWIND(InDev)


DO I = 1, 4
  READ(InDev, *)
ENDDO

Log2 = LOG(2._8); IMAT = 0
DO I = 1, NISO1
  IZ = NIDINP(I) / 10000
  IA = MOD(NIDINP(I), 10000)/10 
  IST = MOD(NIDINP(I), 10)
  IF(IZ .GE. 90) ATOMLIB0(iz)%A(IA)%STAT(IST)%IGRP = 2
  IMAT = IMAT + 1
  AtomLib0(IZ)%A(IA)%STAT(IST)%IMAT = IMAT
  !READ DECAY DATA
  READ(InDev, *) NID, IUNIT, RAMBDA, (AtomLib0(iz)%A(IA)%STAT(IST)%FRAC(K), k= 1, 7)
  IF(RAMBDA .GT. ERR) AtomLib0(IZ)%A(IA)%STAT(IST)%RAMBDA=LOG2/(RAMBDA*TCONV(IUNIT))
  AtomLib0(iz)%A(IA)%STAT(IST)%FRAC(8) = 1._8 - SUM(AtomLib0(iz)%A(IA)%STAT(IST)%FRAC(1:7))
  IF(AtomLib0(iz)%A(IA)%STAT(IST)%FRAC(8) .LT. 0._8) AtomLib0(iz)%A(IA)%STAT(IST)%FRAC(8) = 0
  !One- Group Xs DATA
  
  READ(InDev, *) (AtomLib0(iz)%A(IA)%STAT(IST)%XS(K), K = 1, 6), AtomLib0(iz)%A(IA)%STAT(IST)%kappa
  AtomLib0(iz)%A(IA)%STAT(IST)%XS(0) = SUM(AtomLib0(iz)%A(IA)%STAT(IST)%XS(1:6))
ENDDO

READ(InDev, *)
!FISSION PRODUCT YIELD DATA
DO I = NISO1 + 1, NISODEP
  IZ = NIDINP(I) / 10000
  IA = MOD(NIDINP(I), 10000)/10 
  IST = MOD(NIDINP(I), 10)
  AtomLib0(IZ)%A(IA)%STAT(IST)%IGRP = 3
  IMAT = IMAT +1
  AtomLib0(IZ)%A(IA)%STAT(IST)%IMAT = IMAT
  CALL DMALLOC(AtomLib0(IZ)%A(IA)%STAT(IST)%Y, NFISYLD)
  READ(InDev, *) N, (AtomLib0(IZ)%A(IA)%STAT(IST)%Y(K), K = 1, NFISYLD)
  AtomLib0(IZ)%A(IA)%STAT(IST)%Y = 1.E-2_8 * AtomLib0(IZ)%A(IA)%STAT(IST)%Y
ENDDO

DO I = 1, NFISYLD
  FisYield(I)%AtomNum = FisYield(I)%MatID / 10000
  FisYield(I)%AtomWeight = MOD(FisYield(I)%MatID, 10000)/10
  FisYield(I)%Stat = MOD(FisYield(I)%MatID, 10)
  FisYield(I)%MatID = AtomLib0(FisYield(I)%AtomNum)%A(FisYield(I)%AtomWeight)%STAT(FisYield(I)%Stat)%IMAT
ENDDO

!      DO I=1,NFISYLD
!        IFISYD(1,I)=IFISYD(0,I)/1000                             !! Atomic number
!        IFISYD(2,I)=MOD(IFISYD(0,I),1000)                        !! Atomic weight
!        IFISYD(0,I)=Z0(IFISYD(1,I))%A(IFISYD(2,I))%STAT(0)%IMAT  !! matrix ID
!      ENDDO

DO I = 1, NATOMDEP
  IF(AtomLib0(I)%IE .LT. 500) CYCLE
  DO J = 500, AtomLib0(I)%IE
    IF(AtomLib0(I)%A(J)%NSTATE .EQ. -1) CYCLE
    DO K = 0, AtomLib0(I)%A(J)%NSTATE
      IF(AtomLib0(I)%A(J)%STAT(K)%IMAT .EQ. 0) CYCLE
      JORG = J - 500
      AtomLib0(I)%A(J)%STAT(K)%XS(0:6) = AtomLib0(I)%A(JORG)%STAT(K)%XS(0:6)
      AtomLib0(I)%A(J)%STAT(K)%FRAC(1:8) = AtomLib0(I)%A(JORG)%STAT(K)%FRAC(1:8)
      AtomLib0(I)%A(J)%STAT(K)%RAMBDA = AtomLib0(I)%A(JORG)%STAT(K)%RAMBDA
      CONTINUE
    ENDDO
  ENDDO
ENDDO
DEALLOCATE(NIDINP)
CLOSE(InDeV)
NULLIFY(AtomLib0);NULLIFY(FisYield)
END SUBROUTINE