! Routines for define depletion chains
#include <Depletion.h>
SUBROUTINE DeplLibInit(libIO, Filename, DeplLib)
  USE HPDeplType
  IMPLICIT NONE
  INTERFACE
    SUBROUTINE DefineDeplChainType(DeplLib)
    USE HPDeplType
    IMPLICIT NONE
    TYPE(DeplLib_Type) :: DeplLib
    END SUBROUTINE DefineDeplChainType

    INTEGER FUNCTION IsoidSrch(IsoArr, NofIso, Idiso)
    USE HPDeplType
    INTEGER, INTENT(IN) :: NofIso, Idiso
    INTEGER, INTENT(IN) :: IsoArr(*)
    END FUNCTION
  END INTERFACE
  INTEGER, INTENT(IN) :: libIO
  CHARACTER(*) :: Filename
  TYPE(DeplLib_Type) :: DeplLib
  CHARACTER*256 :: Line
  INTEGER :: i, j, N, dum1, dum2, idf
  REAL(8) :: lambda
  INTEGER :: NofIso, NofYld, NofFis
  INTEGER :: NofRct, NofDec, IdLast
  REAL(8), POINTER :: DecFrac(:), RctXs(:), FisFrac(:)
  INTEGER, ALLOCATABLE :: delIdRct(:), delIdDec(:), IdAftRct(:), IdAftDec(:)
  REAL(8), PARAMETER :: TCONV(9) = (/1._8, 60._8, 3600._8, 86400._8, 3.15576E+7_8, &
                                1._8, 3.15576E+10_8, 3.15576E+13_8, 3.15576E+16_8/)
  LOGICAL, ALLOCATABLE :: YldMapIsNZ(:,:)

  CALL DefineDeplChainType(DeplLib)
  NofRct = DeplLib%NofRct; NofDec = DeplLib%NofDec
  ALLOCATE(DecFrac(NofDec), RctXs(NofRct))
  OPEN(libIO, FILE = Filename)
  DO WHILE(.TRUE.)
    READ(libIO,*) Line
    IF (TRIM(Line) .EQ. 'GENINFO') EXIT
  END DO
  READ(libIO, *) NofIso, NofYld, NofFis
  N = NofIso
#ifdef nTRACER_Depl
  NofIso = NofIso+NofYld
#endif
  ALLOCATE(FisFrac(NofFis))
  DeplLib%NofIso = NofIso; DeplLib%NofYld = NofYld; DeplLib%NofFis = NofFis
  ALLOCATE(DeplLib%Idiso(NofIso), DeplLib%lActinide(NofIso))
  ALLOCATE(DeplLib%RctXs(NofRct, NofIso), DeplLib%DecFrac(NofDec, NofIso))
  ALLOCaTE(DeplLib%IdAftRct(NofRct, NofIso), DeplLib%IdAftDec(NofDec, NofIso))
  ALLOCATE(DeplLib%IdFis(NofFis), DeplLib%FisFrac(NofFis, NofYld), DeplLib%IdYld(NofYld))
  ALLOCATE(DeplLIb%Kappa(NofIso))

  ALLOCATE(DeplLib%IzDec(NofDec,NofIso), DeplLib%IzRct(NofRct,NofIso))
  ALLOCATE(DeplLib%IzSYD(NofDec,NofIso), DeplLib%IzSYR(NofRct,NofIso), DeplLib%IzSYRAct(NofRct,NofIso))
  ALLOCATE(DeplLib%IzFY(NofFis,NofYld))

  DeplLib%lActinide(:) = .FALSE.
  READ(libIO, *) DeplLib%IdFis(:)

  DO WHILE(.TRUE.)
    READ(libIO,*) Line
    IF (TRIM(Line) .EQ. 'ISOTOPE') EXIT
  END DO
  RctXs = 0.; DecFrac = 0.
  DeplLib%RctXs(:,:) = 0.; DeplLib%DecFrac(:,:) = 0.; DeplLib%Kappa(:) = 0.
  DO i = 1, N
    READ(libIO, *) dum1, dum2, DecFrac(:)
    READ(libIO, *) RctXs(1:NofRct), DeplLib%Kappa(i)
    lambda = DecFrac(1)
    IF(lambda .NE. 0.) lambda = log(2.)/lambda/TCONV(dum2)
    DecFrac(1:7) = DecFrac(2:8)*lambda
    DecFrac(8) = lambda - SUM(DecFrac(1:7))
    IF(DecFrac(8) .LT. 0.) DecFrac(8) = 0.
    DeplLib%RctXs(1:NofRct,i) = RctXs(:)
    DeplLib%DecFrac(1:NofDec,i) = DecFrac(:)

    IF (dum1/10000 .GE. 90) THEN
      DeplLib%lActinide(i) = .TRUE.
    ELSE
      DeplLib%lActinide(i) = .FALSE.
    END IF
    DeplLib%Idiso(i) = dum1
  END DO
  IdLast = dum1

  DO WHILE(.TRUE.)
    READ(libIO,*) Line
    IF (TRIM(Line) .EQ. 'FISPROD') EXIT
  END DO
  FisFrac = 0.
  DO i = 1, NofYld
    READ(libIO, *) dum1, FisFrac(:)
    DeplLib%FisFrac(:,i) = FisFrac(:)
    DeplLib%IdYld(i) = dum1
#ifdef nTRACER_Depl
    DeplLib%Idiso(i+N) = dum1
    dum2 = Isoidsrch(DeplLib%IdIso, NofIso, dum1-5000)
    IF (dum2 .EQ. 0) CYCLE
    DeplLib%DecFrac(:,i+N) = DeplLib%DecFrac(:,dum2)
    DeplLib%RctXs(:,i+N) = DeplLib%RctXs(:,dum2)
#endif
  END DO
  DeplLib%FisFrac(:,:) = DeplLib%FisFrac*1.0e-2
  CLOSE(libIO)

  ALLOCATE(delIdRct(NofRct), delIdDec(NofDec))
  DO i = 1, NofIso
    dum1 = DeplLib%Idiso(i)
    IF (DeplLib%lActinide(i)) THEN
      delIdRct(:) = DeplLib%delIdRctAct(:)
    ELSE
      delIdRct(:) = DeplLib%delIdRct(:)
    END IF
    delIdDec(:) = DeplLib%delIdDec(:)
    DO j = 1, NofRct
      IF (mod(DelIdRct(j),10) .NE. 9) THEN
        dum2 = (dum1/10)*10+DelIdRct(j)
      ELSE
        dum2 = dum1 + delIdRct(j)
      END IF
      IF (dum2 .EQ. dum1) dum2 = 0
      DeplLib%IdAftRct(j,i) = dum2
    END DO
    DO j = 1, NofDec
      IF (mod(DelIdDec(j),10) .NE. 9) THEN
        dum2 = (dum1/10)*10+delIdDec(j)
      ELSE
        dum2 = dum1 + delIdDec(j)
      END IF
      IF (dum2 .EQ. dum1) dum2 = 0
      DeplLib%IdAftDec(j,i) = dum2
    END DO
  END DO
  DEALLOCATE(delIdRct, delIdDec)
  ALLOCATE(IdAftRct(NofRct), IdAftDec(NofDec))
  DO i = 1, NofIso
    IdAftRct(:) = DeplLib%IdAftRct(:,i)
    IdAftDec(:) = DeplLib%IdAftDec(:,i)
    DO j = 1, NofRct
      dum1 = IdAftRct(j)
      dum2 = IsoidSrch(DeplLib%Idiso(:), NofIso, dum1)
      IdAftRct(j) = dum2
    END DO
    DO j = 1, NofDec
      dum1 = IdAftDec(j)
      dum2 = IsoidSrch(DeplLib%Idiso(:), NofIso, dum1)
      IdAftDec(j) = dum2
    END DO
    DeplLib%IdAftRct(:,i) = IdAftRct(:)
    DeplLib%IdAftDec(:,i) = IdAftDec(:)
  END DO
  DO i = 1, NofYld
    DeplLib%IdYld(i) = IsoidSrch(DeplLib%Idiso(:), NofIso, DeplLib%IdYld(i))
  END DO
  DO i = 1, NofRct
    DeplLib%SubYldRct(i) = IsoidSrch(DeplLib%Idiso(:), NofIso, DeplLib%SubYldRct(i))
    DeplLib%SubYldRctAct(i) = IsoidSrch(DeplLib%Idiso(:), NofIso, DeplLib%SubYldRctAct(i))
  END DO
  DO i = 1, NofDec
    DeplLib%SubYldDec(i) = IsoidSrch(DeplLib%Idiso(:), NofIso, DeplLib%SubYldDec(i))
  END DO
  DO i = 1, NofFis
    DeplLib%IdFis(i) = IsoidSrch(DeplLib%Idiso(:), NofIso, DeplLib%IdFis(i))
  END DO
  DEALLOCATE(DecFrac, RctXs)

  ALLOCATE(YldMapIsNZ(NofIso, NofIso))
  YldMapIsNZ(:,:) = .FALSE.
  DO i = 1, NofFis
    DO j = 1, NofYld
      YldMapIsNZ(DeplLib%IdYld(j), DeplLib%IdFis(i)) = .TRUE.
    END DO
  END DO
  DO i = 1, NofIso
    YldMapIsNZ(i,i) = .TRUE.
    RctXs => DeplLib%RctXs(:,i); DecFrac => DeplLib%DecFrac(:,i)
    IdAftRct(:) = DeplLib%IdAftRct(:,i); IdAftDec(:) = DeplLib%IdAftDec(:,i)
    DO j = 1, NofRct
      IF (RctXs(j) .LT. 1.e-30) CYCLE
      dum1 = IdAftRct(j)
      IF (dum1 .NE. 0) THEN
        YldMapIsNZ(dum1, i) = .TRUE.
      END IF
      IF (DeplLib%lActinide(i)) THEN
        dum1 = DeplLib%SubYldRctAct(j)
      ELSE
        dum1 = DeplLib%SubYldRct(j)
      END IF
      IF (dum1 .NE. 0) THEN
        YldMapIsNZ(dum1, i) = .TRUE.
      END IF
    END DO
    DO j = 1, NofDec
      IF (DecFrac(j) .LT. 1.e-30) CYCLE
      dum1 = IdAftDec(j)
      IF (dum1 .NE. 0) THEN
        YldMapIsNZ(dum1, i) = .TRUE.
      END IF
      dum1 = DeplLib%SubYldDec(j)
      IF (dum1 .NE. 0) THEN
        YldMapIsNZ(dum1, i) = .TRUE.
      END IF
    END DO
  END DO

  dum1 = 0
  DO i = 1, NofIso
    DO j = 1, NofIso
      IF (YldMapIsNZ(i,j)) dum1 = dum1 + 1
    END DO
  END DO
  ALLOCATE(DeplLib%YldMapColIdx(dum1), DeplLib%YldMapRowPtr(NofIso+1))
  dum1 = 0
  DO i = 1, NofIso
    DeplLib%YldMapRowPtr(i) = dum1+1
    Do j = 1, NofIso
      IF (YldMapIsNZ(i,j)) THEN
        dum1 = dum1 +1
        DeplLib%YldMapColIdx(dum1) = j
      END IF
    END DO
  END DO
  DeplLib%YldMapRowPtr(NofIso+1) = dum1+1

  DeplLib%IzDec=0; DeplLib%IzRct=0
  DeplLib%IzSYR=0; DeplLib%IzSYD=0; DeplLib%IzSYRAct=0
  DeplLib%IzFY=0

  DO i = 1, NofIso
    DO j = 1, NofDec
      dum1 = DeplLib%IdAftDec(j,i)
      IF (dum1.GT.0) THEN
        DO dum2 = DeplLib%YldMapRowPtr(dum1), DeplLib%YldMapRowPtr(dum1+1)-1
          IF (DeplLib%YldMapColIdx(dum2) .NE. i) CYCLE
          DeplLib%IzDec(j,i) = dum2
          EXIT
        END DO
      END IF
      dum1 = DeplLib%SubYldDec(j)
      IF (dum1.GT.0) THEN
        DO dum2 = DeplLib%YldMapRowPtr(dum1), DeplLib%YldMapRowPtr(dum1+1)-1
          IF (DeplLib%YldMapColIdx(dum2) .NE. i) CYCLE
          DeplLib%IzSYD(j,i) = dum2
          EXIT
        END DO
      END IF
    END DO

    DO j = 1, NofRct
      dum1 = DeplLib%IdAftRct(j,i)
      IF (dum1.GT.0) THEN
        DO dum2 = DeplLib%YldMapRowPtr(dum1), DeplLib%YldMapRowPtr(dum1+1)-1
          IF (DeplLib%YldMapColIdx(dum2) .NE. i) CYCLE
          DeplLib%IzRct(j,i) = dum2
          EXIT
        END DO
      END IF
      IF (.NOT.DeplLib%lActinide(i)) THEN
        dum1 = DeplLib%SubYldRct(j)
        IF (dum1.GT.0) THEN 
          DO dum2 = DeplLib%YldMapRowPtr(dum1), DeplLib%YldMapRowPtr(dum1+1)-1
            IF (DeplLib%YldMapColIdx(dum2) .NE. i) CYCLE
            DeplLib%IzSYR(j,i) = dum2
            EXIT
          END DO
        END IF
      ELSE
        dum1 = DeplLib%SubYldRctAct(j)
        IF (dum1 .GT. 0) THEN
          DO dum2 = DeplLib%YldMapRowPtr(dum1), DeplLib%YldMapRowPtr(dum1+1)-1
            IF (DeplLib%YldMapColIdx(dum2) .NE. i) CYCLE
            DeplLib%IzSYRAct(j,i) = dum2
            EXIT
          END DO
        END IF 
      END IF
    END DO
  END DO

  DO idf = 1, NofFis
    i = DeplLib%IdFis(idf)
    DO j = 1, NofYld
      dum1 = DeplLib%IdYld(j)
      IF (dum1.GT.0) THEN
        DO dum2 = DeplLib%YldMapRowPtr(dum1), DeplLib%YldMapRowPtr(dum1+1)-1
          IF (DeplLib%YldMapColIdx(dum2) .NE. i) CYCLE
          DeplLib%IzFY(idf,j) = dum2
          EXIT
        END DO
      END IF
    END DO
  END DO

  DEALLOCATE(YldMapIsNZ, IdAftDec, IdAftRct)
  NULLIFY(DecFrac, RctXs)
END SUBROUTINE DeplLibInit
  SUBROUTINE DefineDeplChainType(DeplLib)
  USE HPDeplType
  IMPLICIT NONE
  TYPE(DeplLib_Type) :: DeplLib

  DeplLib%NofRct = 6; DeplLib%NofDec = 8; DeplLib%FisRctId = 4
  ALLOCATE(DeplLib%delIdRct(6), DeplLib%delIdRctAct(6), DeplLib%delIdDec(8))
  ALLOCATE(DeplLib%SubYldRct(6), DeplLib%SubYldRctAct(6), DeplLib%SubYldDec(8))
  DeplLib%delIdRct = (/10, -10, -20030, -10000, 11, -10+1/)
  DeplLib%delIdRctAct = (/10, -10, -20, 0, 11, -10+1/)
  DeplLib%delIdDec = (/10001, -10000, -10000+1, -20040, -1, 0, 10000-10, 10000/)
  DeplLib%SubYldRct = (/0, 0, 20040, 10010, 0, 0/)
  DeplLib%SubYldRctAct = (/0, 0, 0, 0, 0, 0/)
  DeplLib%SubYldDec = (/0, 0, 0, 20040, 0, 0, 0, 0/)
END SUBROUTINE

INTEGER FUNCTION IsoidSrch(IsoArr, NofIso, Idiso)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NofIso, Idiso
  INTEGER, INTENT(IN) :: IsoArr(*)
  INTEGER :: i

  IsoidSrch = 0
  DO i = 1, NofIso
    IF(IsoArr(i) .EQ. Idiso) THEN
      IsoidSrch = i; EXIT
    END IF
  END DO
END FUNCTION
