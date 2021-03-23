#include <Depletion.h>
#ifdef __PGI
  MODULE CUDA_Workspace
    USE CUDAFOR
    USE CUDA_CONST, ONLY : rowptr, colIdx, eyeIdx, CRAM_ORDER, Pole, Res, Res0
    IMPLICIT NONE
    INTEGER, PARAMETER :: Grid = 64, Bloc = 256
!    INTEGER, PARAMETER :: Grid = 1, Bloc = 1

    LOGICAL :: lcsrT = .FALSE., lFMv = .TRUE., lGS = .TRUE., lILP = .FALSE.
    LOGICAL :: lAlloc = .FALSE., lMatInfo = .FALSE., lAllocWS = .FALSE.
    INTEGER :: NsysMax=0, Nsys=0
    INTEGER :: NNZMax=0, NRMax=0, NNZ=0, NR=0
    COMPLEX(8), ALLOCATABLE, DEVICE :: BatchMat(:,:), BatchDiag(:,:)
    REAL(8), ALLOCATABLE, DEVICE :: BatchOffDiag(:,:), BatchVec(:,:), BatchSolVec(:), trVec(:,:), BatchDiagReal(:,:)
    COMPLEX(8), ALLOCATABLE, DEVICE :: vec1(:,:), vec2(:,:), r(:,:), p(:,:), phat(:,:), s(:,:), shat(:,:), t(:,:), v(:,:), rtilde(:,:)

    INTERFACE CopyOutBatchMem
      MODULE PROCEDURE CopyOutBatchMem_CSR
      MODULE PROCEDURE CopyOutBatchMem_woCSR
    END INTERFACE

    CONTAINS

    SUBROUTINE WriteConstCRAM()
    IMPLICIT NONE
#ifdef CRAM_14
    CRAM_Order = 14;
    Pole = (/ (-8.897773186468888, 16.630982619902085), (-3.703275049423448, 13.656371871483268), &
      (-0.208758638250130, 10.991260561901260), (3.993369710578568, 6.004831642235037), &
      (5.089345060580624, 3.588824029027006), (5.623142572745977, 1.194069046343966), &
      (2.269783829231112, 8.461797973040221) /);
    Res = (/ (-7.154288063589067e-5, 1.436104334854130e-4), (9.439025310736168e-3, -1.718479195848301e-2)  ,&
      (-3.763600387822696e-1, 3.351834702945010e-1), (-2.349823209108270e+1, -5.808359129714207)    ,&
      (4.693327448883129e+1, 4.564364976882776e+1), (-2.787516194014564e+1, -1.0214733999015645e+2) ,&
      (4.807112098832508, -1.320979383742872) /);
    Res0 = 1.832174378254041e-14
#endif
#ifdef CRAM_16
    CRAM_Order = 16;
    Pole = (/(-1.0843917078696988026e1, 1.9277446167181652284e1), (-5.2649713434426468895, 1.6220221473167927305e1)  ,&
      (5.9481522689511774808, 3.5874573620183222829), (3.5091036084149180974, 8.4361989858843750826)          ,&
      (6.4161776990994341923, 1.1941223933701386874), (1.4193758971856659786, 1.0925363484496722585e1)        ,&
      (4.9931747377179963991, 5.9968817136039422260), (-1.4139284624888862114, 1.3497725698892745389e1)/);
    Res = (/(-5.0901521865224915650e-7, -2.4220017652852287970e-5), (2.1151742182466030907e-4, 4.3892969647380673918e-3) ,&
      (1.1339775178483930527e2, 1.0194721704215856450e2), (1.5059585270023467528e1, -5.7514052776421819979)        ,&
      (-6.4500878025539646595e1, -2.2459440762652096056e2), (-1.4793007113557999718, 1.7686588323782937906)        ,&
      (-6.2518392463207918892e1, -1.1190391094283228480e1), (4.1023136835410021273e-2, -1.5743466173455468191e-1)/);
    Res0 = 2.1248537104952237488e-16
#endif
    END SUBROUTINE

    SUBROUTINE SetMatInfo(nnzin, nrin, rowptrin, colIdxin)
    IMPLICIT NONE
    INTEGER :: nnzin, nrin
    INTEGER :: rowptrin(:), colIdxin(:)
    INTEGER :: i,j,k,l
    NNZ = nnzin; NR = nrin
    IF (NNZMax .LT. NNZ .OR. NRMax .LT. NR) THEN
      NNZMax = NNZ; NRMax = NR;
      IF (lAlloc) THEN
        IF (lcsrT) THEN
          DEALLOCATE(BatchMat, BatchVec, BatchSolVec)
          ALLOCATE(BatchMat(NsysMax, NNZ), BatchVec(NsysMax, NR), BatchSolVec(NsysMax*NR))
        ELSE
!          DEALLOCATE(BatchDiag,BatchOffDiag,BatchVec,BatchSolVec)
          DEALLOCATE(BatchDiagReal,BatchOffDiag,BatchVec,BatchSolVec)
!          ALLOCATE(BatchDiag(NsysMtchOffDiag,BatchVec,BatchSolVec)
          ALLOCATE(BatchDiagReal(NsysMax,NRMax),BatchOffDiag(NsysMax,NNZMax-NRMax),BatchVec(NR,NsysMax),BatchSolVec(NsysMax*NRMax))
!          print*, NsysMax, NRMax, NNZMax-NRMax
        END IF
      END IF
      IF (lAllocWs) THEN
        CALL DeleteBlockMem()
        CALL AllocBlockMem()
      END IF
    END IF

    IF (lcsrT) THEN
      DO i = 1, NR
        rowptr(i) = rowptrin(i)
        DO j = rowptrin(i),rowptrin(i+1)-1
          IF (colIdxin(j) .EQ. i) THEN
            eyeIdx(i) = j; CYCLE
          END IF
        END DO
      END DO
      rowptr(NR+1) = rowptrin(NR+1)
      DO i = 1, NNZ
        colIdx(i) = colIdxin(i)
      END DO
    ELSE
      l = 0;
      DO i = 1, NR
        rowptr(i) = rowptrin(i)-i+1
        DO j = rowptrin(i),rowptrin(i+1)-1
          k = colIdxin(j)
          IF (k .NE. i) THEN
            l = l+1
            colIdx(l)=k
          END IF
        END DO
      END DO
      rowptr(NR+1) = rowptrin(NR+1)-NR
    END IF

    lMatInfo = .TRUE.
    END SUBROUTINE

    SUBROUTINE AllocBlockMem()
    IMPLICIT NONE
    IF (lAllocWS) RETURN
    IF (lGS) THEN
      ALLOCATE(trVec(Grid*Bloc,NRMax))
!      ALLOCATE(vec1(Grid*Bloc,NRMax), vec2(Grid*Bloc,NRMax))
      IF (lILP) THEN
#ifdef CRAM_14
        ALLOCATE(vec1(Grid*Bloc,7*NRMax))
#endif
#ifdef CRAM_16
        ALLOCATE(vec1(Grid*Bloc,8*NRMax))
#endif
      ELSE
        ALLOCATE(vec1(Grid*Bloc,NRMax))
!        ALLOCATE(vec1(Grid*Bloc,NRMax), vec2(Grid*Bloc,NRMax))
      END IF
!      print*, NRMax, Grid*Bloc
    ELSE
      ALLOCATE(vec2(Grid*Bloc,NRMax), r(Grid*Bloc,NRMax), p(Grid*Bloc,NRMax), phat(Grid*Bloc,NRMax))
      ALLOCATE(s(Grid*Bloc,NRMax), shat(Grid*Bloc,NRMax), t(Grid*Bloc,NRMax), v(Grid*Bloc,NRMax), rtilde(Grid*Bloc,NRMax))
    END IF
    lAllocWS = .TRUE.
    END SUBROUTINE

    SUBROUTINE AllocBatchMem(NsysIn)
    IMPLICIT NONE
    INTEGER :: NsysIn
    IF (NsysMax .LT. NsysIn) THEN
      NsysMax = NsysIn
      IF (lcsrT) THEN
        IF (lAlloc) THEN
          DEALLOCATE(BatchMat, BatchVec, BatchSolVec)
        END IF
        ALLOCATE(BatchMat(NsysIn, NNZ), BatchVec(NsysIn, NR), BatchSolVec(NsysIn*NR))
!        ALLOCATE(BatchMat(NNZ,NsysIn), BatchVec(NR,NsysIn), BatchSolVec(NsysIn*NR))
      ELSE
        IF (lAlloc) THEN
!          DEALLOCATE(BatchDiag,BatchOffDiag,BatchVec,BatchSolVec)
          DEALLOCATE(BatchDiagReal,BatchOffDiag,BatchVec,BatchSolVec)
        END IF
!        ALLOCATE(BatchDiag(NsysIn,NR), BatchOffDiag(NsysIn,NNZ-NR), BatchVec(NR,NsysIn), BatchSolVec(NsysIn*NR))
        ALLOCATE(BatchDiagReal(NsysIn,NR), BatchOffDiag(NsysIn,NNZ-NR), BatchVec(NR,NsysIn), BatchSolVec(NsysIn*NR))
!        print*, NsysIn, NR, NNZ-NR
      END IF
    ELSE IF (.NOT. lAlloc) THEN
      IF (lcsrT) THEN
        ALLOCATE(BatchMat(NsysMax, NNZ), BatchVec(NsysMax, NR), BatchSolVec(NsysMax*NR))
      ELSE
!        ALLOCATE(BatchDiag(NsysMax,NR), BatchOffDiag(NsysMax,NNZ-NR), BatchVec(NR,NsysMax), BatchSolVec(NsysMax*NR))
        ALLOCATE(BatchDiagReal(NsysMax,NR), BatchOffDiag(NsysMax,NNZ-NR), BatchVec(NR,NsysMax), BatchSolVec(NsysMax*NR))
!        print*, NsysMax, NR, NNZ-NR
      END IF
    END IF
    Nsys = NsysIn
    lAlloc = .TRUE.
    END SUBROUTINE

    SUBROUTINE CopyOutBatchMem_CSR(batchvecin,batchmatin)
    IMPLICIT NONE
    REAL(8) :: batchvecin(:,:)
    COMPLEX(8) :: batchmatin(:,:)
    BatchVec = batchvecin
    BatchMat = batchmatin
    END SUBROUTINE

    SUBROUTINE CopyOutBatchMem_woCSR(batchvecin,batchdiagin,batchoffdiagin)
    IMPLICIT NONE
    REAL(8) :: batchvecin(:,:)
!    COMPLEX(8) :: batchdiagin(:,:)
    REAL(8) :: batchdiagin(:,:)
    REAL(8) :: batchoffdiagin(:,:)
    INTEGER :: istat

!    print*, size(batchvecin,1), size(batchvecin,2)
!    print*, size(batchdiagin,1), size(batchdiagin,2)
!    print*, size(batchoffdiagin,1), size(batchoffdiagin,2)
!    BatchVec = batchvecin
!    BatchDiag = batchdiagin
!    BatchOffDiag = batchoffdiagin
!    istat = cudaMemcpy(BatchVec,batchvecin,NR*Nsys)
!    istat = cudaMemcpy(BatchDiag,batchdiagin,Nsys*NR)
!    istat = cudaMemcpy(BatchOffDiag,batchoffdiagin,Nsys*(NNZ-NR))
    istat = cudaMemcpy2D(BatchVec(1,1),NR,batchvecin(1,1),NR,NR,Nsys)
!    istat = cudaMemcpy2D(BatchDiag(1,1),NsysMax,batchdiagin(1,1),Nsys,Nsys,NR)
    istat = cudaMemcpy2D(BatchDiagReal(1,1),NsysMax,batchdiagin(1,1),Nsys,Nsys,NR)
    istat = cudaMemcpy2D(BatchOffDiag(1,1),NsysMax,batchoffdiagin(1,1),Nsys,Nsys,(NNZ-NR))
    END SUBROUTINE

    SUBROUTINE DeleteBatchMem()
    IMPLICIT NONE
    IF (lcsrT) THEN
      DEALLOCATE(BatchMat, BatchVec, BatchSolVec)
    ELSE
!      DEALLOCATE(BatchDiag, BatchOffDiag, BatchVec, BatchSolVec)
      DEALLOCATE(BatchDiagReal, BatchOffDiag, BatchVec, BatchSolVec)
    END IF
    lAlloc = .FALSE.
    END SUBROUTINE

    SUBROUTINE DeleteBlockMem()
    IMPLICIT NONE
    IF (lGS) THEN
      DEALLOCATE(trVec)
!      DEALLOCATE(vec1,vec2)
      DEALLOCATE(vec1)
    ELSE
      DEALLOCATE(vec2, r, p, phat, s, shat, t, v, rtilde)
    END IF
    lAllocWS = .FALSE.
    END SUBROUTINE
  END MODULE
#endif
