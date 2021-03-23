#include <Depletion.h>

#ifdef CRAM_16
#define CRAM_ord 8
#endif
#ifdef CRAM_14
#define CRAM_ord 7
#endif

#ifdef __PGI
MODULE cuMatExponential
  USE CSRMATRIX
  USE CUDA_MASTER
  IMPLICIT NONE

  CONTAINS

  SUBROUTINE MatExpCRAM_Iter(x1,stream)
  !USE AuxilCSR
  USE CUDA_Workspace
  IMPLICIT NONE
  REAL(8), POINTER :: x1(:)
  INTEGER(KIND = cuda_stream_kind) :: stream

  INTEGER :: i, j, k

  CALL AllocBlockMem()

!  CALL cuCRAM_BiCGSTAB<<<Grid, Bloc, stream>>>(BatchMat, BatchVec, BatchSolVec, &
!    vec2, r, p, phat, s, shat, t, v, rtilde, NR, NNZ, Nsys)
!  CALL cuCRAM_BiCGSTAB_FM<<<Grid, Bloc, stream>>>(BatchMat, BatchVec, BatchSolVec, &
!    vec2, r, p, phat, s, shat, t, v, rtilde, NR, NNZ, Nsys)
!  CALL cuCRAM_BiCGSTAB_woCSR<<<Grid, Bloc, stream>>>(BatchDiag, BatchOffDiag, BatchVec, BatchSolVec, &
!    vec2, r, p, phat, s, shat, t, v, rtilde, NR, NNZ, Nsys)
!  CALL cuCRAM_GS<<<Grid, Bloc, stream>>>(BatchDiag, BatchOffDiag, BatchVec, BatchSolVec, &
!    trVec, vec1, vec2, NR, NNZ, Nsys)
!  CALL cuCRAM_GS<<<Grid, Bloc, stream>>>(BatchDiagReal, BatchOffDiag, BatchVec, BatchSolVec, &
!    trVec, vec1, vec2, NR, NNZ, Nsys)
  CALL cuCRAM_GS_LOW<<<Grid, Bloc, stream>>>(BatchDiagReal, BatchOffDiag, BatchVec, BatchSolVec, &
    trVec, vec1, NR, NNZ, Nsys)
!  CALL cuCRAM_GS_ILP<<<Grid, Bloc, stream>>>(BatchDiagReal, BatchOffDiag, BatchVec, BatchSolVec, &
!    trVec, vec1, vec2, NR, NNZ, Nsys)
  i = cudaMemcpyAsync(x1, BatchSolVec, Nsys*NR, cudaMemcpyDeviceToHost, stream)
  !i = cudaMemcpy(x1, BatchSolVec, Nsys*NR, cudaMemcpyDeviceToHost)
  !i = cudaStreamSynchronize(stream); print*, 'stat', i

  END SUBROUTINE

  ATTRIBUTES(GLOBAL) SUBROUTINE cuCRAM_BiCGSTAB(Mat, Vec, SolV, Vec2, r, p, phat, s, shat, t, v, rtilde, nr, nnz, nsys)
  USE CUDA_CONST, ONLY : rowptr, colIdx, eyeIdx, CRAM_Order, Pole, Res, Res0
  IMPLICIT NONE
  COMPLEX(8), DIMENSION(:,:), DEVICE :: Mat, Vec2, r, p, phat, s, shat, t, v, rtilde
  REAL(8), DIMENSION(:), DEVICE :: SolV
  REAL(8), DIMENSION(:,:), DEVICE :: Vec
  INTEGER, VALUE :: nr, nnz, Nsys

  INTEGER :: sizeBlk
  INTEGER :: tid, sysid
  INTEGER :: iord, iter, i, j, k, lb, ub
  INTEGER :: HalfOrder

  COMPLEX(8) :: polei, resi, temp, temp1
  COMPLEX(8) :: rho0, rho1, rho2, a, b, w

  INTEGER :: SizeSys

 ! REAL(4) :: Niter

  HalfOrder = CRAM_Order/2; SizeSys = Nsys

  sizeBlk = blockDim%x
  tid = (blockIdx%x-1)*sizeBlk+threadIdx%x
  sizeBlk = sizeBlk*GridDim%x
  sysid = tid
  !RETURN
  DO WHILE(.TRUE.)
    !IF (tid.EQ.1) WRITE(*,*) 'A', sysid, nsys
    IF (sysid .GT. nsys) EXIT
    !IF (tid.EQ.1) WRITE(*,*) 'B', sysid, Vec(sysid,343), Vec(sysid,350)
    DO i = 1, nr
      SolV(sysid+(i-1)*SizeSys) = Vec(sysid, i)*Res0
    END DO
    polei = 0.
    DO iord = 1, HalfOrder
      polei = Pole(iord)-polei; resi = Res(iord); rho1 = 0.;
      DO i = 1, nr
        Vec2(tid, i) = 0.
        j = eyeidx(i)
        Mat(sysid,j) = Mat(sysid,j)-polei
        temp = 2.*resi*Vec(sysid, i)
        r(tid,i) = temp; rtilde(tid,i) = temp
        p(tid,i) = temp;
        rho1 = rho1+conjg(temp)*temp;
      END DO
      polei = Pole(iord)
      rho0 = rho1;
      !IF(tid.EQ.1) WRITE(*,*) iord, rho0
      DO iter = 1, 10
        IF (iter.NE.1) THEN
          b = rho1/rho2*a/w
          DO i = 1, nr
            p(tid,i) = r(tid,i)+b*(p(tid,i)-w*v(tid,i))
          END DO
        END IF
        DO i = 1, nr
          j = eyeidx(i)
          phat(tid,i) = p(tid,i)/Mat(sysid,j)
        END DO
        temp = 0.;
        DO i = 1, nr
          lb = rowptr(i); ub = rowptr(i+1)-1
          !v(tid,i) = 0.;
          temp1 = 0.;
          DO j = lb, ub
            k = colIdx(j)
            temp1 = temp1+Mat(sysid,j)*phat(tid,k)
          END DO
          v(tid,i) = temp1;
          temp = temp + conjg(rtilde(tid,i))*v(tid,i)
        END DO
        a = rho1/temp
        !IF(tid.EQ.1) WRITE(*,*) 'TEMP', iord, iter, p(tid,350), Mat(sysid,eyeidx(350))!phat(tid, 350), v(tid,350)
        rho2 = 0.
        DO i = 1, nr
          s(tid,i) = r(tid,i)-a*v(tid,i)
          rho2 = rho2 + conjg(rtilde(tid,i))*s(tid,i)
        END DO
        IF (.NOT.(ABS(rho2) .LT. 0. .OR. ABS(rho2) .GE. 0.)) THEN
          !IF(tid.EQ.1) WRITE(*,*) 'Diverge 1'
          !Niter = iter-1
          EXIT
        END IF
        DO i = 1, nr
          Vec2(tid,i) = Vec2(tid,i)+a*phat(tid,i)
        eND DO
        IF (ABS(rho2/rho0).LT. 1.e-35) THEN
          !IF(tid.EQ.1) WRITE(*,*) 'Converge 1'
          !Niter = iter-0.5
          EXIT
        END IF
        ! ------------------- half iteration --------------------- !
        DO i = 1, nr
          j = eyeidx(i)
          shat(tid,i) = s(tid,i)/Mat(sysid,j)
        END DO
        temp = 0.; w = 0.;
        DO i = 1, nr
          lb = rowptr(i); ub = rowptr(i+1)-1
          temp1 = 0.; !t(tid,i) = 0.;
          DO j = lb, ub
            k = colIdx(j)
            temp1 = temp1+Mat(sysid,j)*shat(tid,k)
          END DO
          t(tid,i) = temp1
          w = w + conjg(temp1)*s(tid,i)
          temp = temp + conjg(temp1)*temp1
        END DO
        w = w/temp
        rho1 = 0.;
        DO i = 1, nr
          r(tid,i) = s(tid,i)-w*t(tid,i)
          rho1 = rho1 + conjg(rtilde(tid,i))*r(tid,i)
        END DO
        !IF (tid.EQ.1) WRITE(*,*) iord, iter, a, w, rho1, rho2
        rho2 = rho1
        IF (.NOT.(ABS(rho2) .LT. 0. .OR. ABS(rho2) .GE. 0.)) THEN
          !IF(tid.EQ.1) WRITE(*,*) 'Diverge 2'
          !Niter = iter-0.5
          EXIT
        END IF
        DO i = 1, nr
          Vec2(tid,i) = Vec2(tid,i)+w*shat(tid,i)
        END DO
        IF (ABS(rho2/rho0).LT. 1.e-35) THEN
          !IF(tid.EQ.1) WRITE(*,*) 'Converge 2'
          !Niter = iter
          EXIT
        END IF
      END DO
      !IF (tid.EQ.1) WRITE(*,*) 'Z', iord, sysid, iter
!      CALL syncthreads()
      DO i = 1, nr
        SolV(sysid+(i-1)*SizeSys) = SolV(sysid+(i-1)*SizeSys)+REAL(Vec2(tid,i))
      END DO
      !IF (tid.EQ.1) WRITE(*,*) 'V', sysid, iord, Niter, Vec2(tid,343), Vec2(tid,350)
    END DO
    !IF (tid.EQ.1) WRITE(*,*) 'S', sysid, SolV(sysid,343), SolV(sysid,350)
    sysid = sysid+sizeBlk
  END DO
  END SUBROUTINE

  ATTRIBUTES(GLOBAL) SUBROUTINE cuCRAM_BiCGSTAB_FM(Mat, Vec, SolV, Vec2, r, p, phat, s, shat, t, v, rtilde, nr, nnz, nsys)
  USE CUDA_CONST, ONLY : rowptr, colIdx, eyeIdx, CRAM_Order, Pole, Res, Res0
  IMPLICIT NONE
  COMPLEX(8), DIMENSION(:,:), DEVICE :: Mat, Vec2, r, p, phat, s, shat, t, v, rtilde
  REAL(8), DIMENSION(:), DEVICE :: SolV
  REAL(8), DIMENSION(:,:), DEVICE ::  Vec
  INTEGER, VALUE :: nr, nnz, Nsys

  INTEGER :: sizeBlk
  INTEGER :: tid, sysid
  INTEGER :: iord, iter, i, j, k, lb, ub
  INTEGER :: HalfOrder

  COMPLEX(8) :: polei, resi, temp, temp1
  COMPLEX(8) :: rho0, rho1, rho2, a, b, w

  INTEGER :: SizeSys

 ! REAL(4) :: Niter

  HalfOrder = CRAM_Order/2; SizeSys = Nsys

  sizeBlk = blockDim%x
  tid = (blockIdx%x-1)*sizeBlk+threadIdx%x
  sizeBlk = sizeBlk*GridDim%x
  sysid = tid
  !RETURN
  DO WHILE(.TRUE.)
    !IF (tid.EQ.1) WRITE(*,*) 'A', sysid, nsys
    IF (sysid .GT. nsys) EXIT
    !IF (tid.EQ.1) WRITE(*,*) 'B', sysid, Vec(sysid,343), Vec(sysid,350)
    DO i = 1, nr
!      SolV(sysid+(i-1)*SizeSys) = Vec(sysid, i)*Res0
      SolV(i+(sysid-1)*nr) = Vec(i,sysid)*Res0
    END DO
    polei = 0.
    DO iord = 1, HalfOrder
      polei = Pole(iord)-polei; resi = Res(iord); rho1 = 0.;
      DO i = 1, nr
        Vec2(tid, i) = 0.
        j = eyeidx(i)
!        Mat(sysid,j) = Mat(sysid,j)-polei
!        temp = 2.*resi*Vec(sysid, i)
        Mat(j,sysid) = Mat(j,sysid)-polei
        temp = 2.*resi*Vec(i,sysid)
        r(tid,i) = temp; rtilde(tid,i) = temp
        p(tid,i) = temp;
        rho1 = rho1+conjg(temp)*temp;
      END DO
      polei = Pole(iord)
      rho0 = rho1;
      !IF(tid.EQ.1) WRITE(*,*) iord, rho0
      DO iter = 1, 10
        IF (iter.NE.1) THEN
          b = rho1/rho2*a/w
          DO i = 1, nr
            p(tid,i) = r(tid,i)+b*(p(tid,i)-w*v(tid,i))
          END DO
        END IF
        DO i = 1, nr
          j = eyeidx(i)
!          phat(tid,i) = p(tid,i)/Mat(sysid,j)
          phat(tid,i) = p(tid,i)/Mat(j,sysid)
        END DO
        temp = 0.;
        DO i = 1, nr
          lb = rowptr(i); ub = rowptr(i+1)-1
          !v(tid,i) = 0.;
          temp1 = 0.;
          DO j = lb, ub
            k = colIdx(j)
!            temp1 = temp1+Mat(sysid,j)*phat(tid,k)
            temp1 = temp1+Mat(j,sysid)*phat(tid,k)
          END DO
          v(tid,i) = temp1;
          temp = temp + conjg(rtilde(tid,i))*v(tid,i)
        END DO
        a = rho1/temp
        !IF(tid.EQ.1) WRITE(*,*) 'TEMP', iord, iter, p(tid,350), Mat(sysid,eyeidx(350))!phat(tid, 350), v(tid,350)
        rho2 = 0.
        DO i = 1, nr
          s(tid,i) = r(tid,i)-a*v(tid,i)
          rho2 = rho2 + conjg(rtilde(tid,i))*s(tid,i)
        END DO
        IF (.NOT.(ABS(rho2) .LT. 0. .OR. ABS(rho2) .GE. 0.)) THEN
          !IF(tid.EQ.1) WRITE(*,*) 'Diverge 1'
          !Niter = iter-1
          EXIT
        END IF
        DO i = 1, nr
          Vec2(tid,i) = Vec2(tid,i)+a*phat(tid,i)
        eND DO
        IF (ABS(rho2/rho0).LT. 1.e-35) THEN
          !IF(tid.EQ.1) WRITE(*,*) 'Converge 1'
          !Niter = iter-0.5
          EXIT
        END IF
        ! ------------------- half iteration --------------------- !
        DO i = 1, nr
          j = eyeidx(i)
!          shat(tid,i) = s(tid,i)/Mat(sysid,j)
          shat(tid,i) = s(tid,i)/Mat(j,sysid)
        END DO
        temp = 0.; w = 0.;
        DO i = 1, nr
          lb = rowptr(i); ub = rowptr(i+1)-1
          temp1 = 0.; !t(tid,i) = 0.;
          DO j = lb, ub
            k = colIdx(j)
!            temp1 = temp1+Mat(sysid,j)*shat(tid,k)
            temp1 = temp1+Mat(j,sysid)*shat(tid,k)
          END DO
          t(tid,i) = temp1
          w = w + conjg(temp1)*s(tid,i)
          temp = temp + conjg(temp1)*temp1
        END DO
        w = w/temp
        rho1 = 0.;
        DO i = 1, nr
          r(tid,i) = s(tid,i)-w*t(tid,i)
          rho1 = rho1 + conjg(rtilde(tid,i))*r(tid,i)
        END DO
        !IF (tid.EQ.1) WRITE(*,*) iord, iter, a, w, rho1, rho2
        rho2 = rho1
        IF (.NOT.(ABS(rho2) .LT. 0. .OR. ABS(rho2) .GE. 0.)) THEN
          !IF(tid.EQ.1) WRITE(*,*) 'Diverge 2'
          !Niter = iter-0.5
          EXIT
        END IF
        DO i = 1, nr
          Vec2(tid,i) = Vec2(tid,i)+w*shat(tid,i)
        END DO
        IF (ABS(rho2/rho0).LT. 1.e-35) THEN
          !IF(tid.EQ.1) WRITE(*,*) 'Converge 2'
          !Niter = iter
          EXIT
        END IF
      END DO
      !IF (tid.EQ.1) WRITE(*,*) 'Z', iord, sysid, iter
!      CALL syncthreads()
      DO i = 1, nr
!        SolV(sysid+(i-1)*SizeSys) = SolV(sysid+(i-1)*SizeSys)+REAL(Vec2(tid,i))
        SolV(i+(sysid-1)*nr) = SolV(i+(sysid-1)*nr)+REAL(Vec2(tid,i))
      END DO
      !IF (tid.EQ.1) WRITE(*,*) 'V', sysid, iord, Niter, Vec2(tid,343), Vec2(tid,350)
    END DO
    !IF (tid.EQ.1) WRITE(*,*) 'S', sysid, SolV(sysid,343), SolV(sysid,350)
    sysid = sysid+sizeBlk
  END DO
  END SUBROUTINE

  ATTRIBUTES(GLOBAL) SUBROUTINE cuCRAM_BiCGSTAB_woCSR(Diag, OffDiag, Vec, SolV, Vec2, r, p, phat, s, shat, t, v, rtilde, nr, nnz, nsys)
  USE CUDA_CONST, ONLY : rowptr, colIdx, CRAM_Order, Pole, Res, Res0
  IMPLICIT NONE
  COMPLEX(8), DIMENSION(:,:), DEVICE :: Diag, Vec2, r, p, phat, s, shat, t, v, rtilde
  REAL(8), DIMENSION(:), DEVICE :: SolV
  REAL(8), DIMENSION(:,:), DEVICE :: OffDiag, Vec
  INTEGER, VALUE :: nr, nnz, Nsys

  INTEGER :: sizeBlk
  INTEGER :: tid, sysid
  INTEGER :: iord, iter, i, j, k, lb, ub
  INTEGER :: HalfOrder

  COMPLEX(8) :: polei, resi, temp, temp1
  COMPLEX(8) :: rho0, rho1, rho2, a, b, w

  INTEGER :: SizeSys

 ! REAL(4) :: Niter

  HalfOrder = CRAM_Order/2; SizeSys = Nsys

  sizeBlk = blockDim%x
  tid = (blockIdx%x-1)*sizeBlk+threadIdx%x
  sizeBlk = sizeBlk*GridDim%x
  sysid = tid
  !RETURN
  DO WHILE(.TRUE.)
    !IF (tid.EQ.1) WRITE(*,*) 'A', sysid, nsys
    IF (sysid .GT. nsys) EXIT
    !IF (tid.EQ.1) WRITE(*,*) 'B', sysid, Vec(sysid,343), Vec(sysid,350)
    DO i = 1, nr
      SolV(sysid+(i-1)*SizeSys) = Vec(i, sysid)*Res0
    END DO
    polei = 0.
    DO iord = 1, HalfOrder
      polei = Pole(iord)-polei; resi = Res(iord); rho1 = 0.;
      DO i = 1, nr
        Vec2(tid, i) = 0.
!        j = eyeidx(i)
        Diag(sysid,i) = Diag(sysid,i)-polei
        temp = 2.*resi*Vec(i, sysid)
        r(tid,i) = temp; rtilde(tid,i) = temp
        p(tid,i) = temp;
        rho1 = rho1+conjg(temp)*temp;
      END DO
      polei = Pole(iord)
      rho0 = rho1;
      !IF(tid.EQ.1) WRITE(*,*) iord, rho0
      DO iter = 1, 10
        IF (iter.NE.1) THEN
          b = rho1/rho2*a/w
          DO i = 1, nr
            p(tid,i) = r(tid,i)+b*(p(tid,i)-w*v(tid,i))
          END DO
        END IF
        DO i = 1, nr
!          j = eyeidx(i)
          phat(tid,i) = p(tid,i)/Diag(sysid,i)
        END DO
        temp = 0.;
        DO i = 1, nr
          lb = rowptr(i); ub = rowptr(i+1)-1
          !v(tid,i) = 0.;
          temp1 = Diag(sysid,i)*phat(tid,i);
          DO j = lb, ub
            k = colIdx(j)
            temp1 = temp1+OffDiag(sysid,j)*phat(tid,k)
          END DO
          v(tid,i) = temp1;
          temp = temp + conjg(rtilde(tid,i))*v(tid,i)
        END DO
        a = rho1/temp
        !IF(tid.EQ.1) WRITE(*,*) 'TEMP', iord, iter, p(tid,350), Mat(sysid,eyeidx(350))!phat(tid, 350), v(tid,350)
        rho2 = 0.
        DO i = 1, nr
          s(tid,i) = r(tid,i)-a*v(tid,i)
          rho2 = rho2 + conjg(rtilde(tid,i))*s(tid,i)
        END DO
        IF (.NOT.(ABS(rho2) .LT. 0. .OR. ABS(rho2) .GE. 0.)) THEN
          !IF(tid.EQ.1) WRITE(*,*) 'Diverge 1'
          !Niter = iter-1
          EXIT
        END IF
        DO i = 1, nr
          Vec2(tid,i) = Vec2(tid,i)+a*phat(tid,i)
        eND DO
        IF (ABS(rho2/rho0).LT. 1.e-35) THEN
          !IF(tid.EQ.1) WRITE(*,*) 'Converge 1'
          !Niter = iter-0.5
          EXIT
        END IF
        ! ------------------- half iteration --------------------- !
        DO i = 1, nr
!          j = eyeidx(i)
          shat(tid,i) = s(tid,i)/Diag(sysid,i)
        END DO
        temp = 0.; w = 0.;
        DO i = 1, nr
          lb = rowptr(i); ub = rowptr(i+1)-1
          temp1 = Diag(sysid,i)*shat(tid,i); !t(tid,i) = 0.;
          DO j = lb, ub
            k = colIdx(j)
            temp1 = temp1+OffDiag(sysid,j)*shat(tid,k)
          END DO
          t(tid,i) = temp1
          w = w + conjg(temp1)*s(tid,i)
          temp = temp + conjg(temp1)*temp1
        END DO
        w = w/temp
        rho1 = 0.;
        DO i = 1, nr
          r(tid,i) = s(tid,i)-w*t(tid,i)
          rho1 = rho1 + conjg(rtilde(tid,i))*r(tid,i)
        END DO
        !IF (tid.EQ.1) WRITE(*,*) iord, iter, a, w, rho1, rho2
        rho2 = rho1
        IF (.NOT.(ABS(rho2) .LT. 0. .OR. ABS(rho2) .GE. 0.)) THEN
          !IF(tid.EQ.1) WRITE(*,*) 'Diverge 2'
          !Niter = iter-0.5
          EXIT
        END IF
        DO i = 1, nr
          Vec2(tid,i) = Vec2(tid,i)+w*shat(tid,i)
        END DO
        IF (ABS(rho2/rho0).LT. 1.e-35) THEN
          !IF(tid.EQ.1) WRITE(*,*) 'Converge 2'
          !Niter = iter
          EXIT
        END IF
      END DO
      !IF (tid.EQ.1) WRITE(*,*) 'Z', iord, sysid, iter
!      CALL syncthreads()
      DO i = 1, nr
        SolV(sysid+(i-1)*SizeSys) = SolV(sysid+(i-1)*SizeSys)+REAL(Vec2(tid,i))
      END DO
      !IF (tid.EQ.1) WRITE(*,*) 'V', sysid, iord, Niter, Vec2(tid,343), Vec2(tid,350)
    END DO
    !IF (tid.EQ.1) WRITE(*,*) 'S', sysid, SolV(sysid,343), SolV(sysid,350)
    sysid = sysid+sizeBlk
  END DO
  END SUBROUTINE

  ATTRIBUTES(GLOBAL) SUBROUTINE cuCRAM_GS(Diag, OffDiag, Vec, SolV, trVec, Vec1, Vec2, nr, nnz, nsys)
  USE CUDA_CONST, ONLY : rowptr, colIdx, CRAM_Order, Pole, Res, Res0
  IMPLICIT NONE
!  COMPLEX(8), DIMENSION(:,:), DEVICE :: Diag, Vec1, Vec2
  COMPLEX(8), DIMENSION(:,:), DEVICE :: Vec1, Vec2
  REAL(8), DIMENSION(:), DEVICE :: SolV
  REAL(8), DIMENSION(:,:), DEVICE :: Diag, OffDiag, Vec, trVec
  INTEGER, VALUE :: nr, nnz, Nsys

  INTEGER :: sizeBlk
  INTEGER :: tid, sysid
  INTEGER :: iord, iter, i, j, k, lb, ub
  INTEGER :: HalfOrder

  COMPLEX(8) :: polei, resi, temp1!, val0
  REAL(8) :: norm, relerr, temp

  INTEGER :: SizeSys


  HalfOrder = CRAM_Order/2; SizeSys = Nsys

  sizeBlk = blockDim%x
  tid = (blockIdx%x-1)*sizeBlk+threadIdx%x
  sizeBlk = sizeBlk*GridDim%x
  sysid = tid
  DO WHILE(.TRUE.)
    IF (sysid .GT. nsys) EXIT
    DO i = 1, nr
      relerr = Vec(i,sysid) ! temporal for real variable
      trVec(tid,i) = relerr
      SolV(sysid+(i-1)*SizeSys) = relerr*Res0
    END DO
!    polei = 0.
    DO iord = 1, HalfOrder
      polei = Pole(iord); resi = Res(iord);
!      polei = Pole(iord)-polei; resi = Res(iord);
      DO i = 1, nr
        Vec1(tid,i) = trVec(tid,i)
!        Diag(sysid,i) = Diag(sysid,i)-polei
      END DO
!      polei = Pole(iord)
      DO iter = 1, 50
        DO i = 1, nr
          Vec2(tid, i) = Vec1(tid, i)
        END DO
!        norm = 0.; relerr = 0.
        DO i = 1, nr
!          val0 = Vec1(tid,i)
          lb = rowptr(i); ub = rowptr(i+1)-1
          temp1 = 2.*resi*trVec(tid,i)
          DO j = lb, ub
            k = colIdx(j)
            temp1 = temp1-OffDiag(sysid,j)*Vec1(tid,k)
          END DO
          temp1 = temp1/(Diag(sysid,i)-polei)
          Vec1(tid,i) = temp1
        END DO
        norm = 0.; relerr = 0.
        DO i = 1, nr
          temp = ABS(Vec1(tid, i)- Vec2(tid, i))
          relerr = relerr + temp*temp
          temp = ABS(Vec1(tid, i))
          norm = norm + temp*temp
!          Vec1(tid,i) = temp1
        END DO
        relerr = relerr/norm
        IF (relerr<1.e-16) THEN
          EXIT
        END IF
      END DO
      DO i = 1, nr
        SolV(sysid+(i-1)*SizeSys) = SolV(sysid+(i-1)*SizeSys)+REAL(Vec1(tid,i))
      END DO
    END DO
    sysid = sysid+sizeBlk
  END DO
  END SUBROUTINE


  ATTRIBUTES(GLOBAL) SUBROUTINE cuCRAM_GS_LOW(Diag, OffDiag, Vec, SolV, trVec, Vec1, nr, nnz, nsys)
  USE CUDA_CONST, ONLY : rowptr, colIdx, CRAM_Order, Pole, Res, Res0
  IMPLICIT NONE
!  COMPLEX(8), DIMENSION(:,:), DEVICE :: Diag, Vec1, Vec2
  COMPLEX(8), DIMENSION(:,:), DEVICE :: Vec1!, Vec2
  REAL(8), DIMENSION(:), DEVICE :: SolV
  REAL(8), DIMENSION(:,:), DEVICE :: Diag, OffDiag, Vec, trVec
  INTEGER, VALUE :: nr, nnz, Nsys

  INTEGER :: sizeBlk
  INTEGER :: tid, sysid
  INTEGER :: iord, iter, i, j, k, lb, ub
  INTEGER :: HalfOrder

  COMPLEX(8) :: polei, resi, temp1, val0
  REAL(8) :: norm, relerr, temp

  INTEGER :: SizeSys


  HalfOrder = CRAM_Order/2; SizeSys = Nsys

  sizeBlk = blockDim%x
  tid = (blockIdx%x-1)*sizeBlk+threadIdx%x
  sizeBlk = sizeBlk*GridDim%x
  sysid = tid
  DO WHILE(.TRUE.)
    IF (sysid .GT. nsys) EXIT
    DO i = 1, nr
      relerr = Vec(i,sysid) ! temporal for real variable
      trVec(tid,i) = relerr
      SolV(sysid+(i-1)*SizeSys) = relerr*Res0
    END DO
    DO iord = 1, HalfOrder
      polei = Pole(iord); resi = Res(iord);
      DO i = 1, nr
        Vec1(tid,i) = trVec(tid,i)
      END DO
      DO iter = 1, 50
!        DO i = 1, nr
!          Vec2(tid, i) = Vec1(tid, i)
!        END DO
        norm = 0.; relerr = 0.
        DO i = 1, nr
          val0 = Vec1(tid,i)
          lb = rowptr(i); ub = rowptr(i+1)-1
          temp1 = 2.*resi*trVec(tid,i)
          DO j = lb, ub
            k = colIdx(j)
            temp1 = temp1-OffDiag(sysid,j)*Vec1(tid,k)
          END DO
          temp1 = temp1/(Diag(sysid,i)-polei)
          Vec1(tid,i) = temp1

          temp = ABS(temp1- val0)
          relerr = relerr + temp*temp
          temp = ABS(temp1)
          norm = norm + temp*temp
        END DO
        relerr = relerr/norm
        IF (relerr<1.e-16) THEN
          EXIT
        END IF
      END DO
      DO i = 1, nr
        SolV(sysid+(i-1)*SizeSys) = SolV(sysid+(i-1)*SizeSys)+REAL(Vec1(tid,i))
      END DO
    END DO
    sysid = sysid+sizeBlk
  END DO
  END SUBROUTINE

  ATTRIBUTES(GLOBAL) SUBROUTINE cuCRAM_GS_ILP(Diag, OffDiag, Vec, SolV, trVec, Vec1, Vec2, nr, nnz, nsys)
  USE CUDA_CONST, ONLY : rowptr, colIdx, CRAM_Order, Pole, Res, Res0
  IMPLICIT NONE
!  COMPLEX(8), DIMENSION(:,:), DEVICE :: Diag, Vec1, Vec2
  COMPLEX(8), DIMENSION(:,:), DEVICE :: Vec1, Vec2
  REAL(8), DIMENSION(:), DEVICE :: SolV
  REAL(8), DIMENSION(:,:), DEVICE :: Diag, OffDiag, Vec, trVec
  INTEGER, VALUE :: nr, nnz, Nsys

  INTEGER :: sizeBlk
  INTEGER :: tid, sysid
  INTEGER :: iord, iter, i, j, k, lb, ub
  INTEGER :: HalfOrder

  COMPLEX(8) :: temp1(CRAM_ord), val0(CRAM_ord), temp2
  REAL(8) :: norm(CRAM_ord), relerr(CRAM_ord), temp

  INTEGER :: SizeSys


  HalfOrder = CRAM_Order/2; SizeSys = Nsys

  sizeBlk = blockDim%x
  tid = (blockIdx%x-1)*sizeBlk+threadIdx%x
  sizeBlk = sizeBlk*GridDim%x
  sysid = tid
  DO WHILE(.TRUE.)
    IF (sysid .GT. nsys) EXIT
    DO i = 1, nr
      temp = Vec(i,sysid) ! temporal for real variable
      trVec(tid,i) = temp
      SolV(sysid+(i-1)*SizeSys) = temp*Res0
    END DO

      DO i = 1, nr
        temp2 = trVec(tid,i)
!        Vec1(tid,(i-1)*CRAM_ord+1) = temp2
!        Vec1(tid,(i-1)*CRAM_ord+2) = temp2
!        Vec1(tid,(i-1)*CRAM_ord+3) = temp2
!        Vec1(tid,(i-1)*CRAM_ord+4) = temp2
!        Vec1(tid,(i-1)*CRAM_ord+5) = temp2
!        Vec1(tid,(i-1)*CRAM_ord+6) = temp2
!        Vec1(tid,(i-1)*CRAM_ord+7) = temp2
!#ifdef CRAM_16
!        Vec1(tid,(i-1)*CRAM_ord+8) = temp2
!#endif
        Vec1(tid,(i-1)*CRAM_ord+1:i*CRAM_ord) = temp2
      END DO
      DO iter = 1, 50
        DO i = 1, nr
!          val0(1) = Vec1(tid,(i-1)*CRAM_ord+1)
!          val0(2) = Vec1(tid,(i-1)*CRAM_ord+2)
!          val0(3) = Vec1(tid,(i-1)*CRAM_ord+3)
!          val0(4) = Vec1(tid,(i-1)*CRAM_ord+4)
!          val0(5) = Vec1(tid,(i-1)*CRAM_ord+5)
!          val0(6) = Vec1(tid,(i-1)*CRAM_ord+6)
!          val0(7) = Vec1(tid,(i-1)*CRAM_ord+7)
!#ifdef CRAM_16
!          val0(8) = Vec1(tid,(i-1)*CRAM_ord+8)
!#endif
          val0 = Vec1(tid,(i-1)*CRAM_ord+1:i*CRAM_ord)

          lb = rowptr(i); ub = rowptr(i+1)-1

          temp  = trVec(tid,i)
!          temp1(1) = 2.*Res(1)*temp
!          temp1(2) = 2.*Res(2)*temp
!          temp1(3) = 2.*Res(3)*temp
!          temp1(4) = 2.*Res(4)*temp
!          temp1(5) = 2.*Res(5)*temp
!          temp1(6) = 2.*Res(6)*temp
!          temp1(7) = 2.*Res(7)*temp
!#ifdef CRAM_16
!          temp1(8) = 2.*Res(8)*temp
!#endif
          temp1 = 2.*temp*Res
          DO j = lb, ub
            k = colIdx(j)
            temp = OffDiag(sysid,j)
!            temp1(1) = temp1(1)-temp*Vec1(tid,(k-1)*CRAM_ord+1)
!            temp1(2) = temp1(2)-temp*Vec1(tid,(k-1)*CRAM_ord+2)
!            temp1(3) = temp1(3)-temp*Vec1(tid,(k-1)*CRAM_ord+3)
!            temp1(4) = temp1(4)-temp*Vec1(tid,(k-1)*CRAM_ord+4)
!            temp1(5) = temp1(5)-temp*Vec1(tid,(k-1)*CRAM_ord+5)
!            temp1(6) = temp1(6)-temp*Vec1(tid,(k-1)*CRAM_ord+6)
!            temp1(7) = temp1(7)-temp*Vec1(tid,(k-1)*CRAM_ord+7)
!#ifdef CRAM_16
!            temp1(8) = temp1(8)-temp*Vec1(tid,(k-1)*CRAM_ord+8)
!#endif
            temp1 = temp1-temp*Vec1(tid,(k-1)*CRAM_ord+1:k*CRAM_ord)
          END DO
          temp2 = Diag(sysid,i)
!          temp1(1) = temp1(1)/(temp2-Pole(1))
!          temp1(2) = temp1(2)/(temp2-Pole(2))
!          temp1(3) = temp1(3)/(temp2-Pole(3))
!          temp1(4) = temp1(4)/(temp2-Pole(4))
!          temp1(5) = temp1(5)/(temp2-Pole(5))
!          temp1(6) = temp1(6)/(temp2-Pole(6))
!          temp1(7) = temp1(7)/(temp2-Pole(7))
!#ifdef CRAM_16
!          temp1(8) = temp1(8)/(temp2-Pole(8))
!#endif
          temp1 = temp1/(temp2-Pole)

          temp = ABS(temp1(1)- val0(1))
          relerr(1) = relerr(1) + temp*temp
          temp = ABS(temp1(1))
          norm(1) = norm(1)+ temp*temp

          temp = ABS(temp1(2)- val0(2))
          relerr(2) = relerr(2) + temp*temp
          temp = ABS(temp1(2))
          norm(2) = norm(2)+ temp*temp

          temp = ABS(temp1(3)- val0(3))
          relerr(3) = relerr(3) + temp*temp
          temp = ABS(temp1(3))
          norm(3) = norm(3)+ temp*temp

          temp = ABS(temp1(4)- val0(4))
          relerr(4) = relerr(4) + temp*temp
          temp = ABS(temp1(4))
          norm(4) = norm(4)+ temp*temp

          temp = ABS(temp1(5)- val0(5))
          relerr(5) = relerr(5) + temp*temp
          temp = ABS(temp1(5))
          norm(5) = norm(5)+ temp*temp

          temp = ABS(temp1(6)- val0(6))
          relerr(6) = relerr(6) + temp*temp
          temp = ABS(temp1(6))
          norm(6) = norm(6)+ temp*temp


          temp = ABS(temp1(7)- val0(7))
          relerr(7) = relerr(7) + temp*temp
          temp = ABS(temp1(7))
          norm(7) = norm(7)+ temp*temp

#ifdef CRAM_16
          temp = ABS(temp1(8)- val0(8))
          relerr(8) = relerr(8) + temp*temp
          temp = ABS(temp1(8))
          norm(8) = norm(8)+ temp*temp
#endif
!          Vec1(tid,(i-1)*CRAM_ord+1) = temp1(1)
!          Vec1(tid,(i-1)*CRAM_ord+2) = temp1(2)
!          Vec1(tid,(i-1)*CRAM_ord+3) = temp1(3)
!          Vec1(tid,(i-1)*CRAM_ord+4) = temp1(4)
!          Vec1(tid,(i-1)*CRAM_ord+5) = temp1(5)
!          Vec1(tid,(i-1)*CRAM_ord+6) = temp1(6)
!          Vec1(tid,(i-1)*CRAM_ord+7) = temp1(7)
!#ifdef CRAM_16
!          Vec1(tid,(i-1)*CRAM_ord+8) = temp1(8)
!#endif
          Vec1(tid,(i-1)*CRAM_ord+1:i*CRAM_ord) = temp1
        END DO
        temp = 0.0;
!        relerr(1) = relerr(1)/norm(1)
!        relerr(2) = relerr(2)/norm(2)
!        relerr(3) = relerr(3)/norm(3)
!        relerr(4) = relerr(4)/norm(4)
!        relerr(5) = relerr(5)/norm(5)
!        relerr(6) = relerr(6)/norm(6)
!        relerr(7) = relerr(7)/norm(7)
        relerr = relerr/norm
        temp = max(temp,relerr(1))
        temp = max(temp,relerr(2))
        temp = max(temp,relerr(3))
        temp = max(temp,relerr(4))
        temp = max(temp,relerr(5))
        temp = max(temp,relerr(6))
        temp = max(temp,relerr(7))
#ifdef CRAM_16
!        relerr(8) = relerr(8)/norm(8)
        temp = max(temp,relerr(8))
#endif
        IF (temp<1.e-16) THEN
          EXIT
        END IF
      END DO
      DO i = 1, nr
        temp1 = Vec1(tid, (i-1)*CRAM_ord+1:i*CRAM_ord)
        temp = REAL(temp1(1))
        temp = temp + REAL(temp1(2))
        temp = temp + REAL(temp1(3))
        temp = temp + REAL(temp1(4))
        temp = temp + REAL(temp1(5))
        temp = temp + REAL(temp1(6))
        temp = temp + REAL(temp1(7))
#ifdef CRAM_16
        temp = temp + REAL(temp1(8))
#endif
        SolV(sysid+(i-1)*SizeSys) = SolV(sysid+(i-1)*SizeSys)+temp
      END DO

    sysid = sysid+sizeBlk
  END DO
  END SUBROUTINE
END MODULE
#endif
