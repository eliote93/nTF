module inflowMathUtil
    
    contains
    
    subroutine inverseA(A,n)
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: n
        REAL(8),INTENT(INOUT) :: A(n,n)
        REAL(8) :: Ainv(n,n),pivot,temp(n),rpivot,fmult,UdiagInv(n),y(n),b(n),sum
        INTEGER :: P(n),Pmat(n,n),k,ipivot,i,itemp,j
        !! LU Factorization with Pivoting
        ! Pivot Vector
        do i=1,n; P(i)=i; enddo
        ! Column Elimination
        do k=1,n-1
            ! determine pivot for the k-th column elimination
            ipivot=k    
            pivot=abs(A(k,k))
            do i=k+1,n
                if (pivot.lt.abs(A(i,k))) then
                    ipivot=i
                    pivot=abs(A(i,k))
                endif
            enddo
            ! swap
            itemp=p(k); p(k)=p(ipivot); p(ipivot)=itemp;
            temp=A(k,:); A(k,:)=A(ipivot,:); A(ipivot,:)=temp;
            ! multiplier
            rpivot=1._8/A(k,k)
            do i=k+1,n
                fmult=rpivot*A(i,k)
                A(i,k)=fmult
                do j=k+1,n
                    A(i,j)=A(i,j)-fmult*A(k,j)
                enddo
            enddo
        enddo    
        Ainv=0._8
        ! Permutation Matrix
        Pmat=0
        do i=1,n; Pmat(i,P(i))=1; enddo
        !! Find Inverse
        do i=1,n; UdiagInv(i)=1._8/A(i,i); enddo
        ! column by column
        do i=1,n
            ! forward substitution
            y=0._8    
            b=Pmat(:,i)
            y(1)=b(1)
            do k=2,n
                sum=0._8
                do j=1,k-1
                    sum=sum+A(k,j)*y(j)
                enddo
                y(k)=b(k)-sum
            enddo
            ! backward substitution
            Ainv(n,i)=y(n)*UdiagInv(n)
            do k=n-1,1,-1
                sum=0._8
                do j=n,k+1,-1
                    sum=sum+A(k,j)*Ainv(j,i)
                enddo
                Ainv(k,i)=(y(k)-sum)*UdiagInv(k)
            enddo
        enddo
        A=Ainv
    end subroutine
        
    subroutine matvecprod(A,b,c,n)
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: n
        REAL(8),INTENT(IN) :: A(n,n),b(n)
        REAL(8),INTENT(OUT) :: c(n)
        INTEGER :: i,j
        c=0._8
        do i=1,n
            do j=1,n
                c(i)=c(i)+A(i,j)*b(j)
            enddo
        enddo
    end subroutine
    
    function rmserr(y,n1s,n1,n2s,n2)
        IMPLICIT NONE
        INTEGER :: n1s,n1,n2s,n2,i1,i2
        REAL(8) :: rmserr
        REAL(8) :: y(n1s:n1,n2s:n2)
        rmserr=0._8
        do i1=n1s,n1
            do i2=n2s,n2
                rmserr=rmserr+y(i1,i2)*y(i1,i2)
            enddo
        enddo
        rmserr=rmserr/REAL(n1-n1s+1,8)/REAL(n2-n2s+1,8)
        rmserr=sqrt(rmserr)
    end function
 
    function rmserr1D(y,n1s,n1)
        IMPLICIT NONE
        INTEGER :: n1s,n1,i1
        REAL(8) :: rmserr1D
        REAL(8) :: y(n1s:n1)
        rmserr1D=0._8
        do i1=n1s,n1
            rmserr1D=rmserr1D+y(i1)*y(i1)
        enddo
        rmserr1D=rmserr1D/REAL(n1-n1s+1,8)
        rmserr1D=sqrt(rmserr1D)
    end function    
    
end module