
subroutine matadd(A,B,C,LDA,NROW,NCOL)
  implicit none
  integer i,LDA,NROW,NCOL
  real*8  A(LDA,*),B(LDA,*),C(NROW,*)
  do i=1,NCOL
    call vecadd(A(1,i),B(1,i),C(1,i),NROW)
  enddo
end subroutine matadd

subroutine matadd2(A,scale,C,LDA,NROW,NCOL)
  implicit none
  integer i,LDA,NROW,NCOL
  real*8  scale
  real*8  A(LDA,*),C(NROW,*)
  do i=1,NCOL
    call vecadd3(A(1,i),scale,C(1,i),NROW)
  enddo
end subroutine matadd2

subroutine matsub(A,B,C,LDA,NROW,NCOL)
  implicit none
  integer i,LDA,NROW,NCOL
  real*8  A(LDA,*),B(LDA,*),C(NROW,*)
  do i=1,NCOL
    call vecsub(A(1,i),B(1,i),C(1,i),NROW)
  enddo
end subroutine matsub

subroutine MatDiagSquare(A,EigenValue,N,Info)
implicit NONE
real*8 A(*),EigenValue(*)
INTEGER N,J,LiWork,LenWork,Info,LenWork2
CHARACTER JobZ, UpLo
real*8, dimension(:), allocatable :: A2,Work,Work2
INTEGER, dimension(:), allocatable :: IWork

allocate(A2(N*N))
call VecCopy(A2,A,N*N)

Info = 0
JobZ = 'V'
UpLo = 'U'

LenWork = 10
LiWork = -1
allocate(Work(10))
allocate(IWork(10))
call dsyevd(JobZ, UpLo, N, A, N, EigenValue, Work, LenWork, &
           IWork, LiWork,Info)
LenWork = INT( Work(1))
LiWork = IWork(1)
deallocate(Work)
deallocate(IWork)
allocate(Work(LenWork))
allocate(IWork(LiWork))
call dsyevd(JobZ, UpLo, N, A, N, EigenValue, Work, LenWork, &
            IWork, LiWork,Info)
deallocate(Work)
deallocate(IWork)

if (Info .ne. 0) then
   write(6,*) " Info = %d\n",Info
   write(6,*) " Call to dsyevd failed in Diagonalize"
   LenWork2 = 32*N
   allocate(Work2(LenWork2))
   call VecCopy(A,A2,N*N)
   call dsyev(JobZ, UpLo, N, A, N, EigenValue, Work2, LenWork2, &
              Info)
   if (Info .ne. 0 ) then
     write(6,*) " Info = %d\n",Info
     write(6,*) " Call to dsyevd failed in Diagonalize"
     stop
   endif
   deallocate(Work2)
endif

deallocate(A2)
end

subroutine MatMult(A,B,C,M,N,K,LDA,LDB,LDC,METHOD)
implicit real*8 (a-h,o-z)
real*8 A(LDA,*),B(LDB,*),C(LDC,*)
CHARACTER TRANSA,TRANSB

TRANSA = 'N'
TRANSB = 'N'
IFLAG = ABS(METHOD)
IF (IFLAG.eq.2.or.IFLAG.eq.4) TRANSA = 'T'
IF (IFLAG.eq.3.or.IFLAG.eq.4) TRANSB = 'T'

Alpha = 1D0
Beta  = 0D0
call DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

return
end

subroutine MATPRNT(A,NROW,NCOL,N)
real*8 A(NROW,NCOL)
M = NCOL/N
IF(M.GT.0) THEN
   do I=1,M
      IMIN = (I-1)*N + 1
      IMAX = I*N
      WRITE(6,'(/)')
      do J=1,NROW
         WRITE(6,'(I3,8F12.6)') J,A(J,(I-1)*N+1:I*N)
      enddo
   enddo
endIF
L = NCOL-M*N
IF(L.EQ.0) return
WRITE(6,'(/)')
do J=1,NROW
   WRITE(6,'(I3,8F12.6)') J,A(J,M*N+1:NCOL)
enddo
return
end

! Matrix symmetrization and stored in upper triangular matrix
subroutine MatSym(A,N,LM4)
implicit real*8 (a-h,o-z)
real*8, Dimension(:), Allocatable :: B
real*8 A(N,*)

do I=1,N
  do J=I+1,N
    A(I,J) = 0.5 *(A(I,J) + A(J,I))
    A(J,I) = A(I,J)
  enddo
enddo

allocate(B(LM4))
IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    B(IJ) = A(J,I)
  enddo
enddo
call VecCopy(A(1,1),B,LM4)
deallocate(B)

return
end

! Matrix symmetrization
subroutine AddMatSym(A,N)
implicit real*8 (a-h,o-z)
real*8 A(N,*)

do I=1,N
  do J=I,N
    A(I,J) = (A(I,J) + A(J,I))
    A(J,I) = A(I,J)
  enddo
enddo

return
end

subroutine MatTrans(A,B,M,N)
implicit real*8 (a-h,o-z)
real*8 A(M,*),B(N,*)

do I=1,M
  do J=1,N
    A(I,J) = B(J,I)
  enddo
enddo

return
end

subroutine MATTRANS2(A,M,N)
implicit real*8 (a-h,o-z)
real*8 A(N,*)
real*8, dimension(:,:), allocatable :: B

allocate(B(M,N))
call VECCOPY(B,A,M*N)
do I=1,M
  do J=1,N
    A(J,I) = B(I,J)
  enddo
enddo
deallocate(B)

return

end
