      subroutine csrcscp (n,a,ja,ia,ao,jao,iao)
      IMPLICIT NONE
      integer i,k,n,j,next
      integer, DIMENSION (:), POINTER :: jao,iao,ja,ia
      DOUBLE PRECISION, DIMENSION (:), POINTER :: ao,a
c-----------------------------------------------------------------------
c Compressed Sparse Row     to      Compressed Sparse Column
c
c (transposition operation)   Not in place. 
c----------------------------------------------------------------------- 
c -- not in place --
c this subroutine transposes a matrix stored in a, ja, ia format.
c ---------------
c on entry:
c----------
c n     = dimension of A.
c
c a     = real array of length nnz (nnz=number of nonzero elements in input 
c         matrix) containing the nonzero elements.
c ja    = integer array of length nnz containing the column positions
c         of the corresponding elements in a.
c ia    = integer of size n+1. ia(k) contains the position in a, ja of
c         the beginning of the k-th row.
c
c on return:
c ---------- 
c output arguments:
c ao    = real array of size nzz containing the "a" part of the transpose
c jao   = integer array of size nnz containing the column indices.
c iao   = integer array of size n+1 containing the "ia" index array of
c         the transpose. 
c
c----------------------------------------------------------------------- 
c----------------- compute lengths of rows of transp(A) ----------------
      do 1 i=1,n+1
         iao(i) = 0
 1    continue
      do 3 i=1, n
         do 2 k=ia(i), ia(i+1)-1 
            j = ja(k)+1
            iao(j) = iao(j)+1
 2       continue 
 3    continue
c---------- compute pointers from lengths ------------------------------
      iao(1) = 1
      do 4 i=1,n
         iao(i+1) = iao(i) + iao(i+1)
 4    continue
c--------------- now do the actual copying ----------------------------- 
      do 6 i=1,n
         do 62 k=ia(i),ia(i+1)-1 
            j = ja(k) 
            next = iao(j)
            ao(next) = a(k)
            jao(next) = i
            iao(j) = next+1
 62      continue
 6    continue
c-------------------------- reshift iao and leave ---------------------- 
      do 7 i=n,1,-1
         iao(i+1) = iao(i)
 7    continue
      iao(1) = 1
      return
      end
