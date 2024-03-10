      subroutine amubdgp (n,ja,ia,jb,ib,nnz,iw)
      IMPLICIT NONE
      INTEGER n,nnz,last,j,ii,k,ldg,jr,jc,iw(n),ie
      INTEGER, DIMENSION (:), POINTER :: ja,ia,jb,ib
c-----------------------------------------------------------------------
c gets the number of nonzero elements in each row of A*B and the total 
c number of nonzero elements in A*B. 
c-----------------------------------------------------------------------
c on entry:
c -------- 
c
c n  = integer. row dimension of matrix A
c n  = integer. column dimension of matrix A = row dimension of matrix B
c n  = integer. colum dimension of matrix B
c
c ja, ia = row structure of input matrix A
c ja = column indices of the nonzero elements of A stored by rows
c ia = pointer to beginning of each row in ja
c
c jb, ib = row structure of input matrix B
c jb = column indices of the nonzero elements of B stored by rows
c ib = pointer to beginning of each row in jb
c 
c on return:
c ---------
c nnz = total number of nonzero elements found in A * B
c
c work arrays:
c-------------
c iw  = integer work array of length n
c-----------------------------------------------------------------------
      do 1 k=1,n
         iw(k) = 0 
 1    continue
c      
c     method used: Transp(A) * A = sum [over i=1, n]  a(i)^T a(i)
c     where a(i) = i-th row of  A. We must be careful not to add  the
c     elements already accounted for.
c     
c
      nnz = 0
c     
      do 7 ii=1,n 
c     
c     for each row of A
c     
         ldg = 0 
c     
c    end-of-linked list
c     
         last = -1 
         do 6 j = ia(ii),ia(ii+1)-1 
c     
c     row number to be added:
c
            jr = ja(j) 
            do 5 k=ib(jr),ib(jr+1)-1
               jc = jb(k) 
               if (iw(jc) .eq. 0) then 
c     
c     add one element to the linked list 
c     
                  ldg = ldg + 1
                  iw(jc) = last 
                  last = jc
               endif
 5          continue
 6       continue
c
         nnz = nnz+ldg
c     
c     reset iw to zero
c     only for elements active in row ii of A
c     
         do 61 k=1,ldg 
            j = iw(last) 
            iw(last) = 0
            last = j
 61      continue
c-----------------------------------------------------------------------
 7    continue
c     
      return
      end
