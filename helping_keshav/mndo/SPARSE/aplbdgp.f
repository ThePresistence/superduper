      subroutine aplbdgp (nrow,ja,ia,jb,ib,nnz,iw)
      IMPLICIT NONE
      integer nrow,iw(nrow),ldg,last,j,jr,nnz,ii,k,jc,ie
      INTEGER, DIMENSION (:), POINTER :: jb,ib,ja,ia
c-----------------------------------------------------------------------
c gets the number of nonzero elements in each row of A+B and the total 
c number of nonzero elements in A+B. 
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow = integer. Row dimension of A and B
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format
c 
c b, 
c jb, 
c ib   = Matrix B in compressed sparse row format
c
c on return:
c----------
c nnz  = total number of nonzero elements found in A * B
c
c work arrays:
c------------
c iw   = integer work array of length equal to nrow
c
c-----------------------------------------------------------------------
      do 2 k=1, nrow
         iw(k) = 0
 2    continue
c
      nnz = 0
c
      do 7 ii=1,nrow 
         ldg = 0 
c     
c    end-of-linked list
c     
         last = -1 
c     
c     row of A
c     
         do 5 j = ia(ii),ia(ii+1)-1 
            jr = ja(j) 
c     
c     add element to the linked list 
c     
            ldg = ldg + 1
            iw(jr) = last 
            last = jr
 5       continue
c     
c     row of B
c     
         do 6 j=ib(ii),ib(ii+1)-1
            jc = jb(j)
            if (iw(jc) .eq. 0) then 
c     
c     add one element to the linked list 
c     
               ldg = ldg + 1
               iw(jc) = last 
               last = jc
            endif
 6       continue
c     done with row ii. 
c
         nnz = nnz+ldg
c     
c     reset iw to zero
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
