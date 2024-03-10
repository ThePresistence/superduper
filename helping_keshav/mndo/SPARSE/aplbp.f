      subroutine aplbp (nrow,a,ja,ia,b,jb,ib,c,jc,ic,nzmax,iw,ierr)
      IMPLICIT NONE
      integer nrow,nzmax,iw(nrow),j,ii,ka,k,jpos,len,kb,jcol,ierr
      DOUBLE PRECISION, DIMENSION (:), POINTER :: b,c,a
      INTEGER, DIMENSION (:), POINTER :: jb,ib,ic,jc,ja,ia
c-----------------------------------------------------------------------
c performs the matrix sum  C = A+B. 
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow  = integer. The row dimension of A and B
c
c a,
c ja,
c ia    = matrix A in compressed sparse row format.
c 
c b, 
c jb, 
c ib    = matrix B in compressed sparse row format.
c
c nzmax = integer. The length of the arrays c and jc.
c         aplbp will stop if the result matrix C  has a number 
c         of elements that exceeds exceeds nzmax. See ierr.
c 
c on return:
c----------
c c, 
c jc, 
c ic    = resulting matrix C in compressed sparse row sparse format.
c
c ierr  = integer. serving as error message. 
c         ierr = 0 means normal return,
c         ierr .gt. 0 means that aplbp stopped while computing the
c         i-th row  of C with i=ierr, because the number 
c         of elements in C exceeds nzmax.
c
c work arrays:
c------------
c iw    = integer work array of length equal to the number of
c         columns in A.
c
c-----------------------------------------------------------------------
c
c     Comment by WT.
c     The check on overflow in the innermost loop is unnessary 
c     since we always perform a symbolic addition beforehand.
c
      ierr = 0
      len = 0
      ic(1) = 1 
      do 1 j=1,nrow
         iw(j) = 0
 1    continue
c     
      do 500 ii=1, nrow
c     row i 
         do 200 ka=ia(ii), ia(ii+1)-1 
            len = len+1
            jcol    = ja(ka)
c           check on overflow removed.
c           if (len .gt. nzmax) then
c              ierr = ii
c              return
c           endif
            jc(len) = jcol 
            c(len)  = a(ka) 
            iw(jcol)= len
 200     continue
c     
         do 300 kb=ib(ii),ib(ii+1)-1
            jcol = jb(kb)
            jpos = iw(jcol)
            if (jpos .eq. 0) then
               len = len+1
c              check on overflow removed.
c              if (len .gt. nzmax) then
c                 ierr = ii
c                 return
c              endif
               jc(len) = jcol
               c(len)  = b(kb)
               iw(jcol)= len
            else
               c(jpos) = c(jpos) + b(kb)
            endif
 300     continue
         do 301 k=ic(ii), len
            iw(jc(k)) = 0
 301     continue
         ic(ii+1) = len+1
 500  continue
      return
      end
