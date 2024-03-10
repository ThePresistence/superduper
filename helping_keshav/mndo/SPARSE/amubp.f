      subroutine amubp (nrow,a,ja,ia,b,jb,ib,c,jc,ic,nzmax,iw,ierr)
      IMPLICIT NONE
      integer nrow,nzmax,iw(nrow),j,ii,ka,jj,jcol,jpos,k,len,ierr,kb
      DOUBLE PRECISION scal
      DOUBLE PRECISION, DIMENSION (:), POINTER :: a,c,b
      INTEGER, DIMENSION (:), POINTER :: ja,ia,jc,ic,jb,ib
c-----------------------------------------------------------------------
c performs the matrix by matrix product C = A B 
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow  = integer. The row dimension of A = row dimension of C
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c 
c b, 
c jb, 
c ib    =  Matrix B in compressed sparse row format.
c
c nzmax = integer. The  length of the arrays c and jc.
c         amubp will stop if the result matrix C  has a number 
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
c         ierr .gt. 0 means that amubp stopped while computing the
c         i-th row  of C with i=ierr, because the number 
c         of elements in C exceeds nzmax.
c
c work arrays:
c------------
c iw    = integer work array of length equal to the number of
c         columns in A.
c Note: 
c-------
c   The row dimension of B is not needed. However there is no checking 
c   on the condition that ncol(A) = nrow(B). 
c
c----------------------------------------------------------------------- 
c
c     Comment by WT.
c     The check on overflow in the innermost loop is unnessary 
c     since we always perform a symbolic multiplication beforehand.
c     Removing this check reduces the cpu time by about 20 %.
c
      len = 0
      ic(1) = 1 
      ierr = 0
c     Initialize array iw.
      do 1 j=1, nrow
         iw(j) = 0
 1    continue
c
c     Outer loop over rows of matrix a.
      do 500 ii=1, nrow 
c     Inner loop over elements of row ii.
         do 200 ka=ia(ii), ia(ii+1)-1 
            scal = a(ka)
            jj   = ja(ka)
c           Innermost loop for given element a(ii,jj) of matrix a.
c           Loop over all nonzero elements of matrix b in row jj.
            do 100 kb=ib(jj),ib(jj+1)-1
               jcol = jb(kb)
               jpos = iw(jcol)
c              Column index jcol refers to matrix b and c. 
c              Check whether element c(ii,jcol) exists already.
c              If not (jpos=0) generate a new entry (len+1).
               if (jpos .eq. 0) then
                  len = len+1
c                 Check on overflow removed.
c                 if (len .gt. nzmax) then
c                    ierr = ii
c                    return
c                 endif
                  jc(len) = jcol
                  iw(jcol)= len
                  c(len)  = scal*b(kb)
               else
                  c(jpos) = c(jpos) + scal*b(kb)
               endif
 100        continue
 200     continue
c        End of loop over elements of row ii.
c        jc(k) are the columns of b(..) used when processing row ii.
c        Only the the corresponding elements of iw(..) are reset to 0.
         do 201 k=ic(ii), len
            iw(jc(k)) = 0
 201     continue
         ic(ii+1) = len+1
 500  continue
      return
      end
