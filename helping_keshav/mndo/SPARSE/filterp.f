      subroutine filterp (n,job,drptol,a,ja,ia,b,jb,ib,len,ierr)
      IMPLICIT NONE
      DOUBLE PRECISION drptol
      integer n,job,len,ierr
      DOUBLE PRECISION, DIMENSION (:), POINTER :: a,b
      integer, DIMENSION (:), POINTER :: ja,jb,ia,ib
c-----------------------------------------------------------------------
c     This module removes any elements whose absolute value
c     is small from an input matrix A and puts the resulting
c     matrix in B.  The input parameter job selects a definition
c     of small.
c-----------------------------------------------------------------------
c on entry:
c---------
c  n     = integer. row dimension of matrix
c  job   = integer. used to determine strategy chosen by caller to
c         drop elements from matrix A. 
c          job = 1  
c              Elements whose absolute value is less than the
c              drop tolerance are removed.
c          job = 2
c              Elements whose absolute value is less than the 
c              product of the drop tolerance and the Euclidean
c              norm of the row are removed. 
c          job = 3
c              Elements whose absolute value is less that the
c              product of the drop tolerance and the largest
c              element in the row are removed.
c 
c drptol = real. drop tolerance used for dropping strategy.
c a
c ja
c ia     = input matrix in compressed sparse format
c len    = integer. the amount of space available in arrays b and jb.
c
c on return:
c---------- 
c b
c jb
c ib    = resulting matrix in compressed sparse format. 
c 
c ierr  = integer. containing error message.
c         ierr .eq. 0 indicates normal return
c         ierr .gt. 0 indicates that there is'nt enough
c         space is a and ja to store the resulting matrix.
c         ierr then contains the row number where filter stopped.
c note:
c       This module is in place. (b,jb,ib can be the same as 
c       a, ja, ia in which case the result will be overwritten).
c----------------------------------------------------------------------c
c           contributed by David Day,  Sep 19, 1989.                   c
c----------------------------------------------------------------------c
c
c     Comment by WT.
c     The check on overflow in the innermost loop is unnessary 
c     since we always use arrays a and b which are the same.
c     b cannot become larger than a upon filtering.
c
c local variables
      DOUBLE PRECISION norm,loctol
      integer index,row,k,k1,k2 
c
      index = 1
      do 10 row= 1,n
         k1 = ia(row)
         k2 = ia(row+1) - 1
         ib(row) = index
         goto (100,200,300) job
 100     norm = 1.0d0
         goto 400
 200     norm = 0.0d0
         do 22 k = k1,k2
            norm = norm + a(k) * a(k)
 22      continue
         norm = sqrt(norm)
         goto 400
 300     norm = 0.0d0
         do 23 k = k1,k2
            if( abs(a(k))  .gt. norm) then
               norm = abs(a(k))
            endif
 23      continue
 400     loctol = drptol * norm
         do 30 k = k1,k2
            if( abs(a(k)) .gt. loctol)then 
c              check on overflow removed.
c              if (index .gt. len) then
c                 ierr = row 
c                 return
c              endif
               b(index) =  a(k)
               jb(index) = ja(k)
               index = index + 1
            endif
 30      continue
 10   continue
      ib(n+1) = index
      return
      end
