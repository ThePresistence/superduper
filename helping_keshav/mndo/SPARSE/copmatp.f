      subroutine copmatp (nrow,a,ja,ia,ao,jao,iao,mode)
      IMPLICIT NONE
      integer nrow,i,mode
      INTEGER, DIMENSION (:), POINTER :: ja,ia,iao,jao
      DOUBLE PRECISION, DIMENSION (:), POINTER :: a,ao 
c----------------------------------------------------------------------
c copies the matrix a, ja, ia, into the matrix ao, jao, iao. 
c----------------------------------------------------------------------
c on entry:
c---------
c nrow  = row dimension of the matrix 
c a,
c ja,
c ia    = input matrix in compressed sparse row format. 
c
c on return:
c----------
c ao,
c jao,
c iao   = output matrix containing the same data as a, ja, ia.
c-----------------------------------------------------------------------
c           Y. Saad, March 1990. 
c-----------------------------------------------------------------------
      do 100 i = 1, nrow+1
         iao(i) = ia(i)
 100  continue
c     
      if (mode.eq.1) then
         do 201 i=ia(1), ia(nrow+1)-1
            ao(i) = a(i)
            jao(i)= ja(i)
 201     continue
      else
         do 202 i=ia(1), ia(nrow+1)-1
            ao(i) = -a(i)
            jao(i)= ja(i)
 202     continue
      endif
c
      return
      end
