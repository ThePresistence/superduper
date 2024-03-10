      module module3
      INTERFACE
      subroutine aplbdgp (nrow,ja,ia,jb,ib,nnz,iw)
      IMPLICIT NONE
      integer nrow,iw(nrow),nnz
      INTEGER, DIMENSION (:), POINTER :: jb,ib,ja,ia
      end subroutine aplbdgp
!
      subroutine amubp (nrow,a,ja,ia,b,jb,ib,c,jc,ic,nzmax,iw,ierr)
      IMPLICIT NONE
      integer nrow,nzmax,iw(nrow),ierr
      DOUBLE PRECISION, DIMENSION (:), POINTER :: a,c,b
      INTEGER, DIMENSION (:), POINTER :: ja,ia,jc,ic,jb,ib
      end subroutine amubp
!
      subroutine amubdgp (n,ja,ia,jb,ib,nnz,iw)
      IMPLICIT NONE
      INTEGER n,nnz,iw(n)
      INTEGER, DIMENSION (:), POINTER :: ja,ia,jb,ib
      end subroutine amubdgp
!
      subroutine filterp (n,job,drptol,a,ja,ia,b,jb,ib,len,ierr)
      IMPLICIT NONE
      DOUBLE PRECISION drptol
      integer n,job,len,ierr
      DOUBLE PRECISION, DIMENSION (:), POINTER :: a,b
      integer, DIMENSION (:), POINTER :: ja,jb,ia,ib
      end subroutine filterp
!
      subroutine aplbp (nrow,a,ja,ia,b,jb,ib,c,jc,ic,nzmax,iw,ierr)
      IMPLICIT NONE
      integer nrow,nzmax,iw(nrow),ierr
      DOUBLE PRECISION, DIMENSION (:), POINTER :: b,c,a
      INTEGER, DIMENSION (:), POINTER :: jb,ib,ic,jc,ja,ia
      end subroutine aplbp
!
      subroutine transpp (nrow,ncol,a,ja,ia,iwk,ierr)
      integer nrow,ncol,iwk(*),ierr
      integer, DIMENSION (:), POINTER :: ja,ia
      DOUBLE PRECISION, DIMENSION (:), POINTER :: a
      end subroutine transpp
!
      SUBROUTINE TRACENP (NO,A,IAA,KKA,B,IAB,KKB,Y,Tr)
      IMPLICIT NONE
      INTEGER NO
      DOUBLE PRECISION Tr
      DOUBLE PRECISION, DIMENSION (:), POINTER :: A,B,Y
      INTEGER, DIMENSION (:), POINTER :: IAA,KKA,IAB,KKB
      END SUBROUTINE TRACENP
!
      SUBROUTINE TRACEAP (NO,A,IAA,KKA,Tr)
      IMPLICIT NONE
      DOUBLE PRECISION Tr
      INTEGER NO
      DOUBLE PRECISION, DIMENSION (:), POINTER :: A
      INTEGER, DIMENSION (:), POINTER :: IAA,KKA
      END SUBROUTINE TRACEAP
!
      SUBROUTINE SPAUPM (A,B,IIA,IIB,JJA,JJB,N,FACT1,FACT2,UPMIN1, & 
                         UPMAX1,UPMIN2,UPMAX2)
      IMPLICIT NONE
      INTEGER :: N
      DOUBLE PRECISION :: UPMIN1,UPMAX1,FACT1,FACT2,UPMIN2,UPMAX2
      INTEGER, DIMENSION (:), POINTER :: IIB,JJB,IIA,JJA
      DOUBLE PRECISION, DIMENSION (:), POINTER :: B,A
      END SUBROUTINE SPAUPM
!
      subroutine csrcscp (n,a,ja,ia,ao,jao,iao)
      IMPLICIT NONE
      integer n
      integer, DIMENSION (:), POINTER :: jao,iao,ja,ia
      DOUBLE PRECISION, DIMENSION (:), POINTER :: ao,a
      end subroutine csrcscp
!
      integer function idamaxp(n,dx) 
      DOUBLE PRECISION, DIMENSION (:), POINTER :: dx
      integer n
      end function idamaxp
!
      subroutine coicsrp (n,nnz,a,ja,ia,iwk)
      IMPLICIT NONE
      integer nnz,n,iwk(n+1)
      INTEGER, DIMENSION (:), POINTER :: ja,ia
      DOUBLE PRECISION, DIMENSION (:), POINTER :: a
      end subroutine coicsrp
!
      subroutine copmatp (nrow,a,ja,ia,ao,jao,iao,mode)
      IMPLICIT NONE
      integer nrow,mode
      INTEGER, DIMENSION (:), POINTER :: ja,ia,iao,jao
      DOUBLE PRECISION, DIMENSION (:), POINTER :: a,ao 
      end subroutine copmatp
!
      SUBROUTINE SPAPRT (A,IIA,JJA,N)
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION (:), POINTER ::A
      INTEGER, DIMENSION (:), POINTER :: IIA,JJA
      INTEGER N
      END SUBROUTINE SPAPRT
      END INTERFACE
      end module module3
