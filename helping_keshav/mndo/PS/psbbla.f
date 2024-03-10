C     ******************************************************************
C
C     Patch for a (semi-)defective BLAS implementation on IBM/RS6k.
C
C     ******************************************************************
C
C     At least on some versions of AIX (3.2.5), BLAS levels 2 and 3
C     routines DGEMV, DGBMV, DSYMV, DSBMV, DSPMV, DGEMM, DSYMM,
C     DSYRK, and DSYR2K with zero parameter BETA will fail if any
C     of elements in output array/vector are initialized to NAN or
C     INF.
C   
C     Since this interpretation is allowed by the letter of the
C     BLAS definition, similar problems might be expected on other
C     implementation, especially if MADD and ADD can be issued
C     at the same cost.
C   
C     (Hopefully) All cases where the problem might be encountered
C     contain calls to PSZRVB/PSZRMB routines, which should initialize
C     the area to an arbitrary non-NAN, non-INF value.
C
      SUBROUTINE PSZRVB(LEN,X)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(LEN)
CZR   CALL PSZRV(LEN,X)
      RETURN
      END
C
      SUBROUTINE PSZRMB(L1,L2,X,LDX)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(LDX,L2)
CZR   CALL PSZRM(L1,L2,X,LDX)
      RETURN
      END
