      SUBROUTINE DYNDIA (IDIAG,N,LWORK,LIWORK)
C     *
C     DEFINE LENGTH OF SCRATCH ARRAYS FOR DIAGONALIZATION.
C     *
C     NOTATION. I=INPUT, O=OUTPUT.
C     IDIAG     DIAGONALIZATION ROUTINE TO BE USED (I).
C               = 0  EVVRSP (EISPACK-BASED), SQUARE MATRIX.
C               = 1  TDIAG  (EISPACK-BASED), PACKED MATRIX.
C               = 2  TQL2/TRED2 (EISPACK)  , SQUARE MATRIX.
C               = 3  DSPEV  (LAPACK), PACKED MATRIX.
C               = 4  DSPEVX (LAPACK), PACKED MATRIX.
C               = 5  DSYEV  (LAPACK), SQUARE MATRIX.
C               = 6  DSYEVX (LAPACK), SQUARE MATRIX.
C               = 7  DSPEVD (LAPACK), PACKED MATRIX, DIVIDE-CONQUER.
C               = 8  DSYEVD (LAPACK), SQUARE MATRIX, DIVIDE-CONQUER.
C     N         ORDER OF THE MATRIX (I).
C     LWORK     LENGTH OF REAL SCRATCH ARRAYS (O).
C     LIWORK    LENGTH OF INTEGER SCRATCH ARRAYS (O).
C     *
C     THE DEFINITIONS FOR THE LAPACK ROUTINES (IDIAG=3-8) ARE BASED ON
C     THE SPECIFICATIONS IN THE LAPACK USERS' GUIDE (2ND EDITION, 1995).
C     THIS GUIDE CONTAINS SPECIAL RECOMMENDATIONS FOR DSYEV AND DSYEVX
C     WHEN USING BLOCK ALGORITHMS TO IMPROVE THE EFFICIENCY: IN THIS
C     CASE, THE LENGTH OF THE REAL SCRATCH ARRAYS SHOULD BE AT LEAST
C     LWORK = (NB+2)*N  FOR DSYEV  (IDIAG=5) AND
C     LWORK = (NB+3)*N  FOR DSYEVX (IDIAG=6)
C     WHERE NB IS THE OPTIMUM BLOCK SIZE FROM LAPACK ROUTINE ILAENV.
C     INSPECTION OF THE STANDARD LAPACK(2.0) VERSION OF ILAENV SHOWS,
C     HOWEVER, THAT NB=1 BY DEFAULT FOR DSYEV AND DSYEVX (UNBLOCKED).
C     HENCE, THE STANDARD LWORK VALUES FROM THE GUIDE ARE ADOPTED.
C     IF FUTURE RELEASES EMPLOY BLOCKING FOR DSYEV AND DSYEVX, THE
C     DEFINITION OF LWORK MUST BE ADJUSTED.
C     *
C     IN THE COMPUTATIONAL ROUTINES, THERE WILL BE NO CHECKS ON THE
C     LENGTH OF THE SCRATCH ARRAYS FOR THE LAPACK ROUTINES SINCE THEY
C     ARE DEFINED HERE IN ACCORDANCE WITH THE LAPACK GUIDE.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C *** BRANCH ACCORDING TO IDIAG VALUE.
      IF(IDIAG.LE.0 .OR. IDIAG.GT.8) THEN
         LWORK  = 8*N
         LIWORK = N
      ELSE IF(IDIAG.EQ.1) THEN
         LWORK  = 2*N
         LIWORK = 0
      ELSE IF(IDIAG.EQ.2) THEN
         LWORK  = N
         LIWORK = 0
      ELSE IF(IDIAG.EQ.3) THEN
         LWORK  = 3*N
         LIWORK = 0
      ELSE IF(IDIAG.EQ.4) THEN
         LWORK  = 8*N
         LIWORK = 6*N
      ELSE IF(IDIAG.EQ.5) THEN
         NB = ILAENV(1, 'DSYEV', 'VU', N, -1, -1, -1)
         IF(NB.GT.1) THEN
            LWORK = (NB+2)*N
         ELSE
            LWORK = 3*N
         ENDIF
         LIWORK = 0
      ELSE IF(IDIAG.EQ.6) THEN
         NB = ILAENV(1, 'DSYEVX', 'VIU', N, -1, -1, -1)
         IF(NB.GT.5) THEN
            LWORK = (NB+3)*N
         ELSE
            LWORK = 8*N
         ENDIF
         LIWORK = 6*N
      ELSE IF(IDIAG.EQ.7) THEN
         LWORK  = N*N+6*N+1
         LIWORK = 5*N+3
      ELSE IF(IDIAG.EQ.8) THEN
         LWORK  = 2*N*N+6*N+1
         LIWORK = 5*N+3
      ENDIF
      RETURN
      END
