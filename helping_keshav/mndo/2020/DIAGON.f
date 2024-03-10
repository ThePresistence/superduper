      SUBROUTINE DIAGON (A,C,E,Q,IWORK,LM2,LM4,N,NROOT,IDIAG,MODE,INFO)
C     *
C     CONTROL ROUTINE FOR DIAGONALIZATIONS.
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     A(*)      SYMMETRIC REAL MATRIX (I).
C               MODE=0: PACKED COLUMNWISE IN A LINEAR ARRAY.
C               MODE=1: FULL TWO-DIMENSIONAL MATRIX.
C     C(LM2,N)  EIGENVECTORS (O).
C     E(N)      EIGENVALUES IN ASCENDING ORDER (O).
C     Q(*)      WORKSPACE, REAL ARRAY (S).
C     IWORK     WORKSPACE, INTEGER ARRAY (S).
C     LM2       LEADING DIMENSION OF ARRAY C (I).
C     LM4       DIMENSION OF ARRAY A (I).
C     N         ORDER OF THE MATRICES A AND C (I).
C     NROOT     NUMBER OF EIGENVECTORS REQUESTED (1...NROOT) (I).
C     IDIAG     DIAGONALIZATION ROUTINE TO BE USED (I).
C               = 0  EWWRSP (EISPACK-BASED), SQUARE MATRIX.
C               = 1  TDIAG  (EISPACK-BASED), PACKED MATRIX.
C               = 2  TQL2/TRED2 (EISPACK)  , SQUARE MATRIX.
C               = 3  DSPEV  (LAPACK), PACKED MATRIX.
C               = 4  DSPEVX (LAPACK), PACKED MATRIX.
C               = 5  DSYEV  (LAPACK), SQUARE MATRIX.
C               = 6  DSYEVX (LAPACK), SQUARE MATRIX.
C               = 7  DSPEVD (LAPACK), PACKED MATRIX, DIVIDE-CONQUER.
C               = 8  DSYEVD (LAPACK), SQUARE MATRIX, DIVIDE-CONQUER.
C     MODE      MATRIX FORMAT(I).
C               = 0  LINEARLY PACKED SYMMETRIC MATRIX.
C               = 1  FULL TWO-DIMENSIONAL MATRIX.
C     INFO      ERROR FLAG (O).
C               = 0  NO ERROR.
C     *
C     THE MEMORY ALLOCATION IS HANDLED IN ROUTINES DYNSCF AND DYNDIA.
C     SEE COMMENTS IN THESE ROUTINES FOR DETAILS. GENERAL REMARKS:
C     - IF THE PACKED INPUT MATRIX A(LM4) IS EXPANDED TO A SQUARE
C       MATRIX (IDIAG=0,2,5,6,8), THE RESULTING MATRIX A(LM2,LM2)
C       MAY OCCUPY THE SAME SPACE (SAME ADDRESS FOR A(1) AND A(1,1)).
C       THE CORRESPONDING MEMORY LOCATIONS ARE RESERVED BY DYNSCF.
C     - THERE ARE NO CHECKS ON Q(1) AND/OR IWORK(1) AFTER RETURNING
C       FROM LAPACK ROUTINES (IDIAG=5,6,7,8) BECAUSE THE CORRESPONDING
C       DIAGNOSTICS SHOULD NOT APPEAR DUE TO THE FACT THAT THE MEMORY
C       HAS BEEN ALLOCATED PROPERLY BY DYNDIA.
C     - DSPEVD AND DSYEVD REQUIRE SIGNIFICANTLY MORE SCRATCH SPACE THAN
C       THE OTHER ROUTINES (SEE DYNDIA). CALCULATIONS THAT FAIL WITH
C       THESE ROUTINES DUE TO MEMORY LIMITATIONS MAY STILL BE FEASIBLE
C       WHEN USING OTHER DIAGONALIZATION ROUTINES. IT IS UP TO THE USER
C       TO CHOOSE THE APPROPRIATE DIAGONALIZER IN THIS CASE.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*1 RANGE
      CHARACTER*6 SUBROU
      COMMON /NBFILE/ NBF(20)
      DIMENSION A(LM4),C(LM2,N),E(N)
      DIMENSION Q(*),IWORK(LM2*6)
      DIMENSION SUBROU(9)
      DATA SUBROU /'EWWRSP','TDIAG ','TQL2  ','DSPEV ','DSPEVX',
     1             'DSYEV ','DSYEVX','DSPEVD','DSYEVD'/
C *** BRANCH ACCORDING TO IDIAG VALUE.
      INFO = 0
      IF(IDIAG.LE.0 .OR. IDIAG.GT.8) THEN
         IF(MODE.EQ.0) CALL SQUARE(A,A,N,LM2,LM4)
         CALL EWWRSP(N,NROOT,LM2,A,Q,IWORK,E,C,INFO)
      ELSE IF(IDIAG.EQ.1) THEN
         IF(MODE.NE.0) CALL LINEAR(A,A,N,LM2,LM4)
         CALL TDIAG (A,C,E,Q,LM4,LM2,N,NROOT)
      ELSE IF(IDIAG.EQ.2) THEN
         IF(MODE.EQ.0) CALL SQUARE(A,A,N,LM2,LM4)
         CALL TRED2 (LM2,N,A,E,Q,C)
         CALL TQL2  (LM2,N,E,Q,C,INFO)
      ELSE IF(IDIAG.EQ.3) THEN
         IF(MODE.NE.0) CALL LINEAR(A,A,N,LM2,LM4)
         CALL DSPEV('V','U',N,A,E,C,LM2,Q,INFO)
      ELSE IF(IDIAG.EQ.4) THEN
         IF(MODE.NE.0) CALL LINEAR(A,A,N,LM2,LM4)
         DUMMY = 0.0D0
         EVTOL = 2.0D0*DLAMCH('S')
         RANGE = 'I'
         IF(NROOT.GE.N) RANGE = 'A'
         CALL DSPEVX('V',RANGE,'U',N,A,DUMMY,DUMMY,1,NROOT,EVTOL,
     1               NEVS,E,C,LM2,Q,IWORK,IWORK(1+5*N),INFO)
      ELSE IF(IDIAG.EQ.5) THEN
         IF(MODE.EQ.0) CALL SQUARE(A,C,N,LM2,LM4)
         NB = ILAENV(1, 'DSYEV', 'VU', N, -1, -1, -1)
         IF(NB.GT.1) THEN
            LWORK = (NB+2)*N
         ELSE
            LWORK = 3*N
         ENDIF
         CALL DSYEV('V','U',N,C,LM2,E,Q,LWORK,INFO)
      ELSE IF(IDIAG.EQ.6) THEN
         IF(MODE.EQ.0) CALL SQUARE(A,A,N,LM2,LM4)
         DUMMY = 0.0D0
         EVTOL = 2.0D0*DLAMCH('S')
C     WE ALWAYS USE RANGE='I' IN THE NEXT FUNCTION CALL
C     TO BE CONSISTENT WITH SUBROUTINE DYNDIA.
         NB = ILAENV(1, 'DSYEVX', 'VIU', N, -1, -1, -1)
         IF(NB.GT.5) THEN
            LWORK = (NB+3)*N
         ELSE
            LWORK = 8*N
         ENDIF
         RANGE = 'I'
         IF(NROOT.GE.N) RANGE = 'A'
         CALL DSYEVX('V',RANGE,'U',N,A,LM2,DUMMY,DUMMY,1,NROOT,EVTOL,
     1               NEVS,E,C,LM2,Q,LWORK,IWORK,IWORK(1+5*N),INFO)
      ELSE IF(IDIAG.EQ.7) THEN
         IF(MODE.NE.0) CALL LINEAR(A,A,N,LM2,LM4)
         LWORK  = N*N+6*N+1
         LIWORK = 5*N+3
         CALL DSPEVD('V','U',N,A,E,C,LM2,Q,LWORK,IWORK,LIWORK,INFO)
      ELSE IF(IDIAG.EQ.8) THEN
         IF(MODE.EQ.0) CALL SQUARE(A,C,N,LM2,LM4)
         LWORK  = 2*N*N+6*N+1
         LIWORK = 5*N+3
         CALL DSYEVD('V','U',N,C,LM2,E,Q,LWORK,IWORK,LIWORK,INFO)
      ENDIF
C *** ERROR SECTION.
C     IN CASE OF AN ERROR, STOP FOR IDIAG=1-6 (STANDARD EISPACK,LAPACK).
C     TRY TO RECOVER FOR IDIAG=0,7,8 (ROUTINES EWWRSP,DSPEVD,DSYEVD)
C     AFTER RECOMPUTING THE INPUT MATRIX A AND REDEFINING IDIAG TO A
C     VALUE BETWEEN 1 AND 6 IN THE CALLING ROUTINE. THIS ATTEMPT IS
C     MOTIVATED BY THE OBSERVATION THAT OPTIONS IDIAG=1-6 ARE LESS
C     PRONE TO FAILURE THAN OPTIONS IDIAG=0,7,8. IN ADDITION, THEY
C     DO NOT USE MORE REAL SCRATCH SPACE SO THAT THERE IS NO NEED TO
C     REALLOCATE MEMORY.
      IF(INFO.NE.0) THEN
         NB6 = NBF(6)
         IDSUB  = MAX(1,IDIAG+1)
         IF(IDSUB.GT.9) IDSUB=1
         WRITE(NB6,500) SUBROU(IDSUB),INFO
         IF(IDIAG.GE.1 .AND. IDIAG.LE.6) THEN
            STOP 'DIAGON'
         ENDIF
         WRITE(NB6,510)
      ENDIF
      RETURN
  500 FORMAT(/1X,A,' FAILED WITH ERROR CODE ',I6,' IN DIAGON.')
  510 FORMAT( 1X,'RETURN TO THE CALLING ROUTINE AND TRY TO RECOVER.')
      END