      SUBROUTINE EWWRSP (N,NVECT,NV,A,B,IWORK,ROOT,VECT,IERR)
C     *
C     DIAGONALIZATION PROCEDURE USING EISPACK-BASED ROUTINES.
C     *
C     ORIGINAL AUTHOR:  S.T.ELBERT, AMES LABORATORY-USDOE, JUNE 1985.
C     REVISED: MAY 1987.
C     THE ORIGINAL ROUTINE REFERS TO A SYMMETRIC PACKED MATRIX.
C     MODIFICATIONS BY WALTER THIEL IN JUNE 1990 AND JANUARY 1991
C     FOR INTEGRATION INTO THE MNDO PROGRAM.
C
C     PURPOSE -
C        FINDS (ALL) EIGENVALUES  AND  (SOME OR ALL) EIGENVECTORS
C        OF A SYMMETRIC MATRIX.
C
C     METHOD -
C        THE METHOD AS PRESENTED IN THIS ROUTINE CONSISTS OF FOUR STEPS:
C        FIRST, THE INPUT MATRIX IS REDUCED TO TRIDIAGONAL FORM BY THE
C        HOUSEHOLDER TECHNIQUE (ORTHOGONAL SIMILARITY TRANSFORMATIONS).
C        SECOND, THE ROOTS ARE LOCATED USING THE RATIONAL QL METHOD.
C        THIRD, THE VECTORS OF THE TRIDIAGONAL FORM ARE EVALUATED BY THE
C        INVERSE ITERATION TECHNIQUE.  VECTORS FOR DEGENERATE OR NEAR-
C        DEGENERATE ROOTS ARE FORCED TO BE ORTHOGONAL.
C        FOURTH, THE TRIDIAGONAL VECTORS ARE ROTATED TO VECTORS OF THE
C        ORIGINAL ARRAY.
C
C        THE ROUTINES USED ARE MODIFICATIONS OF THE EISPACK 3
C        ROUTINES TRED1, TQLRAT, TINVIT AND TRBAK1
C
C     REFERENCE:
C        STEPHEN T. ELBERT,
C        THEOR. CHIM. ACTA (1987) 71:169-186.
C
C        FOR FURTHER DETAILS, SEE EISPACK USERS GUIDE, B. T. SMITH
C        ET AL, SPRINGER-VERLAG, LECTURE NOTES IN COMPUTER SCIENCE,
C        VOL. 6, 2-ND EDITION, 1976.  ANOTHER GOOD REFERENCE IS
C        THE SYMMETRIC EIGENVALUE PROBLEM BY B. N. PARLETT
C        PUBLISHED BY PRENTICE-HALL, INC., ENGLEWOOD CLIFFS, N.J. (1980)
C
C     ON ENTRY -
C        N     - INTEGER
C                ORDER OF MATRIX A.
C        NVECT - INTEGER
C                NUMBER OF VECTORS DESIRED.  0.LE.NVECT AND NVECT.LE.N
C        NV    - INTEGER
C                ROW DIMENSION OF A AND VECT IN CALLING ROUTINE. N.LE.NV
C        A     - WORKING PRECISION (NV,NV)
C                INPUT MATRIX, 2-D SYMMETRIC MATRIX
C        B     - WORKING PRECISION (NV,9)
C                SCRATCH ARRAY, 9*NV ELEMENTS
C        IWORK - INTEGER
C                SCRATCH ARRAY, NV ELEMENTS
C
C     ON EXIT  -
C        A     - DESTROYED.  NOW HOLDS REFLECTION OPERATORS.
C        ROOT  - WORKING PRECISION (N)
C                ALL EIGENVALUES IN ASCENDING ORDER.
C        VECT  - WORKING PRECISION (NV,NVECT)
C                EIGENVECTORS FOR ROOT(1), ..., ROOT(NVECT).
C        IERR  - INTEGER
C                = 0 IF NO ERROR DETECTED,
C                = K IF ITERATION FOR K-TH EIGENVALUE FAILED,
C                =-K IF ITERATION FOR K-TH EIGENVECTOR FAILED.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /NBFILE/ NBF(20)
      DIMENSION A(NV,N),B(N,8),ROOT(NVECT),VECT(NV,NVECT)
      DIMENSION IWORK(N)
C
C *** FILE NUMBERS.
      NB6    = NBF(6)
C *** CHECK FOR SIMPLE ERRORS
      IERR = N-1
      IF(N.LE.0) THEN
         WRITE(NB6,910)
         WRITE(NB6,900) N,NVECT,NV,IERR
         STOP 'EWWRSP'
      ENDIF
      IF(NV.LT.N) THEN
         WRITE(NB6,920)
         WRITE(NB6,900) N,NVECT,NV,IERR
         STOP 'EWWRSP'
      ENDIF
      IERR = N+1
C *** REDUCE SYMMETRIC MATRIX A TO TRIDIAGONAL FORM
      CALL TRED1(NV,N,A,B(1,1),B(1,2),B(1,3))
C *** FIND ALL EIGENVALUES OF TRIDIAGONAL MATRIX
      CALL EQLRAT(N,B(1,1),B(1,2),B(1,3),ROOT,IWORK,IERR,B(1,4))
      IF(IERR.NE.0) THEN
         WRITE(NB6,930) IERR
         WRITE(NB6,900) N,NVECT,NV,IERR
         RETURN
      ENDIF
      IF(NVECT.LE.0) RETURN
C *** FIND EIGENVECTORS OF TRI-DIAGONAL MATRIX VIA INVERSE ITERATION
      IERR   = 6
      B(1,3) = 0.0D0
      CALL EINVIT(NV,N,B(1,1),B(1,2),B(1,3),NVECT,ROOT,IWORK,
     1            VECT,IERR,B(1,4),B(1,5),B(1,6),B(1,7),B(1,8))
      IF(IERR.NE.0) THEN
         WRITE(NB6,940) -IERR
         WRITE(NB6,900) N,NVECT,NV,IERR
         RETURN
      ENDIF
C *** FIND EIGENVECTORS OF SYMMETRIC MATRIX VIA BACK TRANSFORMATION
      CALL TRBAK1(NV,N,A,B(1,2),NVECT,VECT)
      RETURN
  900 FORMAT(/1X,'*** EWWRSP PARAMETERS ***',
     1       /1X,'***      N = ',I8,' ***',
     2       /1X,'***  NVECT = ',I8,' ***',
     3       /1X,'***     NV = ',I8,' ***',
     4       /1X,'***   IERR = ',I8,' ***')
  910 FORMAT(1X,'VALUE OF N IS LESS THAN OR EQUAL ZERO')
  920 FORMAT(1X,'NV IS LESS THAN N')
  930 FORMAT(1X,'EQLRAT HAS FAILED TO CONVERGE FOR ROOT',I5)
  940 FORMAT(1X,'EINVIT HAS FAILED TO CONVERGE FOR VECTOR',I5)
      END
