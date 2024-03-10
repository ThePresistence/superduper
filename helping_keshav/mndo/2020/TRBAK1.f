      SUBROUTINE TRBAK1(NM,N,A,E,M,Z)
C
      INTEGER I,J,K,L,M,N,NM
      DOUBLE PRECISION A(NM,N),E(N),Z(NM,M)
      DOUBLE PRECISION S
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRBAK1,
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A REAL SYMMETRIC
C     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
C     SYMMETRIC TRIDIAGONAL MATRIX DETERMINED BY  TRED1.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANS-
C          FORMATIONS USED IN THE REDUCTION BY  TRED1
C          IN ITS STRICT LOWER TRIANGLE.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
C
C        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED.
C
C        Z CONTAINS THE EIGENVECTORS TO BE BACK TRANSFORMED
C          IN ITS FIRST M COLUMNS.
C
C     ON OUTPUT
C
C        Z CONTAINS THE TRANSFORMED EIGENVECTORS
C          IN ITS FIRST M COLUMNS.
C
C     NOTE THAT TRBAK1 PRESERVES VECTOR EUCLIDEAN NORMS.
C
C     THIS VERSION DATED AUGUST 1983.
C
C     MODIFICATION BY WALTER THIEL, JANUARY 1991.
C        USE UPPER TRIANGLE OF MATRIX A IN ORDER TO LOOP OVER COLUMNS.
C        LOWER TRIANGLE IS COPIED TO UPPER TRIANGLE IN DO 105 LOOP.
C        IN THE REMAINDER OF THE CODE, A(I,J) IS CHANGED TO A(J,I) ETC.
C
C     MODIFICATIONS BY WALTER THIEL, DECEMBER 1999.
C        USE CONTINUE STATEMENT TO END DO LOOPS
C        IN ORDER TO AVOID WARNINGS FROM F90 COMPILER.
C
C     NOTES BY WALTER THIEL, MAY 2003.
C        NO GAINS UPON REPLACING DO-110 LOOP BY CALL TO DDOT.
C        NO GAINS UPON REPLACING DO-120 LOOP BY CALL TO DAXPY.
C        OPENMP PARALLELIZATION FOR DO-130 LOOP APPEARS FEASIBLE.
C
      IF (M .EQ. 0) GO TO 200
      IF (N .EQ. 1) GO TO 200
      DO 105 I = 2, N
CDIR$ IVDEP
C$DIR NO_RECURRENCE
*VDIR NODEP
*VOCL LOOP,NOVREC
      DO 100 J = 1, I-1
         A(J,I) = A(I,J)
  100 CONTINUE
  105 CONTINUE
C     .......... MAIN LOOP ...................
      DO 140 I = 2, N
         L = I - 1
         IF (E(I) .EQ. ZERO) GO TO 140
C.OMP PARALLEL DO SCHEDULE(STATIC)
C.OMP.PRIVATE(S)
         DO 130 J = 1, M
            S = ZERO
            DO 110 K = 1, L
               S = S + A(K,I) * Z(K,J)
  110       CONTINUE
C           S = DDOT (L,A(1,I),1,Z(1,J),1)
C     .......... DIVISOR BELOW IS NEGATIVE OF H FORMED IN TRED1.
C                DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
            S = (S / A(L,I)) / E(I)
            DO 120 K = 1, L
               Z(K,J) = Z(K,J) + S * A(K,I)
  120       CONTINUE
C           CALL DAXPY (L,S,A(1,I),1,Z(1,J),1)
  130    CONTINUE
C.OMP END PARALLEL DO
  140 CONTINUE
  200 RETURN
      END
