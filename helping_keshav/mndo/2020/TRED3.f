      SUBROUTINE TRED3(N,NV,A,D,E,E2)
C
      INTEGER I,J,K,L,N,II,IZ,JK,NV,JM1
      DOUBLE PRECISION A(NV),D(N),E(N),E2(N)
      DOUBLE PRECISION F,G,H,HH,SCALE
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED3,
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX, STORED AS
C     A ONE-DIMENSIONAL ARRAY, TO A SYMMETRIC TRIDIAGONAL MATRIX
C     USING ORTHOGONAL SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT
C
C        N IS THE ORDER OF THE MATRIX.
C
C        NV MUST BE SET TO THE DIMENSION OF THE ARRAY PARAMETER A
C          AS DECLARED IN THE CALLING PROGRAM DIMENSION STATEMENT.
C
C        A CONTAINS THE LOWER TRIANGLE OF THE REAL SYMMETRIC
C          INPUT MATRIX, STORED ROW-WISE AS A ONE-DIMENSIONAL
C          ARRAY, IN ITS FIRST N*(N+1)/2 POSITIONS.
C
C     ON OUTPUT
C
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL
C          TRANSFORMATIONS USED IN THE REDUCTION.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO.
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.
C
C     THIS VERSION DATED AUGUST 1983.
C
C     MODIFICATIONS BY WALTER THIEL, DECEMBER 1999.
C        USE CONTINUE STATEMENT TO END DO LOOPS
C        IN ORDER TO AVOID WARNINGS FROM F90 COMPILER.
C
C     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
      DO 300 II = 1, N
         I = N + 1 - II
         L = I - 1
         IZ = (I * L) / 2
         H = ZERO
         SCALE = ZERO
         IF (L .LT. 1) GO TO 130
C     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
         DO 120 K = 1, L
            IZ = IZ + 1
            D(K) = A(IZ)
            SCALE = SCALE + ABS(D(K))
  120    CONTINUE
         IF (SCALE .NE. ZERO) GO TO 140
  130    E(I) = ZERO
         E2(I) = ZERO
         GO TO 290
  140    DO 150 K = 1, L
            D(K) = D(K) / SCALE
            H = H + D(K) * D(K)
  150    CONTINUE
C
         E2(I) = SCALE * SCALE * H
         F = D(L)
         G = -SIGN(SQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         D(L) = F - G
         A(IZ) = SCALE * D(L)
         IF (L .EQ. 1) GO TO 290
         JK = 1
         DO 240 J = 1, L
            F = D(J)
            G = ZERO
            JM1 = J - 1
            IF (JM1 .LT. 1) GO TO 220
            DO 200 K = 1, JM1
               G = G + A(JK) * D(K)
               E(K) = E(K) + A(JK) * F
               JK = JK + 1
  200       CONTINUE
  220       E(J) = G + A(JK) * F
            JK = JK + 1
  240    CONTINUE
C     .......... FORM P ..........
         F = ZERO
         DO 245 J = 1, L
            E(J) = E(J) / H
            F = F + E(J) * D(J)
  245    CONTINUE
         HH = F / (H + H)
C     .......... FORM Q ..........
         DO 250 J = 1, L
             E(J) = E(J) - HH * D(J)
  250    CONTINUE
         JK = 1
C     .......... FORM REDUCED A ..........
         DO 280 J = 1, L
            F = D(J)
            G = E(J)
            DO 260 K = 1, J
               A(JK) = A(JK) - F * E(K) - G * D(K)
               JK = JK + 1
  260       CONTINUE
  280    CONTINUE
  290    D(I) = A(IZ+1)
         A(IZ+1) = SCALE * SQRT(H)
  300 CONTINUE
      RETURN
      END
