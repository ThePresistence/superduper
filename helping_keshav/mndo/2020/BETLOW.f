      SUBROUTINE BETLOW (BORT,S,BETA,LM2,LM4)
C     *
C     BETLOW PERFORMES LOEWDIN TRANSFORMATION OF THE RESONANCE
C     INTEGRALS (IN EV) AS DEFINED IN THE OM4 METHOD.
C     BETLOW IS CALLED ONCE FOR A GIVEN GEOMETRY.
C     *
C     NOTATION. I=INPUT, O=OUTPUT.
C     BORT(LM4) ORTHOGONALIZED MATRIX OF RESONANCE INTEGRALS (O).
C     S(LM4)    OVERLAP INTEGRALS (I).
C     BETA(LM4) RESONANCE INTEGRALS (I).
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./ORBITS/ NUMB,N,NMOS,NALPHA,NBETA
C    ./NBFILE/ NBF(20)
      DIMENSION BORT(LM4),S(LM4),BETA(LM4)
C     AUTOMATIC SCRATCH ARRAYS FOR MATRIX OPERATIONS.
      DIMENSION A(LM2,LM2),B(LM2,LM2),C(LM2,LM2),D(LM2)
C *** TRANSFER OVERLAP MATRIX TO BORT.
C     THE CONTENTS OF THIS ARRAY WILL BE OVERWRITTEN IN TDIAG.
      DO 10 I=1,LM4
      BORT(I) = S(I)
   10 CONTINUE
C *** DIAGONAL ELEMENTS OF OVERLAP MATRIX MUST BE EQUAL TO ONE.
      DO 20 I=1,N
      II     = (I*(I+1))/2
      BORT(II) = ONE
   20 CONTINUE
C *** DIAGONALIZE OVERLAP MATRIX.
      CALL TDIAG (BORT,C,D,A,LM4,LM2,N,N)
C *** COMPUTE INVERSE SQUARE ROOT OF EIGENVALUES.
      DO 30 I=1,N
      D(I)   = ONE/SQRT(ABS(D(I)))
   30 CONTINUE
C *** COMPUTE S**(-1/2) MATRIX.
      DO 55 I=1,N
      DO 50 J=1,I
      A(I,J)  = ZERO
      DO 40 K=1,N
      A(I,J)  = A(I,J)+C(I,K)*D(K)*C(J,K)
   40 CONTINUE
      IF(J.NE.I) A(J,I) = A(I,J)
   50 CONTINUE
   55 CONTINUE
C *** TEST: S-1/2 * S * S-1/2 = 1
C *** BRING LINEARLY PACKED OVERLAP MATRIX IN SQUARE FORM
C     CALL SQUARE(S,B,N,LM2,LM4)
C     DO 60 I=1,N
C     B(I,I)  = ONE
C  60 CONTINUE
C *** FORM S-1/2 * S * S-1/2
C     CALL DGEMM('N','N',N,N,N,ONE,A,LM2,B,LM2,ZERO,C,LM2)
C     CALL DGEMM('N','N',N,N,N,ONE,C,LM2,A,LM2,ZERO,B,LM2)
C     CALL MATPRT(B,N,N,LM2,LM2)
C     DIFF = ZERO
C     DIFF1 = ZERO
C     DO 75 I=1,N
C     DIFF = DIFF + ABS(ONE-B(I,I))
C     DO 70 J=1,I-1
C     DIFF1 = DIFF1 + ABS(B(I,J))
C  70 CONTINUE
C  75 CONTINUE
C     print 700,diff,diff1
C 700 format(1x,'diagonal deviation',f16.15,
C    1     ' off-diagonal deviation',f16.15)
C     NB6 = NBF(6)
C     WRITE(NB6,*) '****** S-1/2 ******'
C     CALL MATPRT(A,N,N,LM2,LM2)
C *** BRING LINEARLY PACKED BETA MATRIX IN SQUARE FORM
      CALL SQUARE(BETA,B,N,LM2,LM4)
C *** PERFORM LOEWDIN-TRANSFORMATION: S-1/2 * BETA * S-1/2
      CALL DGEMM('N','N',N,N,N,ONE,A,LM2,B,LM2,ZERO,C,LM2)
      CALL DGEMM('N','N',N,N,N,ONE,C,LM2,A,LM2,ZERO,B,LM2)
C *** FORM INVERSE OF OVERLAP MATRIX
C     CALL DGEMM('N','N',N,N,N,ONE,A,LM2,A,LM2,ZERO,C,LM2)
C     WRITE(NB6,*) '****** S-1 ******'
C     CALL MATPRT(C,N,N,LM2,LM2)
C     WRITE(NB6,*) '****** BORT ******'
C     CALL MATPRT(B,N,N,LM2,LM2)
C *** BRING ORTHOGONALIZED BETA MATRIX IN LINEARLY PACKED FORM
      IJ=0
      DO 85 I=1,N
      DO 80 J=1,I
      IJ     = IJ+1
      BORT(IJ)  = PT5*(B(I,J)+B(J,I))
   80 CONTINUE
   85 CONTINUE
      IJ = IJ + 1
      DO 90 I=IJ,LM4
      BORT(I) = ZERO
   90 CONTINUE
      RETURN
      END