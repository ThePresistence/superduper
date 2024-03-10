      SUBROUTINE ESCF (EE,N,P,F,LM4)
C     *
C     COMPUTE CONTRIBUTIONS TO THE ELECTRONIC ENERGY.
C     *
C     NOTATION. I=INPUT,O=OUTPUT.
C     EE        ELECTRONIC ENERGY (O).
C     N         NUMBER OF BASIS FUNCTIONS (I).
C     P(LM4)    DENSITY MATRIX (I).
C     F(LM4)    FOCK MATRIX (I) OR CORE HAMILTONIAN (I).
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION P(LM4),F(LM4)
C     KMAX   = (N*(N+1))/2
C     EE     = 0.0D0
C     DO 10 K=1,KMAX
C     EE     = EE + P(K)*F(K)
C  10 CONTINUE
C     DO 20 I=1,N
C     K      = (I*(I+1))/2
C     EE     = EE - 0.5D0*P(K)*F(K)
C  20 CONTINUE
      KMAX   = (N*(N+1))/2
      EE     = DDOT(KMAX,P,1,F,1)
      EEDIAG = 0.0D0
      DO 20 I=1,N
      K      = (I*(I+1))/2
      EEDIAG = EEDIAG + P(K)*F(K)
   20 CONTINUE
      EE     = EE - 0.5D0*EEDIAG
      RETURN
      END
