      SUBROUTINE TRDEF (AMS,AMSUM,B,PC,PM,I3N,N,NEG,NTR)
C     *
C     DEFINITION OF TRANSLATIONAL AND ROTATIONAL MODES
C     WITH RESPECT TO PRINCIPAL AXES.
C     *
C     NOTATION. I=INPUT, O=OUTPUT.
C     AMS       ATOMIC MASSES IN AMU (I).
C     AMSUM     MOLECULAR WEIGHT IN AMU (I).
C     B         MODES IN MASS-WEIGHTED CARTESIAN COORDINATES,
C               WITH RESPECT TO PRINCIPAL AXES (O).
C     PC        CARTESIAN COORDINATES FOR PRINCIPAL AXES IN ANGSTROM(I).
C     PM        PRINCIPAL MOMENTS OF INERTIA IN AMU*ANGSTROM**2 (I).
C     I3N       NUMBER OF CARTESIAN COORDINATES (I).
C     N         NUMBER OF ATOMS (I).
C     NEG       NUMBER OF MODES WITH NEGATIVE EIGENVALUES (I).
C     NTR       NUMBER OF TRANSLATIONAL AND ROTATIONAL MODES (I).
C               NTR=5 FOR LINEAR MOLECULES, NTR=6 OTHERWISE.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AMS(N),B(I3N,I3N),PC(3,N),PM(3)
      DIMENSION JJJ(3),RPM(3)
      DATA JJJ/2,3,1/
      DATA ZERO,ONE/0.0D0,1.0D0/
C     *
C *** TRANSLATIONAL MODES.
      IA     = NEG+1
      IB     = NEG+3
      DO 20 I=IA,IB
      DO 10 J=1,I3N
      B(J,I) = ZERO
   10 CONTINUE
   20 CONTINUE
      RSQM   = ONE/SQRT(AMSUM)
      DO 40 K=1,N
      KA     = 3*K-3
      SQAMS  = SQRT(AMS(K))
      DO 30 I=IA,IB
      J      = KA+I-IA+1
      B(J,I) = SQAMS*RSQM
   30 CONTINUE
   40 CONTINUE
C     *
C     ROTATIONAL MODES.
      IA     = NEG+4
      IB     = NEG+NTR
      RPM(1) = ZERO
      RPM(2) = ZERO
      RPM(3) = ZERO
      IF(PM(1).GT.ZERO) RPM(1) = ONE/SQRT(PM(1))
      IF(PM(2).GT.ZERO) RPM(2) = ONE/SQRT(PM(2))
      IF(PM(3).GT.ZERO .AND. NTR.EQ.6) RPM(3) = ONE/SQRT(PM(3))
      DO 60 K=1,N
      KA     = 3*K-3
      SQAMS  = SQRT(AMS(K))
      DO 50 I=IA,IB
      II     = I-IA+1
      JJ     = JJJ(II)
      KK     = 6-II-JJ
      B(KA+II,I) = ZERO
      B(KA+JJ,I) =-SQAMS*RPM(II)*PC(KK,K)
      B(KA+KK,I) = SQAMS*RPM(II)*PC(JJ,K)
   50 CONTINUE
   60 CONTINUE
      RETURN
      END
