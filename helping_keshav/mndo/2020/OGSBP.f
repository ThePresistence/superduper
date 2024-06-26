      SUBROUTINE OGSBP (P1,P2,LM4)
C     *
C     CALCULATE AND STORE MULLIKEN CHARGES
C     USED FOR THE GENERALIZED SOLVENT BOUNDARY POTENTIAL (GSBP)
C     *
C     P1(LM4)   ALPHA DENSITY MATRIX (I)
C     P2(LM4)   ALPHA (RHF) OR BETA (UHF) DENSITY MATRIX (I)
C     LM4       NUMBER OF BASIS FUNCTION PAIRS (I)
      USE LIMIT, ONLY: LM1, LMX, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./INDEX / INDX(LMX)
     ./PARDER/ CORE(LMZ),EHEAT(LMZ),EISOL(LMZ)
      DIMENSION P1(LM4),P2(LM4),CRG(NUMAT)

      NB88 = 88
      OPEN (NB88, FILE='mndo.gsbp.mulliken', STATUS= 'UNKNOWN')

C     CALCULATE MULLIKEN CHARGES FIRST
      DO 10 I=1,NUMAT
       NI   = NAT(I)
       IA   = NFIRST(I)
       IB   = NLAST(I)
       QI   = ZERO
       DO 20 J=IA,IB
         LL   = INDX(J)+J
         QI   = QI + P1(LL)+P2(LL)
   20  CONTINUE
       CRG(I) = CORE(NI)-QI
       WRITE(NB88,*) CRG(I)
   10 CONTINUE
      CLOSE(NB88)


      RETURN
      END


