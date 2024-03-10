      SUBROUTINE FOCK1(F,PA,PB,LM4)
C     *
C     ONE-CENTER TWO-ELECTRON CONTRIBUTIONS TO THE FOCK MATRIX.
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     F(LM4)    FOCK MATRIX (I,O).
C     PA(LM4)   RHF OR UHF-ALPHA DENSITY MATRIX (I).
C     PB(LM4)   UHF-BETA DENSITY MATRIX (I).
C     *
      USE LIMIT, ONLY: LM1, LMI, LMX, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL UHF
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./DPARM4/ REPD(52,LMZ)
     ./DPARM5/ INTIJ(243),INTKL(243),INTREP(243),INTRF1(243),INTRF2(243)
     ./FINDX1/ IP(LMI),IP1(LMI),IP2(LMI)
     ./INDEX / INDX(LMX)
     ./INDEXW/ NW(LM1)
     ./REP   / GSS(LMZ),GPP(LMZ),GSP(LMZ),GP2(LMZ),HSP(LMZ),HPP(LMZ)
     ./UHF   / UHF
      DIMENSION F(LM4),PA(LM4),PB(LM4)
      DIMENSION PP(45),W(243)
C *** LOOP OVER ATOMS.
      DO 100 II=1,NUMAT
      IA     = NFIRST(II)
      IORBS  = NLAST(II)-IA+1
      NI     = NAT(II)
C *** CODE FOR SP BASIS.
C     F(S,S) - CONTRIBUTION FROM GSS.
      KS     = INDX(IA)+IA
      F(KS)  = F(KS)+PB(KS)*GSS(NI)
      IF(IORBS.EQ.1) GO TO 100
      KX     = KS+IA+1
      KY     = KX+IA+2
      KZ     = KY+IA+3
C     F(S,S) - CONTRIBUTIONS FROM GSP AND HSP.
      QA     = PA(KX)+PA(KY)+PA(KZ)
      QAB    = QA+PB(KX)+PB(KY)+PB(KZ)
      F(KS)  = F(KS)+QAB*GSP(NI)-QA*HSP(NI)
C     F(P,P) - CONTRIBUTIONS FROM GSP,HSP,GPP,GP2,HPP.
      FROMSS =(PA(KS)+PB(KS))*GSP(NI)-PA(KS)*HSP(NI)
      QXX    = PA(KX)+PB(KX)
      QYY    = PA(KY)+PB(KY)
      QZZ    = PA(KZ)+PB(KZ)
      F(KX)  = F(KX)+FROMSS+PB(KX)*GPP(NI)+(QYY+QZZ)*GP2(NI)
     1              -(PA(KY)+PA(KZ))*HPP(NI)
      F(KY)  = F(KY)+FROMSS+PB(KY)*GPP(NI)+(QXX+QZZ)*GP2(NI)
     1              -(PA(KX)+PA(KZ))*HPP(NI)
      F(KZ)  = F(KZ)+FROMSS+PB(KZ)*GPP(NI)+(QXX+QYY)*GP2(NI)
     1              -(PA(KX)+PA(KY))*HPP(NI)
C     F(S,P) - CONTRIBUTIONS FROM GSP AND HSP.
      TWOHSP = TWO*HSP(NI)
      HSPGSP = HSP(NI)+GSP(NI)
      M      = KS+IA
      F(M)   = F(M)+(PA(M)+PB(M))*TWOHSP-PA(M)*HSPGSP
      M      = KX+IA
      F(M)   = F(M)+(PA(M)+PB(M))*TWOHSP-PA(M)*HSPGSP
      M      = KY+IA
      F(M)   = F(M)+(PA(M)+PB(M))*TWOHSP-PA(M)*HSPGSP
C     F(P,P*)- CONTRIBUTIONS FROM GP2 AND HPP.
      TWOHPP = TWO*HPP(NI)
      HPPGP2 = HPP(NI)+GP2(NI)
      M      = KX+IA+1
      F(M)   = F(M)+(PA(M)+PB(M))*TWOHPP-PA(M)*HPPGP2
      M      = KY+IA+1
      F(M)   = F(M)+(PA(M)+PB(M))*TWOHPP-PA(M)*HPPGP2
      M      = KY+IA+2
      F(M)   = F(M)+(PA(M)+PB(M))*TWOHPP-PA(M)*HPPGP2
C *** CODE FOR SPD BASIS.
C     CONTRIBUTIONS FROM INTEGRALS INVOLVING D ORBITALS.
C     IN THE RHF CASE, THE EXCHANGE CONTRIBUTIONS ARE INCLUDED
C     BY USING MODIFIED RAFFENETTI-TYPE INTEGRALS.
      IF(IORBS.EQ.4) GO TO 100
      IF(UHF) THEN
         DO 10 INTG=1,243
         W(INTG) = REPD(INTREP(INTG),NI)
   10    CONTINUE
      ELSE
         DO 20 INTG=1,243
         INT1 = INTRF1(INTG)
         INT2 = INTRF2(INTG)
         W(INTG) = REPD(INTREP(INTG),NI)
         IF(INT1.GT.0) W(INTG) = W(INTG)-PT25*REPD(INT1,NI)
         IF(INT2.GT.0) W(INTG) = W(INTG)-PT25*REPD(INT2,NI)
   20    CONTINUE
      ENDIF
C     ONE-CENTER MATRIX ELEMENTS OF TOTAL DENSITY.
      NWII   = NW(II)-1
      DO 30 KLW=1,45
      KL     = KLW+NWII
      PP(KLW)= PA(IP(KL))+PB(IP(KL))
      IF(IP1(KL).NE.IP2(KL)) PP(KLW)=PP(KLW)*TWO
   30 CONTINUE
C     COULOMB CONTRIBUTIONS TO THE FOCK MATRIX.
      DO 40 INTG=1,243
      IJ     = INTIJ(INTG)+NWII
      KLW    = INTKL(INTG)
      F(IP(IJ)) = F(IP(IJ))+PP(KLW)*W(INTG)
   40 CONTINUE
      IF(.NOT.UHF) GO TO 100
C     EXCHANGE CONTRIBUTIONS TO THE FOCK MATRIX IN THE UHF CASE.
      DO 50 INTG=1,243
      IJ     = INTIJ(INTG)+NWII
      KL     = INTKL(INTG)+NWII
      I      = IP1(IJ)
      J      = IP2(IJ)
      K      = IP1(KL)
      L      = IP2(KL)
      IF(I.LT.L) GO TO 50
      IF(I.GE.K) THEN
         IK  = INDX(I)+K
         JL  = INDX(MAX(J,L))+MIN(J,L)
         F(IK) = F(IK)-PA(JL)*W(INTG)
      ENDIF
      IF(I.GE.L .AND. K.NE.L) THEN
         IL  = INDX(I)+L
         JK  = INDX(MAX(J,K))+MIN(J,K)
         F(IL) = F(IL)-PA(JK)*W(INTG)
      ENDIF
      IF(J.GE.K .AND. I.NE.J) THEN
         JK  = INDX(J)+K
         IL  = INDX(MAX(I,L))+MIN(I,L)
         F(JK) = F(JK)-PA(IL)*W(INTG)
      ENDIF
      IF(J.GE.L .AND. I.NE.J .AND. K.NE.L) THEN
         JL  = INDX(J)+L
         IK  = INDX(MAX(I,K))+MIN(I,K)
         F(JL) = F(JL)-PA(IK)*W(INTG)
      ENDIF
   50 CONTINUE
  100 CONTINUE
      RETURN
      END