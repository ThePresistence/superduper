      SUBROUTINE BMAT3 (B,LMQ,NPS)
C     *
C     CALCULATION OF THE B-MATRIX WHICH REPRESENTS THE INTERACTIONS
C     BETWEEN THE SOLUTE AND THE SCREENING CHARGES IN COSMO.
C     *
C     NOTATION. I=INPUT, O=OUTPUT.
C     B         INTERACTION MATRIX (O), DIMENSION B(LMQ,NPS).
C     LMQ       TOTAL NUMBER OF SOLUTE CHARGES, ELECTRONS AND CORES (I).
C     NPS       NUMBER OF SCREENING CHARGES ON COSMO SURFACE (I).
C     *
      USE LIMIT, ONLY: LM1, LMZ, LMNPS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (SMALL=1.0D-07)
      PARAMETER (SQ3=1.7320508075689D0)
      COMMON
     ./ATOMC / COORD(3,LM1)
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./COSMO2/ NCOSMO,NPSG(4)
     ./COSMO6/ CAVITA(LMZ*4),DELTA(LMZ),OMEGA(LMZ)
     ./COSMO8/ COSURF(3,LMNPS),AREASG(LMNPS),QSG(LMNPS)
     ./INOPT2/ IN2(300)
     ./NBFILE/ NBF(20)
     ./MULTIP/ DD(6,0:LMZ),PO(9,0:LMZ)
      DIMENSION B(LMQ,NPS)
      DIMENSION P(8)
C *** INITIALIZATION.
      NPRINT = IN2(72)
      IPOT   = IN2(233)
      IDEN   = 0
C *** COMPUTE THE B-MATRIX.
C     LOOP OVER ALL ATOMS.
      DO 80 I=1,NUMAT
      NI     = NAT(I)
      IORBS  = NLAST(I)-NFIRST(I)+1
      I1     = IDEN
      I2     = LMQ-NUMAT+I
C     ATOMS WITH S-BASIS
      IF(IORBS.EQ.1) THEN
        IF(IPOT.LT.3) THEN
         DO 20 IPS=1,NPS
         R2  = (COORD(1,I)-COSURF(1,IPS))**2
         R2  = R2+(COORD(2,I)-COSURF(2,IPS))**2
         R2  = R2+(COORD(3,I)-COSURF(3,IPS))**2
         R   = SQRT(R2)
         B(I1+1,IPS) = ONE/R
         B(I2,IPS)   = ONE/R
   20    CONTINUE
        ELSE
         RHO1 = PO(1,NI)*PO(1,NI)*A0*A0
         DO 30 IPS=1,NPS
         R2  = (COORD(1,I)-COSURF(1,IPS))**2
         R2  = R2+(COORD(2,I)-COSURF(2,IPS))**2
         R2  = R2+(COORD(3,I)-COSURF(3,IPS))**2
         R   = SQRT(R2)
         B(I1+1,IPS)=ONE/SQRT(R2+RHO1)
         B(I2,IPS)=(ONE+EXP(-OMEGA(NI)*(R-DELTA(NI))))/SQRT(R2+RHO1)
   30    CONTINUE
        ENDIF
        IDEN=IDEN+1
C     ATOMS WITH SP-BASIS
      ELSE IF(IORBS.EQ.4) THEN
        IF(IPOT.LT.3) THEN
         DA  = DD(2,NI)*A0
         QA  = DD(3,NI)*DD(3,NI)*A0*A0
         DO 40 IPS=1,NPS
         X11 = COSURF(1,IPS)-COORD(1,I)
         X22 = COSURF(2,IPS)-COORD(2,I)
         X33 = COSURF(3,IPS)-COORD(3,I)
         R   = SQRT(X11*X11+X22*X22+X33*X33)
         RM1 = ONE/R
         RM3 = RM1**3
         RM5 = RM1**5
         B(I1+1,IPS) = RM1
         B(I1+2,IPS) = X11*DA*RM3
         B(I1+4,IPS) = X22*DA*RM3
         B(I1+7,IPS) = X33*DA*RM3
         B(I1+3,IPS) = RM1+THREE*X11*X11*QA*RM5-QA*RM3
         B(I1+6,IPS) = RM1+THREE*X22*X22*QA*RM5-QA*RM3
         B(I1+10,IPS)= RM1+THREE*X33*X33*QA*RM5-QA*RM3
         B(I1+5,IPS) = THREE*X11*X22*QA*RM5
         B(I1+8,IPS) = THREE*X11*X33*QA*RM5
         B(I1+9,IPS) = THREE*X22*X33*QA*RM5
         B(I2,IPS)   = RM1
   40    CONTINUE
        ELSE
         DA    = DD(2,NI)*A0
         TWOQA = TWO*DD(3,NI)*A0
         RHO1  = PO(1,NI)*PO(1,NI)*A0*A0
         RHO2  = PO(2,NI)*PO(2,NI)*A0*A0
         RHO3  = PO(3,NI)*PO(3,NI)*A0*A0
         DO 50 IPS=1,NPS
         X11 = COSURF(1,IPS)-COORD(1,I)
         X22 = COSURF(2,IPS)-COORD(2,I)
         X33 = COSURF(3,IPS)-COORD(3,I)
         BB  = X11*X11+X22*X22
         R2  = BB+X33*X33
         R   = SQRT(R2)
         SQB = SQRT(BB)
         SB  = SQB/R
C        INTEGRALS IN LOCAL COORDINATES
         B(I1+1,IPS) = ONE/SQRT(R2+RHO1)
         SPZ = PT5/SQRT((R-DA)*(R-DA)+RHO2)-PT5/SQRT((R+DA)*(R+DA)+RHO2)
         TERM1 = PT25/SQRT((R+TWOQA)*(R+TWOQA)+RHO3)
         TERM2 = PT25/SQRT((R-TWOQA)*(R-TWOQA)+RHO3)
         TERM3 = PT5 /SQRT(R2+TWOQA*TWOQA+RHO3)
         PZPZ  = TERM1+TERM2-PT5/SQRT(R2+RHO3)+ONE/SQRT(R2+RHO1)
         PXPX  = TERM3-PT5/SQRT(R2+RHO3)+ONE/SQRT(R2+RHO1)
C        CONSTRUCT ROTATION MATRIX
         IF(SB.GT.SMALL) THEN
            CA = X11/SQB
            SA = X22/SQB
            CB = X33/R
         ELSE
            SA = ZERO
            SB = ZERO
            IF(X33.LT.ZERO) THEN
               CA =-ONE
               CB =-ONE
            ELSE
               CA = ONE
               CB = ONE
            ENDIF
         ENDIF
         P(1) = CA*SB
         P(2) = SA*SB
         P(3) = CB
         P(4) = CA*CB
         P(5) = SA*CB
         P(6) = -SB
         P(7) = -SA
         P(8) = CA
C        ROTATE INTEGRALS TO MOLECULAR COORDINATES
         B(I1+2,IPS) = SPZ*P(1)
         B(I1+4,IPS) = SPZ*P(2)
         B(I1+7,IPS) = SPZ*P(3)
         B(I1+3,IPS) = PZPZ*P(1)*P(1)+PXPX*(P(4)*P(4)+P(7)*P(7))
         B(I1+6,IPS) = PZPZ*P(2)*P(2)+PXPX*(P(5)*P(5)+P(8)*P(8))
         B(I1+10,IPS)= PZPZ*P(3)*P(3)+PXPX* P(6)*P(6)
         B(I1+5,IPS) = (PZPZ*P(1)*P(2)+PXPX*(P(4)*P(5)+P(7)*P(8)))
         B(I1+8,IPS) = (PZPZ*P(1)*P(3)+PXPX*P(4)*P(6))
         B(I1+9,IPS) = (PZPZ*P(2)*P(3)+PXPX*P(5)*P(6))
         B(I2,IPS) = (ONE+EXP(-OMEGA(NI)*(R-DELTA(NI))))/SQRT(R2+RHO1)
   50    CONTINUE
        ENDIF
        IDEN=IDEN+10
C     ATOMS WITH SPD-BASIS
C     NOTE: THE QUADRUPOLE LENGTHS DD(4,NI) AND DD(6,NI) ALREADY
C     INCLUDE A FACTOR OF SQRT(2). SEE SUBROUTINES DDPOHY AND REPPD.
      ELSE
         DA  = DD(2,NI)*A0
         QA  = DD(3,NI)*DD(3,NI)*A0*A0
         QB  = DD(4,NI)*DD(4,NI)*A0*A0*PT5
         DB  = DD(5,NI)*A0
         QC  = DD(6,NI)*DD(6,NI)*A0*A0*PT5
         DO 70 IPS=1,NPS
         X11 = COSURF(1,IPS)-COORD(1,I)
         X22 = COSURF(2,IPS)-COORD(2,I)
         X33 = COSURF(3,IPS)-COORD(3,I)
         R   = SQRT(X11*X11+X22*X22+X33*X33)
         RM1 = ONE/R
         RM3 = RM1**3
         RM5 = RM1**5
         B(I1+1,IPS) = RM1
         B(I1+2,IPS) = X11*DA*RM3
         B(I1+4,IPS) = X22*DA*RM3
         B(I1+7,IPS) = X33*DA*RM3
         B(I1+3,IPS) = RM1+THREE*X11*X11*QA*RM5-QA*RM3
         B(I1+6,IPS) = RM1+THREE*X22*X22*QA*RM5-QA*RM3
         B(I1+10,IPS)= RM1+THREE*X33*X33*QA*RM5-QA*RM3
         B(I1+5,IPS) = THREE*X11*X22*QA*RM5
         B(I1+8,IPS) = THREE*X11*X33*QA*RM5
         B(I1+9,IPS) = THREE*X22*X33*QA*RM5
         DIPOLX = X11*DB*RM3
         DIPOLY = X22*DB*RM3
         DIPOLZ = X33*DB*RM3
         QUADRU = QC*RM5
         B(I1+11,IPS) = 1.5D0*(X11-X22)*(X11+X22)*QB*RM5
         B(I1+16,IPS) = THREE*X11*X33*QB*RM5
         B(I1+22,IPS) = PT5*SQ3*(TWO*X33*X33-X11*X11-X22*X22)*QB*RM5
         B(I1+29,IPS) = THREE*X22*X33*QB*RM5
         B(I1+37,IPS) = THREE*X11*X22*QB*RM5
         B(I1+12,IPS) = DIPOLX
         B(I1+13,IPS) =-DIPOLY
         B(I1+14,IPS) = ZERO
         B(I1+17,IPS) = DIPOLZ
         B(I1+18,IPS) = ZERO
         B(I1+19,IPS) = DIPOLX
         B(I1+23,IPS) =-DIPOLX/SQ3
         B(I1+24,IPS) =-DIPOLY/SQ3
         B(I1+25,IPS) = TWO*DIPOLZ/SQ3
         B(I1+30,IPS) = ZERO
         B(I1+31,IPS) = DIPOLZ
         B(I1+32,IPS) = DIPOLY
         B(I1+38,IPS) = DIPOLY
         B(I1+39,IPS) = DIPOLX
         B(I1+40,IPS) = ZERO
         B(I1+15,IPS) = RM1-THREE*X33*X33*QUADRU+QC*RM3
         B(I1+21,IPS) = RM1-THREE*X22*X22*QUADRU+QC*RM3
         B(I1+28,IPS) = RM1+THREE*X33*X33*QUADRU-QC*RM3
         B(I1+36,IPS) = RM1-THREE*X11*X11*QUADRU+QC*RM3
         B(I1+45,IPS) = RM1-THREE*X33*X33*QUADRU+QC*RM3
         B(I1+20,IPS) = THREE*X11*X33*QUADRU
         B(I1+26,IPS) =-SQ3*(X11-X22)*(X11+X22)*QUADRU
         B(I1+27,IPS) = SQ3*X11*X33*QUADRU
         B(I1+33,IPS) =-THREE*X22*X33*QUADRU
         B(I1+34,IPS) = THREE*X11*X22*QUADRU
         B(I1+35,IPS) = SQ3*X22*X33*QUADRU
         B(I1+41,IPS) = ZERO
         B(I1+42,IPS) = THREE*X22*X33*QUADRU
         B(I1+43,IPS) =-TWO*SQ3*X11*X22*QUADRU
         B(I1+44,IPS) = THREE*X11*X33*QUADRU
         B(I2,IPS) = RM1
   70    CONTINUE
         IDEN = IDEN+45
      ENDIF
   80 CONTINUE
C *** DEBUG PRINT OF B-MATRIX.
      IF(NCOSMO.EQ.1 .AND. NPRINT.GT.5) THEN
         NB6 = NBF(6)
         WRITE(NB6,500)
         CALL MATPRT (B,LMQ,NPS,LMQ,NPS)
      ENDIF
      RETURN
  500 FORMAT(///1X,'COSMO B-MATRIX: SOLUTE/SURFACE INTERACTIONS',
     1          1X,'(ANGSTROM-1).'/)
      END
