      SUBROUTINE OVERLP (NI,NJ,RIJ,Z,N1DEL,N2DEL)
C     *
C     OVERLAP INTEGRALS.
C     *
C     NOTATION. I=INPUT, O=OUTPUT.
C     NI,NJ         ATOMIC NUMBERS OF ATOMS I,J (I).
C     RIJ           INTERNUCLEAR DISTANCE IN ATOMIC UNITS (I).
C     Z(I)          OVERLAP INTEGRALS (O).
C     N1DEL,N2DEL   INCREMENTS FOR MAIN QUANTUM NUMBERS OF ATOMS I,J:
C                   NONZERO ONLY FOR DERIVATIVES (I).
C     *
      USE LIMIT, ONLY: LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (RT3=0.57735026918962576451D0)
C     PARAMETER (RT3=1.0D0/SQRT(3.0D0))
      LOGICAL DIFF
      COMMON
     ./CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./DNBND / III(LMZ),IIID(LMZ)
     ./DPARM / LORBS(LMZ)
     ./DPARM1/ UDD(LMZ),ZD(LMZ),BETAD(LMZ)
     ./PAROPT/ USS(LMZ),UPP(LMZ),ZS(LMZ),ZP(LMZ),BETAS(LMZ),BETAP(LMZ),
     .         ALP(LMZ)
C    ./NBFILE/ NBF(20)
      DIMENSION Z(14),A(15),B(15)
C *** FILE NUMBERS.
C     NB6    = NBF(6)
C *** INITIALIZE OVERLAP ARRAY TO ZERO.
      DO 10 I=1,14
      Z(I)   = ZERO
   10 CONTINUE
C *** CHECK FOR IMMEDIATE RETURN (OVERLAPS BELOW THRESHOLD).
      IORBS  = LORBS(NI)
      JORBS  = LORBS(NJ)
      ZSI    = ZS(NI)
      ZPI    = ZP(NI)
      ZSJ    = ZS(NJ)
      ZPJ    = ZP(NJ)
      ZIMIN  = MIN(ZSI,ZPI)
      ZJMIN  = MIN(ZSJ,ZPJ)
      IF(IORBS.GE.9) ZIMIN = MIN(ZIMIN,ZD(NI))
      IF(JORBS.GE.9) ZJMIN = MIN(ZJMIN,ZD(NJ))
      IF((PT5*(ZIMIN+ZJMIN)*RIJ).GT.BIGEXP) RETURN
C *** FURTHER INITIALIZATION.
C     N1 ,N2 : MAIN QUANTUM NUMBERS FOR S ORBITALS.
C     N1P,N2P: MAIN QUANTUM NUMBERS FOR P ORBITALS.
      DIFF   = ZSI.NE.ZPI .OR. ZSJ.NE.ZPJ
      N1     = III(NI)+N1DEL
      N2     = III(NJ)+N2DEL
      N1P    = MAX(III(NI),2)+N1DEL
      N2P    = MAX(III(NJ),2)+N2DEL
      NT     = N1+N2
C *** BRANCHING DEPENDING ON MAIN QUANTUM NUMBERS.
      IF(N1.GT.3 .OR. N2.GT.3) GO TO 80
      IF(N1.LE.N2) THEN
         II  = N2*(N2-1)/2+N1
         ISP = 2
         IPS = 3
         IOR = IORBS
         JOR = JORBS
         FAC =-ONE
         SA  = ZSI
         PA  = ZPI
         SB  = ZSJ
         PB  = ZPJ
      ELSE
         II  = N1*(N1-1)/2+N2
         ISP = 3
         IPS = 2
         IOR = JORBS
         JOR = IORBS
         FAC = ONE
         SA  = ZSJ
         PA  = ZPJ
         SB  = ZSI
         PB  = ZPI
      ENDIF
      GO TO (20,30,50,40,60,70),II
C *** THE ORDERING OF THE ELEMENTS WITHIN Z IS
C     Z(1)   = S(I)       / S(J)
C     Z(2)   = S(I)       / P-SIGMA(J)
C     Z(3)   = P-SIGMA(I) / S(J)
C     Z(4)   = P-SIGMA(I) / P-SIGMA(J)
C     Z(5)   = P-PI(I)    / P-PI(J)
C     Z(6)   = D-SIGMA(J) / S(I)
C     Z(7)   = S(J)       / D-SIGMA(I)
C     Z(8)   = D-SIGMA(J) / P-SIGMA(I)
C     Z(9)   = P-SIGMA(J) / D-SIGMA(I)
C     Z(10)  = D-PI(J)    / P-PI(I)
C     Z(11)  = P-PI(J)    / D-PI(I)
C     Z(12)  = D-SIGMA(J) / D-SIGMA(I)
C     Z(13)  = D-PI(J)    / D-PI(I)
C     Z(14)  = D-DELTA(J) / D-DELTA(I)
C *** FIRST ROW - FIRST ROW OVERLAPS
   20 CALL SET(NT,ZSI,ZSJ,A,B,RIJ)
      IF(NI.EQ.NJ) THEN
         W   = (ZSI*RIJ)**3
      ELSE
         W   = SQRT((ZSI*ZSJ*RIJ*RIJ)**3)
      ENDIF
      Z(1)   = PT25*W*(A(3)*B(1)-B(3)*A(1))
      IF(IORBS.EQ.1 .AND. JORBS.EQ.1) RETURN
C     SPECIAL CASE: THERE ARE 2P POLARIZATION FUNCTIONS.
      IF(JORBS.EQ.4) THEN
         CALL SET(NT+1,ZSI,ZPJ,A,B,RIJ)
         W   = SQRT((ZSI**3)*(ZPJ**5))*(RIJ**4)*0.125D0
         Z(2)=-W*(A(3)*B(1)-B(3)*A(1)-A(4)*B(2)+B(4)*A(2))
      ENDIF
      IF(IORBS.EQ.4) THEN
         CALL SET(NT+1,ZSJ,ZPI,A,B,RIJ)
         W   = SQRT((ZSJ**3)*(ZPI**5))*(RIJ**4)*0.125D0
         Z(3)= W*(A(3)*B(1)-B(3)*A(1)-A(4)*B(2)+B(4)*A(2))
      ENDIF
      IF(IORBS.EQ.4 .AND. JORBS.EQ.4) THEN
         CALL SET(NT+2,ZPI,ZPJ,A,B,RIJ)
         IF(NI.EQ.NJ) THEN
           W = 0.0625D0*(ZPI*RIJ)**5
         ELSE
           W = 0.0625D0*SQRT((ZPI*ZPJ*RIJ*RIJ)**5)
         ENDIF
         Z(4)= W*(B(3)*(A(5)+A(1))-A(3)*(B(5)+B(1)))
         Z(5)= PT5*W*(A(5)*(B(1)-B(3))-B(5)*(A(1)-A(3))
     1               -A(3)*B(1)+B(3)*A(1))
      ENDIF
      RETURN
C *** FIRST ROW - SECOND ROW OVERLAPS
   30 CALL SET(NT,SA,SB,A,B,RIJ)
      W      = SQRT((SA**3)*(SB**5))*(RIJ**4)*0.125D0
      Z(1)   = W*RT3*(A(4)*B(1)+B(4)*A(1)-A(3)*B(2)-B(3)*A(2))
      IF(DIFF) THEN
         CALL SET(NT,SA,PB,A,B,RIJ)
         W   = SQRT((SA**3)*(PB**5))*(RIJ**4)*0.125D0
      ENDIF
      Z(ISP) = W*(A(3)*B(1)-B(3)*A(1)-A(4)*B(2)+B(4)*A(2))*FAC
C     SPECIAL CASE: THERE ARE 2P POLARIZATION FUNCTIONS.
      IF(IOR.EQ.4) THEN
         CALL SET(NT+1,PA,SB,A,B,RIJ)
         W   = 0.0625D0*SQRT((PA*SB*RIJ*RIJ)**5)
         D   = A(4)*(B(1)-B(3))-A(2)*(B(3)-B(5))
         E   = B(4)*(A(1)-A(3))-B(2)*(A(3)-A(5))
         Z(IPS) =-W*RT3*(D+E)*FAC
         CALL SET(NT+1,ZPI,ZPJ,A,B,RIJ)
         W   = 0.0625D0*SQRT((ZPI*ZPJ*RIJ*RIJ)**5)
         Z(4)= W*(B(3)*(A(5)+A(1))-A(3)*(B(5)+B(1)))
         Z(5)= PT5*W*(A(5)*(B(1)-B(3))-B(5)*(A(1)-A(3))
     1               -A(3)*B(1)+B(3)*A(1))
      ENDIF
      IF(JOR.GE.9) GO TO 90
      RETURN
C *** FIRST ROW - THIRD ROW OVERLAPS
   40 CALL SET(NT,SA,SB,A,B,RIJ)
      W      = SQRT((SA**3)*(SB**7)/7.5D0)*(RIJ**5)*0.0625D0
      Z(1)   = W*RT3*(A(5)*B(1)-B(5)*A(1)-TWO*(A(4)*B(2)-B(4)*A(2)))
      IF(DIFF) THEN
         CALL SET(NT,SA,PB,A,B,RIJ)
         W   = SQRT((SA**3)*(PB**7)/7.5D0)*(RIJ**5)*0.0625D0
      ENDIF
      Z(ISP) = W*FAC*(A(4)*(B(1)+B(3))+B(4)*(A(1)+A(3))
     1               -B(2)*(A(3)+A(5))-A(2)*(B(3)+B(5)))
C     SPECIAL CASE: THERE ARE 2P POLARIZATION FUNCTIONS.
      IF(IOR.EQ.4) THEN
         CALL SET(NT+1,PA,SB,A,B,RIJ)
         W   = SQRT((PA**5)*(SB**7)/7.5D0)*(RIJ**6)*0.03125D0
         Z(IPS) =-W*RT3*FAC*(A(5)*(TWO*B(3)-B(1))-B(5)*(TWO*A(3)-A(1))
     1                      +A(2)*(B(6)-TWO*B(4))-B(2)*(A(6)-TWO*A(4)))
         CALL SET(NT+1,PA,PB,A,B,RIJ)
         W   = SQRT((PA**5)*(PB**7)/7.5D0)*(RIJ**6)*0.03125D0
         Z(4)= W*(-B(4)*(A(1)+A(5))-A(4)*(B(1)+B(5))
     1            +B(3)*(A(2)+A(6))+A(3)*(B(2)+B(6)))
         Z(5)= PT5*W*(A(6)*(B(1)-B(3))+B(6)*(A(1)-A(3))-A(5)*(B(2)-B(4))
     1               -B(5)*(A(2)-A(4))-A(4)*B(1)-B(4)*A(1)
     2               +A(3)*B(2)+B(3)*A(2))
      ENDIF
      IF(JOR.GE.9) GO TO 90
      RETURN
C *** SECOND ROW - SECOND ROW OVERLAPS
   50 CALL SET(NT,ZSI,ZSJ,A,B,RIJ)
      IF(NI.EQ.NJ) THEN
         W   = 0.0625D0*(ZSI*RIJ)**5
      ELSE
         W   = 0.0625D0*SQRT((ZSI*ZSJ*RIJ*RIJ)**5)
      ENDIF
      Z(1)   = W*(A(5)*B(1)+B(5)*A(1)-TWO*A(3)*B(3))/THREE
      IF(DIFF) THEN
         CALL SET(NT,ZSI,ZPJ,A,B,RIJ)
         W   = 0.0625D0*SQRT((ZSI*ZPJ*RIJ*RIJ)**5)
      ENDIF
      D      = A(4)*(B(1)-B(3))-A(2)*(B(3)-B(5))
      E      = B(4)*(A(1)-A(3))-B(2)*(A(3)-A(5))
      Z(2)   =-W*RT3*(D-E)
      IF(DIFF) THEN
         CALL SET(NT,ZPI,ZSJ,A,B,RIJ)
         W   = 0.0625D0*SQRT((ZPI*ZSJ*RIJ*RIJ)**5)
         D   = A(4)*(B(1)-B(3))-A(2)*(B(3)-B(5))
         E   = B(4)*(A(1)-A(3))-B(2)*(A(3)-A(5))
      ENDIF
      Z(3)   = W*RT3*(D+E)
      IF(DIFF) THEN
         CALL SET(NT,ZPI,ZPJ,A,B,RIJ)
         IF(NI.EQ.NJ) THEN
           W = 0.0625D0*(ZPI*RIJ)**5
         ELSE
           W = 0.0625D0*SQRT((ZPI*ZPJ*RIJ*RIJ)**5)
         ENDIF
      ENDIF
      Z(4)   = W*(B(3)*(A(5)+A(1))-A(3)*(B(5)+B(1)))
      Z(5)   = PT5*W*(A(5)*(B(1)-B(3))-B(5)*(A(1)-A(3))
     1               -A(3)*B(1)+B(3)*A(1))
      IF(IORBS.GE.9 .OR. JORBS.GE.9) GO TO 90
      RETURN
C *** SECOND ROW - THIRD ROW OVERLAPS
   60 CALL SET(NT,SA,SB,A,B,RIJ)
      W      = SQRT((SA**5)*(SB**7)/7.5D0)*(RIJ**6)*0.03125D0
      Z(1)   = W*(A(6)*B(1)-A(5)*B(2)-TWO*(A(4)*B(3)-A(3)*B(4))
     1           +A(2)*B(5)-A(1)*B(6))/THREE
      IF(DIFF) THEN
         CALL SET(NT,SA,PB,A,B,RIJ)
         W   = SQRT((SA**5)*(PB**7)/7.5D0)*(RIJ**6)*0.03125D0
      ENDIF
      Z(ISP) = W*RT3*FAC*(-A(6)*B(2)+A(5)*B(1)-TWO*(A(3)*B(3)-A(4)*B(4))
     1                    -A(2)*B(6)+A(1)*B(5))
      IF(DIFF) THEN
         CALL SET(NT,PA,SB,A,B,RIJ)
         W   = SQRT((PA**5)*(SB**7)/7.5D0)*(RIJ**6)*0.03125D0
      ENDIF
      Z(IPS) = W*RT3*FAC*(A(5)*(TWO*B(3)-B(1))-B(5)*(TWO*A(3)-A(1))
     1                   +A(2)*(B(6)-TWO*B(4))-B(2)*(A(6)-TWO*A(4)))
      IF(DIFF) THEN
         CALL SET(NT,PA,PB,A,B,RIJ)
         W   = SQRT((PA**5)*(PB**7)/7.5D0)*(RIJ**6)*0.03125D0
      ENDIF
      Z(4)   = W*(-B(4)*(A(1)+A(5))-A(4)*(B(1)+B(5))
     1            +B(3)*(A(2)+A(6))+A(3)*(B(2)+B(6)))
      Z(5)   = PT5*W*(A(6)*(B(1)-B(3))+B(6)*(A(1)-A(3))-A(5)*(B(2)-B(4))
     1               -B(5)*(A(2)-A(4))-A(4)*B(1)-B(4)*A(1)
     2               +A(3)*B(2)+B(3)*A(2))
      IF(IORBS.GE.9 .OR. JORBS.GE.9) GO TO 90
      RETURN
C *** THIRD ROW - THIRD ROW OVERLAPS
   70 CALL SET(NT,ZSI,ZSJ,A,B,RIJ)
      IF(NI.EQ.NJ) THEN
         W   = ((ZSI*RIJ)**7)/480.0D0
      ELSE
         W   = SQRT((ZSI*ZSJ*RIJ*RIJ)**7)/480.0D0
      ENDIF
      Z(1)   = W*(A(7)*B(1)-THREE*(A(5)*B(3)-A(3)*B(5))-A(1)*B(7))/THREE
      IF(DIFF) THEN
         CALL SET(NT,ZSI,ZPJ,A,B,RIJ)
         W   = SQRT((ZSI*ZPJ*RIJ*RIJ)**7)/480.0D0
      ENDIF
      D      = A(6)*(B(1)-B(3))-TWO*A(4)*(B(3)-B(5))+A(2)*(B(5)-B(7))
      E      = B(6)*(A(1)-A(3))-TWO*B(4)*(A(3)-A(5))+B(2)*(A(5)-A(7))
      Z(2)   =-W*RT3*(D+E)
      IF(DIFF) THEN
         CALL SET(NT,ZPI,ZSJ,A,B,RIJ)
         W   = SQRT((ZPI*ZSJ*RIJ*RIJ)**7)/480.0D0
         D   = A(6)*(B(1)-B(3))-TWO*A(4)*(B(3)-B(5))+A(2)*(B(5)-B(7))
         E   = B(6)*(A(1)-A(3))-TWO*B(4)*(A(3)-A(5))+B(2)*(A(5)-A(7))
      ENDIF
      Z(3)   = W*RT3*(D-E)
      IF(DIFF) THEN
         CALL SET(NT,ZPI,ZPJ,A,B,RIJ)
         IF(NI.EQ.NJ) THEN
           W = ((ZPI*RIJ)**7)/480.0D0
         ELSE
           W = SQRT((ZPI*ZPJ*RIJ*RIJ)**7)/480.0D0
         ENDIF
      ENDIF
      Z(4)   = W*(A(3)*(B(7)+TWO*B(3))-A(5)*(B(1)+TWO*B(5))
     1           -B(5)*A(1)+A(7)*B(3))
      Z(5)   = PT5*W*(A(7)*(B(1)-B(3))+B(7)*(A(1)-A(3))
     1               +A(5)*(B(5)-B(3)-B(1))+B(5)*(A(5)-A(3)-A(1))
     2               +TWO*A(3)*B(3))
      IF(IORBS.GE.9 .OR. JORBS.GE.9) GO TO 90
      RETURN
C *** OVERLAPS INVOLVING HIGHER ROWS.
   80 CONTINUE
      CALL SET(N1 +N2 ,ZSI,ZSJ,A,B,RIJ)
      Z(1)   = SS(N1 ,0,0,N2 ,0,ZSI*RIJ,ZSJ*RIJ,A,B)
      IF(JORBS.GE.4) THEN
         IF(DIFF) CALL SET(N1 +N2P,ZSI,ZPJ,A,B,RIJ)
         Z(2) = SS(N1 ,0,0,N2P,1,ZSI*RIJ,ZPJ*RIJ,A,B)
      ENDIF
      IF(IORBS.GE.4) THEN
         IF(DIFF) CALL SET(N1P+N2 ,ZPI,ZSJ,A,B,RIJ)
         Z(3) = SS(N1P,1,0,N2 ,0,ZPI*RIJ,ZSJ*RIJ,A,B)
      ENDIF
      IF(IORBS.GE.4 .AND. JORBS.GE.4) THEN
         IF(DIFF) CALL SET(N1P+N2P,ZPI,ZPJ,A,B,RIJ)
         Z(4) = SS(N1P,1,0,N2P,1,ZPI*RIJ,ZPJ*RIJ,A,B)
         Z(5) = SS(N1P,1,1,N2P,1,ZPI*RIJ,ZPJ*RIJ,A,B)
      ENDIF
      IF(IORBS.LE.4 .AND. JORBS.LE.4) RETURN
C *** OVERLAPS INVOLVING D ORBITALS.
   90 CONTINUE
      ZDI    = ZD(NI)
      ZDJ    = ZD(NJ)
      N1D    = IIID(NI)+N1DEL
      N2D    = IIID(NJ)+N2DEL
C     WRITE(NB6,1001) (Z(I),I=1,5)
      IF(IORBS.GE.9 .AND. JORBS.LE.4) THEN
         CALL SET(N1D+N2 ,ZDI,ZSJ,A,B,RIJ)
         Z(6)  = SS(N1D,2,0,N2 ,0,ZDI*RIJ,ZSJ*RIJ,A,B)
         IF(JORBS.EQ.4) THEN
            CALL SET(N1D+N2P,ZDI,ZPJ,A,B,RIJ)
            Z(8)  = SS(N1D,2,0,N2P,1,ZDI*RIJ,ZPJ*RIJ,A,B)
            Z(10) = SS(N1D,2,1,N2P,1,ZDI*RIJ,ZPJ*RIJ,A,B) * THREE
         ENDIF
      ELSE IF(IORBS.LE.4 .AND. JORBS.GE.9) THEN
         CALL SET(N1 +N2D,ZSI,ZDJ,A,B,RIJ)
         Z(7)  = SS(N1 ,0,0,N2D,2,ZSI*RIJ,ZDJ*RIJ,A,B)
         IF(IORBS.EQ.4) THEN
            CALL SET(N1P+N2D,ZPI,ZDJ,A,B,RIJ)
            Z(9)  = SS(N1P,1,0,N2D,2,ZPI*RIJ,ZDJ*RIJ,A,B)
            Z(11) = SS(N1P,1,1,N2D,2,ZPI*RIJ,ZDJ*RIJ,A,B) * THREE
         ENDIF
      ELSE IF(IORBS.GE.9 .AND. JORBS.GE.9) THEN
         CALL SET(N1D+N2 ,ZDI,ZSJ,A,B,RIJ)
         Z(6)  = SS(N1D,2,0,N2 ,0,ZDI*RIJ,ZSJ*RIJ,A,B)
         CALL SET(N1D+N2P,ZDI,ZPJ,A,B,RIJ)
         Z(8)  = SS(N1D,2,0,N2P,1,ZDI*RIJ,ZPJ*RIJ,A,B)
         Z(10) = SS(N1D,2,1,N2P,1,ZDI*RIJ,ZPJ*RIJ,A,B) * THREE
         CALL SET(N1 +N2D,ZSI,ZDJ,A,B,RIJ)
         Z(7)  = SS(N1 ,0,0,N2D,2,ZSI*RIJ,ZDJ*RIJ,A,B)
         CALL SET(N1P+N2D,ZPI,ZDJ,A,B,RIJ)
         Z(9)  = SS(N1P,1,0,N2D,2,ZPI*RIJ,ZDJ*RIJ,A,B)
         Z(11) = SS(N1P,1,1,N2D,2,ZPI*RIJ,ZDJ*RIJ,A,B) * THREE
         CALL SET(N1D+N2D,ZDI,ZDJ,A,B,RIJ)
         Z(12) = SS(N1D,2,0,N2D,2,ZDI*RIJ,ZDJ*RIJ,A,B)
         Z(13) = SS(N1D,2,1,N2D,2,ZDI*RIJ,ZDJ*RIJ,A,B) * 9.0D0
         Z(14) = SS(N1D,2,2,N2D,2,ZDI*RIJ,ZDJ*RIJ,A,B)
      ENDIF
C     WRITE(NB6,1004) (Z(I),I=6,14)
      RETURN
C1001 FORMAT(10X,'SP-SP INTEGRALS',/2X,14F7.3)
C1002 FORMAT(10X,' D-SP INTEGRALS',/2X,14F7.3)
C1003 FORMAT(10X,'SP-D  INTEGRALS',/2X,14F7.3)
C1004 FORMAT(10X,' D-D  INTEGRALS',/2X,14F7.3)
      END
