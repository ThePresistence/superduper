      SUBROUTINE PPECP (IATOM,JATOM,ISHELL,JSHELL,NI,NJ,R,CORPP)
C     *
C     ECP INTEGRALS FOR SP BASIS IN LOCAL COORDINATE SYSTEM.
C     *
C     NOTATION. I=INPUT, O=OUTPUT.
C     IATOM     NUMBER OF FIRST  ATOM (I), LOCATED AT ORIGIN (0,0,0).
C     JATOM     NUMBER OF SECOND ATOM (I), LOCATED IN Z-AXIS (0,0,R).
C     ISHELL    NUMBER OF SHELL LOCATED AT IATOM (I).
C     JSHELL    NUMBER OF SHELL LOCATED AT JATOM (I).
C     NI,NJ     ATOMIC NUMBERS OF IATOM,JATOM (I).
C     R         INTERNUCLEAR DISTANCE IN AU (I).
C     CORPP     ECP INTEGRALS IN EV (O).
C     *
C     ALL ECP INTEGRALS FOR A GIVEN COMBINATION OF SHELLS ARE COMPUTED
C     WHICH CONTRIBUTE TO ONE-CENTER CORE HAMILTONIAN MATRIX ELEMENTS:
C     INTERACTIONS (ISHELL/JATOM/ISHELL), (JSHELL/IATOM/JSHELL).
C     *
C     NOTE ON COMMON BLOCKS IN PPECP AND ROUTINES CALLED BY PPECP.
C     - COMMON BLOCKS THAT CONTAIN CONSTANTS FROM INITIALIZATION:
C       CONSTF,CONSTN,DAWFCM,DERFCM,ERRFCM.
C     - COMMON BLOCKS THAT PROVIDE SPECIFIC INPUT DATA:
C       ECPCM,GAUSS1,GAUSS3,INOPT2.
C     - COMMON BLOCKS THAT ARE USED ONLY FOR INTERNAL COMMUNICATION:
C       ECPCM1,ECPCM2,ECPCM3,ECPCM4,FICMN,FJCMN,FSICMN,PSLIM.
C     *
C     THIS CODE IS RESTRICTED TO S AND SP SHELLS.
C     *
      USE LIMIT, ONLY: LMGP, LMGS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./ECPCM1/ ZETC,ZETB,CBDA(3),ZFN,ACZ
     ./ECPCM2/ IAT,JAT,II,JJ,IG,JG
     ./ECPCM4/ IJTOT,LF,LMX
     ./GAUSS1/ KSTART(LMGS),KNG(LMGS),KTYPE(LMGS),NSHELL,NBASIS
     ./GAUSS3/ EXX(LMGP),C1(LMGP),C2(LMGP),C3(LMGP)
     ./PSLIM / ALIM,ABLIM
C     SAVE /PSLIM/ ! OpenMP
      DIMENSION CORPP(10,2)
      DIMENSION G1(5)
C$OMP THREADPRIVATE (/ECPCM1/,/ECPCM2/,/ECPCM4/,/PSLIM/)
C *** INITIALIZATION.
      ALIM   = 0.317D0
      ABLIM  = 0.1D0
      IAT    = IATOM
      JAT    = JATOM
      II     = ISHELL
      JJ     = JSHELL
      IF(KTYPE(II).EQ.0) THEN
         ID  = 1
      ELSE
         ID  = 4
      ENDIF
      IF(KTYPE(JJ).EQ.0) THEN
         JD  = 1
      ELSE
         JD  = 4
      ENDIF
      IF(ID.EQ.1 .AND. JD.EQ.1) THEN
         IJTOT = 1
      ELSE
         IJTOT = 5
      ENDIF
C *** LOOP OVER THE TWO ATOMS INVOLVED.
      DO 50 K=1,2
      IF(K.EQ.1) THEN
         KAT = IAT
         KK  = II
         KD  = ID
         LAT = JAT
         NL  = NJ
      ELSE
         KAT = JAT
         KK  = JJ
         KD  = JD
         LAT = IAT
         NL  = NI
      ENDIF
      DO 10 I=1,4
      CORPP(I,K) = ZERO
   10 CONTINUE
      IF(NL.LE.2) GO TO 50
      K1     = KSTART(KK)
      K2     = K1+KNG(KK)-1
C *** CALCULATE ANGULAR COEFFICIENTS.
      CALL PPECPA (KD,KAT,LAT,NL,R)
C *** LOOP OVER PRIMITIVES OF SHELL II.
      DO 30 IG=K1,K2
      ZETC   = EXX(IG)
C *** LOOP OVER PRIMITIVES OF SHELL JJ.
      DO 20 JG=K1,IG
      ZETB   = EXX(JG)
C *** CALCULATE ECP INTEGRALS OVER PRIMITIVES.
      CALL PPECPS (KD,KAT,LAT,NL,G1)
C *** MULTIPLY BY CONTRACTION COEFFICIENTS.
      CORPP(1,K) = CORPP(1,K)+C1(IG)*C1(JG)*G1(1)
      IF(KD.GT.1) THEN
         CORPP(2,K) = CORPP(2,K)+C1(IG)*C2(JG)*G1(2)
         CORPP(3,K) = CORPP(3,K)+C2(IG)*C2(JG)*G1(4)
         CORPP(4,K) = CORPP(4,K)+C2(IG)*C2(JG)*G1(5)
      ENDIF
      IF(IG.EQ.JG) GO TO 20
      CORPP(1,K) = CORPP(1,K)+C1(JG)*C1(IG)*G1(1)
      IF(KD.GT.1) THEN
         CORPP(2,K) = CORPP(2,K)+C1(JG)*C2(IG)*G1(3)
         CORPP(3,K) = CORPP(3,K)+C2(JG)*C2(IG)*G1(4)
         CORPP(4,K) = CORPP(4,K)+C2(JG)*C2(IG)*G1(5)
      ENDIF
   20 CONTINUE
   30 CONTINUE
C *** CONVERT TO EV.
      DO 40 I=1,4
      CORPP(I,K) = CORPP(I,K)*EV
   40 CONTINUE
   50 CONTINUE
      RETURN
      END
C     ******************************************************************
      SUBROUTINE PPECPA (ID,IAT,LAT,NL,R)
C     *
C     CALCULATE ANGULAR PART OF ECP INTEGRALS.
C     *
C     INPUT DATA FROM ARGUMENT LIST.
C     ID    TYPE OF I-SHELL. S 1, P 2, D 3, SP 4.
C     IAT   ATOM WHERE SHELL I IS LOCATED.
C     LAT   NUMBER OF ECP CENTER IN LIST OF ATOMS.
C     NL    ATOMIC NUMBER OF ECP CENTER.
C     R     DISTANCE BETWEEN ATOMS (AU).
C     *
C     INPUT DATA FROM COMMON BLOCKS.
C     ICNTR FIRST  ATOM OF CURRENT PAIR (AS IN CALL PPECP).
C     JCNTR SECOND ATOM OF CURRENT PAIR (AS IN CALL PPECP).
C     *
      USE LIMIT, ONLY: LMZ, LMK, LML
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./ECPCM / CLP(LMK),ZLP(LMK),NLP(LMK),KFIRST(LML),KLAST(LML),
     .         LFIRST(LMZ),LMAX(LMZ)
     ./ECPCM1/ ZETC,ZETB,CA,BA,DA,ZFN,ACZ
     ./ECPCM2/ ICNTR,JCNTR,II,JJ,IG,JG
     ./ECPCM3/ NIJ,LIJ,LIJC,LIJB,NIJ2,NLIJ
     ./ECPCM4/ IJTOT,LF,LMX
C$OMP THREADPRIVATE (/ECPCM1/,/ECPCM2/,/ECPCM3/,/ECPCM4/)
      DIMENSION NXYZ(4)
C     NXYZ  L QUANTUM NUMBER OF A SHELL: S 0, P 1, D 2, SP 1.
      DATA NXYZ / 0, 1, 2, 1/
C *** INITIALIZATION FOR SHELLS.
      NXYZI  = NXYZ(ID)
      NIJ    = NXYZI+NXYZI+1
      NIJ2   = NIJ*NIJ
      LIJ    = 1
C *** INITIALIZATION FOR GEOMETRY.
      CAZ    = 0.0D0
      IF(IAT.EQ.JCNTR) CAZ = R
      IF(LAT.EQ.JCNTR) CAZ =-R
      ACZ    =-CAZ
      CA     = R
      BA     = R
      DA     = R
C *** LIMITS FOR L QUANTUM NUMBERS IN THE POTENTIAL.
      LMX    = LMAX(NL)-1
      LF     = 0
C *** COMPUTE DIRECTIONS OF UNIT VECTORS FOR C-A.
      ZFN    = CAZ/CA
C *** LAT (WITH ECP) DIFFERS FROM IAT
      LIJC   = LMX+NXYZI+1
      LIJB   = LMX+NXYZI+1
      NLIJ   = NIJ*LIJC*LIJB
      RETURN
      END
C     ******************************************************************
      SUBROUTINE PPECPS (ID,IAT,LAT,NL,G1)
C     *
C     CALCULATE ECP INTEGRALS OVER ALL PRIMITIVES.
C     *
C     INPUT DATA FROM ARGUMENT LIST.
C     ID    TYPE OF I-SHELL. S 1, P 2, D 3, SP 4.
C     IAT   ATOM WHERE SHELL I IS LOCATED.
C     LAT   NUMBER OF ECP CENTER IN LIST OF ATOMS.
C     NL    ATOMIC NUMBER OF ECP CENTER.
C     *
C     OUTPUT DATA
C     G1    ALL ECP INTEGRALS OVER PRIMITIVES.
C           ITYPE JTYPE AO(ISHELL) AO(JSHELL)
C     G1(1)   1     1       S          S
C     G1(2)   1     2       S        P-SIGMA
C     G1(3)   2     1     P-SIGMA      S
C     G1(4)   2     2     P-SIGMA    P-SIGMA
C     G1(5)   3     3     P-PI       P-PI
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.0D0)
      PARAMETER (A3=0.333333333333333D0)
      PARAMETER (A4=0.666666666666667D0)
      PARAMETER (FPI=12.5663706143592D0)
C     PARAMETER (FPI=4.0D0*PI)
      COMMON
     ./ECPCM1/ ZETC,ZETB,CA,BA,DA,ZFN,ACZ
     ./ECPCM2/ ICNTR,JCNTR,II,JJ,IG,JG
     ./ECPCM3/ NIJ,LIJ,LIJC,LIJB,NIJ2,NLIJ
     ./ECPCM4/ IJTOT,LF,LMX
     ./INOPT2/ IN2(300)
     ./NBFILE/ NBF(20)
      DIMENSION R1(49),R2(343,4)
      DIMENSION G1(5)
C$OMP THREADPRIVATE (/ECPCM1/,/ECPCM2/,/ECPCM3/,/ECPCM4/)
C *** INITIALIZATION.
      NPRINT = IN2(72)
      DO 10 I=1,5
      G1(I)  = ZERO
   10 CONTINUE
C *** COMPUTE RADIAL INTEGRALS
      CALL RCABK(NL,LF,LMX,R1,R2)
C *** COMBINE RADIAL AND ANGULAR PARTS
      G1(1) = R1(1)*FPI+R2(1,1)*FPI
      IF(IJTOT.EQ.5) THEN
         G1(2) = (ACZ*R1(1)  +ZFN*R1(5)  )*FPI
     1          +(ACZ*R2(1,1)+ZFN*R2(8,1))*FPI
         G1(3) = (ACZ*R1(1)  +ZFN*R1(5)  )*FPI
     1          +(ACZ*R2(1,1)+ZFN*R2(5,1))*FPI
         ACB1  =  ACZ*ACZ
         ACB5  = (ACZ+ACZ)*ZFN
         ACZZ  =  ACZ*ZFN
         G1(4) = (ACB1*R1(1)+ACB5*R1(5)+A3*R1(3)+A4*R1(9)     )*FPI
     1          +(ACB1*R2(1,1)+ACZZ*(R2(5,1)+R2(8,1))+R2(12,1))*FPI
         G1(5) = (A3*(R1(3)-R1(9)))*FPI
      ENDIF
C *** DEBUG PRINT.
      IF(NPRINT.GT.5) THEN
         NB6 = NBF(6)
         WRITE(NB6,500) IAT,IAT,LAT,II,JJ,IG,JG,NIJ2,NLIJ,ID
         WRITE(NB6,510) (R1(I),I=1,NIJ2)
         DO 20 L=0,LMX
         WRITE(NB6,520) L,(R2(I,L+1),I=1,NLIJ)
   20    CONTINUE
         WRITE(NB6,530) (G1(I),I=1,5)
      ENDIF
      RETURN
  500 FORMAT(1X,'PPECPS:',2X,12I4)
  510 FORMAT(1X,'R1    :',2X,5D14.5,/(10X,5D14.5))
  520 FORMAT(1X,'R2(L) :',I2,5D14.5,/(10X,5D14.5))
  530 FORMAT(1X,'G1    :',2X,5D14.5)
      END
C     ******************************************************************
      SUBROUTINE RCABK (NCNT,LF,LMX,R1,R2)
C     *
C     RADIAL ECP INTEGRALS OVER PRIMITIVES (KAHN FORMALISM).
C     *
C     INPUT DATA.
C     NCNT     NUMBER OF ECP CENTER
C     LF       MINIMUM L QUANTUM NUMBER IN POTENTIAL
C     LMX      MAXIMUM L QUANTUM NUMBER IN POTENTIAL
C     *
C     OUTPUT DATA.
C     R1(K)    RADIAL INTEGRALS OF TYPE 1
C              K = NIJ*LAMDA + KAPPA+1
C              K CORRESPONDS TO THE DIMENSION
C              (NIJ,NIJ) IN EXTERNAL ROUTINES
C     R2(K,L1) RADIAL INTEGRALS OF TYPE 2
C              K = NIJ*LIJC*LAMDAB+NIJ*LAMDAC+KAPPA+1
C              L1= L+1
C              K CORRESPONDS TO THE DIMENSIONS
C              (NIJ,LIJC,LIJB) OR (NIJ,LIJ)
C              IN EXTERNAL ROUTINES
C     *
      USE LIMIT, ONLY: LMZ, LMK, LML
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CSET  CONTRIBUTIONS TO RADIAL ECP INTEGRALS ARE NEGLECTED IF THE
CSET  PREFACTORS ARE LESS THAN 10**(-ITOL). DEFAULT ITOL=12.
      PARAMETER (ITOL=12)
      PARAMETER (TOL=ITOL*2.302585093D0)
      COMMON
     ./CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./ECPCM / CLP(LMK),ZLP(LMK),NLP(LMK),KFIRST(LML),KLAST(LML),
     .         LFIRST(LMZ),LMAX(LMZ)
     ./ECPCM1/ ZETC,ZETB,CA,BA,DA,ZFN,ACZ
     ./ECPCM3/ NIJ,LIJ,LIJC,LIJB,NIJ2,NLIJ
     ./FICMN / ALF,XI,XP0,XP1
     ./FJCMN / ALEF,BEIT,XXI,XPLS,XMNS,XP
     ./PSLIM / ALIM,ABLIM
C$OMP THREADPRIVATE (/ECPCM1/,/ECPCM3/,/FICMN/,/FJCMN/,/PSLIM/)
C     SAVE /FICMN/, /FJCMN/, /PSLIM/ !OpenMP
      DIMENSION R1(49),R2(343,4)
C
C *** INITIALIZATION.
      ZETCB  = ZETC+ZETB
      DO 10 I=1,NIJ2
      R1(I)  = ZERO
   10 CONTINUE
      DO 30 L1=1,LMX+1
      DO 20 I=1,NLIJ
      R2(I,L1) = ZERO
   20 CONTINUE
   30 CONTINUE
C
C *** FIRST SECTION  <U(LMAX)> TYPE OF INTEGRALS.
      LK     = LFIRST(NCNT)
      KF     = KFIRST(LK)
      KL     = KLAST(LK)
C *** CASE CA.NE.0, BA.NE.0, DA.NE.0, ICAB=4.
      ALFA   = ZETCB*DA
      XALFA  = ALFA*DA
C *** COMPUTE THE RADIAL INTEGRALS FOR ICAB=2,3,4.
      ALFI   = PT5/ALFA
      XP0    = ZERO
      IF(XALFA.LT.BIGEXP) XP0 = EXP(-XALFA)
      DO 40 K=KF,KL
      XI     = 1.0D0/SQRT(ZETCB+ZLP(K))
      ALF    = ALFA*XI
      DUM    = XALFA-ALF*ALF
      IF(DUM.GT.TOL) GO TO 40
      XP1    = EXP(-DUM)
      NLPK   = NLP(K)
      CLPK   = CLP(K)
      CALL FIPREP(NLPK,R1,NIJ,NIJ,CLPK)
   40 CONTINUE
      IF(NIJ.GT.2) CALL RECUR1(0,R1,NIJ,NIJ,ALFI)
C
C *** SECOND SECTION  <LM[U(L)-U(LMAX)]LM> INTEGRALS.
      IF(LMX.LT.0) RETURN
C *** CASE CA.NE.0, BA.NE.0, ICAB=4 OR ICAB=5.
      ALFA   = ZETC*CA
      XALFA  = ALFA*CA
      BETA   = ZETB*BA
      XBETA  = BETA*BA
      ALFBET = XALFA+XBETA
      ALFI   = PT5/ALFA
      BETI   = PT5/BETA
      XP     = ZERO
      IF(ALFBET.LT.BIGEXP) XP = EXP(-ALFBET)
      XTOL   = (ALFA+BETA)**2 / ZETCB
      DO 60 L=0,LMX
      L1     = L+1
      LK     = LFIRST(NCNT)+L1
      KF     = KFIRST(LK)
      KL     = KLAST(LK)
      DO 50 K=KF,KL
      ZETA   = ZETCB+ZLP(K)
      DUMTOL = ZLP(K)*XTOL/ZETA
      IF(DUMTOL.GT.TOL) GO TO 50
      XXI    = 1.0D0/SQRT(ZETA)
      ALEF   = ALFA*XXI
      BEIT   = BETA*XXI
      XPLS   = ZERO
      XMNS   = ZERO
      IF(ALEF*BEIT.GT.ABLIM) THEN
         DUM1 = ALFBET-(ALEF+BEIT)**2
         DUM2 = ALFBET-(ALEF-BEIT)**2
      ELSE
         DUM1 = ALFBET-ALEF*ALEF
         DUM2 = ALFBET-BEIT*BEIT
      ENDIF
      IF(DUM1.GE.BIGEXP .AND. DUM2.GE.BIGEXP) GO TO 50
      IF(DUM1.LT.BIGEXP) XPLS = EXP(-DUM1)
      IF(DUM2.LT.BIGEXP) XMNS = EXP(-DUM2)
      NLPK   = NLP(K)
      CLPK   = CLP(K)
      CALL FJPREP(NLPK,L,R2(1,L1),NIJ,LIJC,LIJB,CLPK)
   50 CONTINUE
      IF(NIJ.LT.2 .OR. (LIJC+LIJB).LT.4) GO TO 60
      CALL RECUR2(L,R2(1,L1),NIJ,LIJC,LIJB,ALFI,BETI)
   60 CONTINUE
      RETURN
      END
C     ******************************************************************
      SUBROUTINE FIPREP(NLPK,RR,NIJ,LIJ,CLPK)
C     *
C     PREPARE TABLE OF I-INTEGRALS BY USE OF A RECURSION RELATIONSHIP.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (LIMA=10)
      PARAMETER (SQPI=1.77245385090552D0)
      COMMON
     ./FICMN / ALF,XI,XP0,XP1
     ./PSLIM / ALIM,ABLIM
C$OMP THREADPRIVATE (/FICMN/,/PSLIM/)
C     SAVE /FICMN/, /PSLIM/ !OpenMP
      DIMENSION FIT(10,2),RR(NIJ,LIJ)
C *** INITIALIZATION
      NMX    = NIJ-1+NLPK
      YI     = 0.5D0*XI
      A1     = XI*ALF
      A2     = XI*YI
      X      = ALF*ALF
      IF(MOD(NLPK,2).EQ.1) GO TO 40
C *** GENERATE A TABLE OF ELEMENTS OF EVEN LATTICE
      IF(ALF.LE.ALIM) THEN
         T11 = 1.0D0
         T22 = 1.0D0/3.0D0
         S11 = T11
         S22 = T22
         DO 10 K=1,LIMA
         T11 = (T11*X*(K+K-1))/(K*(K+K+1))
         T22 = (T22*X*(K+K+1))/(K*(K+K+3))
         S11 = S11+T11
         S22 = S22+T22
   10    CONTINUE
         F11  = S11*SQPI*XP0
         F22  = S22*SQPI*XP0*2.0D0*ALF
      ELSE
         DAWFS = DAWF(ALF)
         F11   = SQPI*DAWFS*XP1/ALF
         F22   = SQPI*(ALF-DAWFS)*XP1/X
      ENDIF
      FIT(1,1) = F11*YI
      IF(NMX.LT.1) GO TO 70
      FIT(2,2) = F22*YI*YI
CDIR$ NOVECTOR
C$DIR SCALAR
      DO 20 N=2,NMX,2
      FIT(N+1,1) = A1*FIT(N,2)  +(N-1)*A2*FIT(N-1,1)
      FIT(N+2,2) = A1*FIT(N+1,1)+(N-2)*A2*FIT(N,2)
   20 CONTINUE
      GO TO 70
C *** GENERATE TABLE OF ELEMENTS OF ODD-LATTICE
   40 CONTINUE
      IF(ALF.LE.ALIM) THEN
         X   = 2.0D0*X
         T12 = 1.0D0/3.0D0
         T21 = 1.0D0
         S12 = T12
         S21 = T21
         DO 50 K=1,LIMA
         T12 = T12*X/(K+K+3)
         T21 = T21*X/(K+K+1)
         S12 = S12+T12
         S21 = S21+T21
   50    CONTINUE
         F12 = S12*XP0*2.0D0*ALF
         F21 = S21*XP0*2.0D0
      ELSE
         ERRFS = SQPI*ERRF(ALF)*XP1
         F21   = ERRFS/ALF
         F12   = (0.5D0*ERRFS/ALF-XP0)/ALF
      ENDIF
      FIT(1,2) = F12*YI
      FIT(2,1) = F21*YI*YI
C$DIR SCALAR
      DO 60 N=2,NMX,2
      FIT(N+1,2) = A1*FIT(N,1)+(N-3)*A2*FIT(N-1,2)
      FIT(N+2,1) = A1*FIT(N+1,2)+N*A2*FIT(N,1)
   60 CONTINUE
C *** ACCUMULATE RESULTS
C$DIR SCALAR
   70 CONTINUE
      DO 80 N=1,NIJ,2
      RR(N,1) = RR(N,1)+FIT(N+NLPK,1)*CLPK
   80 CONTINUE
      IF(NIJ.LE.1) RETURN
C$DIR SCALAR
      DO 90 N=2,NIJ,2
      RR(N,2) = RR(N,2)+FIT(N+NLPK,2)*CLPK
   90 CONTINUE
      RETURN
      END
C     ******************************************************************
      SUBROUTINE FJPREP(NLPK,L,RR,NIJ,LIJC,LIJB,CLPK)
C     *
C     PREPARE TABLE OF J-INTEGRALS BY USE OF A RECURSION RELATIONSHIP.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     SAVE !OpenMP
      COMMON
     ./FJCMN / ALEF,BET,XI,XPLS,XMNS,XP
     ./PSLIM / ALIM,ABLIM
C$OMP THREADPRIVATE (/FJCMN/,/PSLIM/)
      DIMENSION F(10,2,2),RR(NIJ,LIJC,LIJB)
C *** INITIALIZATION.
      LAMMAX = MAX(2,L+1)
      NMX    = NIJ-1+NLPK
CDIR$ SHORTLOOP
C$DIR MAX_TRIPS(128)
      DO 10 N=1,NMX+2
      F(N,1,1) = 0.0D0
      F(N,2,1) = 0.0D0
      F(N,1,2) = 0.0D0
      F(N,2,2) = 0.0D0
   10 CONTINUE
      XA     = XI*ALEF
      XB     = XI*BET
      XX     = XI*XI*0.5D0
      IF(MOD(NLPK,2).EQ.1) GO TO 50
C *** GENERATE A TABLE OF ELEMENTS OF THE EVEN LATTICE OF J-INTEGRALS.
      IF(ALEF*BET.LE.ABLIM) CALL SITABL(LAMMAX,-1)
      F(1,1,1) = FJ00(0)
      F(1,2,2) = FJ11(0)
      IF(NMX.LT.1) GO TO 70
      F(2,2,1) = FJ10(1)
      F(2,1,2) = FJ01(1)
      DO 40 N=2,NMX,2
      NP1    = N+1
      F(NP1,1,1) = XA*F(N,2,1)  +XB*F(N,1,2)  +(N-1)*XX*F(N-1,1,1)
      F(NP1,2,2) = XA*F(N,1,2)  +XB*F(N,2,1)  +(N-5)*XX*F(N-1,2,2)
      F(N+2,2,1) = XA*F(NP1,1,1)+XB*F(NP1,2,2)+(N-2)*XX*F(N,2,1)
      F(N+2,1,2) = XB*F(NP1,1,1)+XA*F(NP1,2,2)+(N-2)*XX*F(N,1,2)
   40 CONTINUE
      GO TO 70
C *** GENERATE A TABLE OF ELEMENTS OF THE ODD LATTICE OF J-INTEGRALS.
   50 IF(ALEF*BET.LE.ABLIM) CALL SITABL(-1,LAMMAX)
      F(1,2,1) = FJ10(0)
      F(1,1,2) = FJ01(0)
      F(2,1,1) = FJ00(1)
      F(2,2,2) = FJ11(1)
      DO 60 N=2,NMX,2
      NP1    = N+1
      F(NP1,2,1) = XA*F(N,1,1)  +XB*F(N,2,2)  +(N-3)*XX*F(N-1,2,1)
      F(NP1,1,2) = XB*F(N,1,1)  +XA*F(N,2,2)  +(N-3)*XX*F(N-1,1,2)
      F(N+2,1,1) = XA*F(NP1,2,1)+XB*F(NP1,1,2)+ N   *XX*F(N,1,1)
      F(N+2,2,2) = XA*F(NP1,1,2)+XB*F(NP1,2,1)+(N-4)*XX*F(N,2,2)
   60 CONTINUE
C *** ACCUMULATE RESULTS.
CDIR$ NOVECTOR
C$DIR SCALAR
   70 DO 80 N=1,NIJ,2
      RR(N,1,1) = RR(N,1,1) + F(N+NLPK,1,1)*CLPK
   80 CONTINUE
      IF(LIJC.GT.1 .AND. LIJB.GT.1) THEN
         DO 85 N=1,NIJ,2
         RR(N,2,2) = RR(N,2,2) + F(N+NLPK,2,2)*CLPK
   85    CONTINUE
      ENDIF
C$DIR SCALAR
      IF(LIJB.GT.1) THEN
         DO 90 N=2,NIJ,2
         RR(N,1,2) = RR(N,1,2) + F(N+NLPK,1,2)*CLPK
   90    CONTINUE
      ENDIF
      IF(LIJC.GT.1) THEN
         DO 95 N=2,NIJ,2
         RR(N,2,1) = RR(N,2,1) + F(N+NLPK,2,1)*CLPK
   95    CONTINUE
      ENDIF
      IF(L.GE.2)  RR(1,3,3) = RR(1,3,3) + FJ22(NLPK)*CLPK
      IF(L.GE.3)  RR(1,4,4) = RR(1,4,4) + FJ33(NLPK)*CLPK
      RETURN
      END
C     ******************************************************************
      SUBROUTINE SITABL(LEMAX,LOMAX)
C     *
C     PREPARE A TABLE OF SCALED I INTEGRALS FOR USE IN THE EVALUATION
C     OF THE J INTEGRALS BY A POWER SERIES.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     SAVE !OpenMP
      COMMON
     ./FICMN / ALFI,XI,XP0,XP1
     ./FJCMN / ALFJ,BETJ,XJ,XPLS,XMNS,XP
     ./FSICMN/ SI(26,4)
C$OMP THREADPRIVATE (/FICMN/,/FJCMN/,/FSICMN/)
C *** PREPARE PARAMETERS OF SCALED INTEGRALS DEPENDING ON THE SIGN
C     OF THE ALFA AND BETA VARIABLES.
      IF(ALFJ.LE.BETJ) THEN
         ALFI= BETJ
         X   = BETJ
         XI  = XJ
         XP0 = XP
         XP1 = XMNS
      ELSE
         ALFI= ALFJ
         X   = ALFJ
         XI  = XJ
         XP0 = XP
         XP1 = XPLS
      ENDIF
C *** GENERATE LATTICE OF SCALED I(N,L) INTEGRALS FOR (N+L) EVEN.
      IF(LEMAX.GE.0) THEN
         SI(1,1) = FSI0(0)
         SI(2,2) = FSI1(1)
         IF(LEMAX.GE.2) THEN
            SI(3,3) = FSI2(2)
            IF(LEMAX.GE.3) THEN
               SI(4,4) = FSI3(3)
            ENDIF
         ENDIF
         CALL RECUR(2,22,LEMAX,X)
      ENDIF
C *** GENERATE LATTICE OF SCALED I(N,L) INTEGRALS FOR (N+L) ODD.
      IF(LOMAX.GE.0) THEN
         SI(1,2) = FSI1(0)
         SI(2,1) = FSI0(1)
         SI(3,2) = FSI1(2)
         IF(LOMAX.GE.2) THEN
            SI(4,3) = FSI2(3)
            IF(LOMAX.GE.3) THEN
               SI(5,4) = FSI3(4)
            ENDIF
         ENDIF
         CALL RECUR(3,21,LOMAX,X)
      ENDIF
      RETURN
      END
C     ******************************************************************
      FUNCTION FSI(N)
C     *
C     EVALUATE THE SCALED I(N,L) INTEGRALS FROM THE ANALYTIC FORMULAS.
C     TRANSFER TO A POWER SERIES APPROACH, IF NECESSARY.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (SQPI=1.77245385090552D0)
C     SAVE ! OpenMP
C     SAVE ERRFS,DAWFS ! OpenMP
      COMMON
     ./FICMN / ALFA,XI,XP0,XP1
     ./PSLIM / ALIM,ABLIM
     ./FSICMN/ ERRFS,DAWFS
C$OMP THREADPRIVATE (/FICMN/,/PSLIM/,/FSICMN/)
C
      ENTRY FSI0(N)
      IF(ALFA.LE.ALIM) THEN
         FSI = FSIPS(N,0)
      ELSE IF(N.NE.1) THEN
         DAWFS = DAWF(ALFA)
         FSI = SQPI*DAWFS*XP1/ALFA
      ELSE
         FSI = ERRFS/ALFA
      ENDIF
      FSI0 = FSI
      RETURN
C
      ENTRY FSI1(N)
      IF(ALFA.LE.ALIM) THEN
         FSI = FSIPS(N,1)
      ELSE IF(N.LE.0) THEN
         ERRFS = SQPI*ERRF(ALFA)*XP1
         FSI = (0.5D0*ERRFS/ALFA-XP0)/ALFA
      ELSE IF(N.EQ.1) THEN
         FSI = SQPI*(ALFA-DAWFS)*XP1/(ALFA*ALFA)
      ELSE
         ALFA2 = ALFA*ALFA
         FSI = (2.0D0*ALFA*XP0+(2.0D0*ALFA2-1.0D0)*ERRFS)/ALFA2
      ENDIF
      FSI1 = FSI
      RETURN
C
      ENTRY FSI2(N)
      IF(ALFA.LE.ALIM) THEN
         FSI = FSIPS(N,2)
      ELSE IF(N.EQ.0) THEN
         ALFA2 = ALFA*ALFA
         FSI = 0.25D0*SQPI*(3.0D0*ALFA-(2.0D0*ALFA2+3.0D0)*DAWFS)*XP1/
     1         (ALFA2*ALFA)
      ELSE IF(N.EQ.1) THEN
         ALFA2 = ALFA*ALFA
         FSI =(0.5D0*(2.0D0*ALFA2-3.0D0)*ERRFS+3.0D0*ALFA*XP0)/
     1        (ALFA2*ALFA)
      ELSE IF(N.EQ.2) THEN
         ALFA2 = ALFA*ALFA
         FSI = SQPI*(ALFA*(2.0D0*ALFA2-3.0D0)+3.0D0*DAWFS)*XP1/
     1         (ALFA2*ALFA)
      ELSE IF(N.EQ.3) THEN
         ALFA2 = ALFA*ALFA
         FSI = ((ALFA2*(4.0D0*ALFA2-4.0D0)+3.0D0)*ERRFS+2.0D0*ALFA*
     1          (2.0D0*ALFA2-3.0D0)*XP0)/(ALFA2*ALFA)
      ENDIF
      FSI2 = FSI
      RETURN
C
      ENTRY FSI3(N)
      IF(ALFA.LE.ALIM) THEN
         FSI = FSIPS(N,3)
      ELSE IF(N.EQ.0) THEN
         ALFA2 = ALFA*ALFA
         FSI = (0.25D0*(2.0D0*ALFA2-5.0D0)*ERRFS+ALFA*(4.0D0*ALFA2
     1          +15.0D0)*XP0/6.0D0)/(ALFA2*ALFA2)
      ELSE IF(N.EQ.1) THEN
         ALFA2 = ALFA*ALFA
         FSI = 0.25D0*SQPI*(ALFA*(4.0D0*ALFA2-15.0D0)+3.0D0*
     1         (2.0D0*ALFA2+5.0D0)*DAWFS)*XP1/(ALFA2*ALFA2)
      ELSE IF(N.EQ.2) THEN
         ALFA2 = ALFA*ALFA
         FSI = (0.5D0*(ALFA2*(4.0D0*ALFA2-12.0D0)+15.0D0)*ERRFS+
     1         ALFA*(2.0D0*ALFA2-15.0D0)*XP0)/(ALFA2*ALFA2)
      ELSE IF(N.EQ.3) THEN
         ALFA2 = ALFA*ALFA
         FSI = SQPI*(ALFA*(ALFA2*(4.0D0*ALFA2-10.0D0)+15.0D0)-
     1         15.0D0*DAWFS)*XP1/(ALFA2*ALFA2)
      ELSE IF(N.EQ.4) THEN
         ALFA2 = ALFA*ALFA
         FSI = (2.0D0*ALFA*(ALFA2*(4.0D0*ALFA2-8.0D0) +15.0D0)*XP0+
     1         (ALFA2*(ALFA2*(8.0D0*ALFA2-12.0D0)+18.0D0)-15.0D0)
     2         *ERRFS)/(ALFA2*ALFA2)
      ENDIF
      FSI3 = FSI
      RETURN
      END
C     ******************************************************************
      FUNCTION FSIPS(N,L)
C     *
C     EVALUATE THE SCALED I(N,L) INTEGRALS BY A POWER SERIES METHOD.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (LIMA=10)
      PARAMETER (SQPI=1.77245385090552D0)
C     SAVE !OpenMP
      COMMON
     ./FICMN / ALFA,XI,XP0,XP1
C$OMP THREADPRIVATE (/FICMN/)
      DIMENSION FCTRL(7),DFCTRL(8)
      DATA FCTRL/1.0D0,1.0D0,2.0D0,6.0D0,24.0D0,120.0D0,720.0D0/
      DATA DFCTRL/1.0D0,1.0D0,3.0D0,15.0D0,105.0D0,945.0D0,10395.0D0,
     1            135135.0D0/
C
C *** TEST WHETHER (N+L) IS EVEN OR ODD INTEGER.
      NL     = N+L
      IF(MOD(NL,2).EQ.1) GO TO 30
C *** LATTICE OF EVEN (N+L) SCALED I(N,L) INTEGRALS.
      LAMBDA = NL/2
      X      = ALFA*ALFA
      L2     = L+L+1
C     SUM THE POWER SERIES.
      T      = DFCTRL(LAMBDA+1)/DFCTRL(L+2)
      SUM    = T
      DO 20 K=1,LIMA
      T      = (T*X*(K+K+NL-1))/(K*(K+K+L2))
      SUM    = SUM+T
   20 CONTINUE
      FSIPS  = SUM*SQPI*XP0*(2.0D0**LAMBDA)*(ALFA**L)
      RETURN
C *** LATTICE OF ODD (N+L) SCALED I(N,L) INTEGRALS.
   30 LAMBDA = (NL-1)/2
      X      = 2.0D0*ALFA*ALFA
      L2     = L+L+1
C     SUM THE POWER SERIES.
      T      = FCTRL(LAMBDA+1)/DFCTRL(L+2)
      SUM    = T
      DO 40 K=1,LIMA
      T      = (T*X*(K+LAMBDA))/(K*(K+K+L2))
      SUM    = SUM+T
   40 CONTINUE
      FSIPS  = SUM*XP0*(2.0D0**NL)*(ALFA**L)
      RETURN
      END
C     ******************************************************************
      SUBROUTINE RECUR (NMIN,NMAX,LMAX,X)
C     *
C     USE RECURSION RELATIONS TO GENERATE A TABLE OF SCALED I(N,L)
C     INTEGRALS.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     SAVE !OpenMP
      COMMON
     ./FSICMN/ SI(26,4)
C$OMP THREADPRIVATE (/FSICMN/)
      TX = X+X
      IF(LMAX.LE.1) THEN
         DO 10 N=NMIN,NMAX,2
         TFNM2     = N+N-4
         SI(N+1,1) = (TFNM2+2.0D0)*SI(N-1,1)+TX*SI(N,2)
         SI(N+2,2) = TFNM2*SI(N,2)+TX*SI(N+1,1)
   10    CONTINUE
      ELSE IF(LMAX.EQ.2) THEN
         DO 20 N=NMIN,NMAX,2
         TFNM2     = N+N-4
         SI(N+1,1) = (TFNM2+2.0D0)*SI(N-1,1)+TX*SI(N,2)
         SI(N+2,2) = TFNM2*SI(N,2)+TX*SI(N+1,1)
         SI(N+3,3) = TFNM2*SI(N+1,3)+TX*SI(N+2,2)
   20    CONTINUE
      ELSE
         DO 30 N=NMIN,NMAX,2
         TFNM2     = N+N-4
         SI(N+1,1) = (TFNM2+2.0D0)*SI(N-1,1)+TX*SI(N,2)
         SI(N+2,2) = TFNM2*SI(N,2)+TX*SI(N+1,1)
         SI(N+3,3) = TFNM2*SI(N+1,3)+TX*SI(N+2,2)
         SI(N+4,4) = TFNM2*SI(N+2,4)+TX*SI(N+3,3)
   30    CONTINUE
      ENDIF
      RETURN
      END
C     ******************************************************************
      FUNCTION FJ(N)
C     *
C     EVALUATE THE J(N,L1,L2) INTEGRALS FROM THEIR ANALYTIC FORMULAS.
C     USE A POWER SERIES  FOR SMALL VALUES OF THE ARGUMENTS.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (SQPI=1.77245385090552D0)
C     SAVE !OpenMP
C     SAVE TAB,EM,EP,HM,DM,DP ! OpenMP
      COMMON
     ./FJCMN / A,B,XI,XPLS,XMNS,XP
     ./PSLIM / ALIM,ABLIM
     ./FJCMN2/ TAB,EM,EP,HM,DM,DP
C$OMP THREADPRIVATE (/FJCMN/,/PSLIM/,/FJCMN2/)
C
      ENTRY FJ00(N)
      IF(A*B.LE.ABLIM) THEN
         FJ = FJPS(N,0,0)
      ELSE IF(N.NE.1) THEN
         TP = XPLS*DAWF(A+B)
         TM = XMNS*DAWF(A-B)
         DP = TP+TM
         DM = TP-TM
         TAB= FM(-1)
         FJ = SQPI*XI*(A*DM+B*DP-TAB*FM(0))/(2.0D0*TAB)
      ELSE
         FJ = SQPI*XI*XI*HM/(4.0D0*TAB)
      ENDIF
      FJ00 = FJ
      RETURN
C
      ENTRY FJ10(N)
      IF(A*B.LE.ABLIM) THEN
         FJ = FJPS(N,1,0)
      ELSE IF(N.NE.1) THEN
         TP = XPLS*ERRF(A+B)
         TM = XMNS*ERRF(A-B)
         EP = TP+TM
         EM = TP-TM
         HM = XPLS*DAWERF(A+B)-XMNS*DAWERF(A-B)
         TAB= FM(-1)
         FJ = XI*(SQPI*((1.0D0+2.0D0*(A+B)*(A-B))*HM+B*EP-A*EM)/
     1        (8.0D0*A*TAB)-XP/(4.0D0*A))
      ELSE
         FJ = SQPI*XI*XI*(2.0D0*FM(0)-DP/A)/(8.0D0*A)
      ENDIF
      FJ10 = FJ
      RETURN
C
      ENTRY FJ01(N)
      IF(A*B.LE.ABLIM) THEN
         FJ = FJPS(N,0,1)
      ELSE IF(N.NE.1) THEN
         FJ = XI*(SQPI*((1.0D0-2.0D0*(A+B)*(A-B))*HM-B*EP+A*EM)/
     1        (8.0D0*B*TAB)-XP/(4.0D0*B))
      ELSE
         FJ = SQPI*XI*XI*(2.0D0*FM(0)-DM/B)/(8.0D0*B)
      ENDIF
      FJ01 = FJ
      RETURN
C
      ENTRY FJ11(N)
      IF(A*B.LE.ABLIM) THEN
         FJ = FJPS(N,1,1)
      ELSE IF(N.NE.1) THEN
         A2 = A*A
         B2 = B*B
         FJ = SQPI*XI*(2.0D0*(A2+B2)*FM(0)-TAB*FM(1)-B2*DP/A-A2*DM/B)/
     1        (6.0D0*TAB)
      ELSE
         FJ = XI*XI*(SQPI*(A*EM+B*EP-(1.0D0+2.0D0*(A*A+B*B))*HM)/
     1        (8.0D0*TAB*TAB)-XP/(4.0D0*TAB))
      ENDIF
      FJ11 = FJ
      RETURN
C
      ENTRY FJ22(N)
      IF(A*B.LE.ABLIM) THEN
         FJ = FJPS(N,2,2)
      ELSE IF(N.LE.0) THEN
         A2 = A*A
         B2 = B*B
         A2PB2 = A2+B2
         FJ = SQPI*XI*((A2*B2-A2PB2*A2PB2)*FM(0)+0.5D0*TAB*(A2PB2+1.5D0)
     1        *FM(1)+0.5D0*B2*B2*DP/A+0.5D0*A2*A2*DM/B)/(10.0D0*A2*B2)
      ELSE IF(N.EQ.1) THEN
         A2 = A*A
         B2 = B*B
         A2PB2 = A2+B2
         A2MB2 = A2-B2
         FJ =XI*XI*(SQPI*(1.5D0*B*(A2MB2-1.5D0)*EP-1.5D0*A*(A2MB2+1.5D0)
     1       *EM+3.0D0*(A2PB2*(A2PB2+1.0D0)+0.75D0-TAB*TAB/3.0D0)*HM)/
     2       (8.0D0*TAB**3)+3.0D0*XP*(A2PB2+1.5D0)/(8.0D0*TAB*TAB))
      ELSE
         FJ = 0.25D0*SQPI*XI*XI*XI*FM(2)
      ENDIF
      FJ22 = FJ
      RETURN
C
      ENTRY FJ33(N)
      IF(A*B.LE.ABLIM) THEN
         FJ = FJPS(N,3,3)
      ELSE IF(N.LE.0) THEN
         A2 = A*A
         B2 = B*B
         A2PB2 = A2+B2
         A4 = A2*A2
         B4 = B2*B2
         FJ = SQPI*XI*(A2PB2*(A4+B4)*FM(0)+A*B*(A2*B2-A2PB2*(A2PB2+
     1        1.5D0))*FM(1)+2.50D0*A2*B2*FM(2)-0.5D0*(A4*A2*DM/B+
     2        B4*B2*DP/A))/(14.0D0*A2*A*B2*B)
      ELSE IF(N.EQ.1) THEN
         A2 = A*A
         B2 = B*B
         A2PB2 = A2+B2
         A2MB2 = A2-B2
         A2TB2 = A2*B2
         FJ =XI*XI*(SQPI*(B*(A2TB2-0.625D0*(3.0D0*A2-2.0D0*B2-3.75D0+
     1       A2PB2*A2MB2))*EP+A*(A2TB2-0.625D0*(3.0D0*B2-2.0D0*A2-
     2       3.75D0-A2PB2*A2MB2))*EM+(3.0D0*A2TB2*(A2PB2+0.5D0)-
     3       1.25D0*(1.875D0+A2PB2*(2.25D0+A2PB2*(1.5D0+A2PB2))))*HM)/
     4       (32.0D0*A2TB2*A2TB2)-0.625D0*(A2PB2*(A2PB2+2.0D0)-
     5       16.0D0*A2*B2/15.0D0+3.75D0)*XP/(TAB**3))
      ELSE
         FJ = SQPI*XI*XI*XI*(TAB*FM(1)-5.0D0*FM(2))/(4.0D0*TAB)
      ENDIF
      FJ33 = FJ
      RETURN
      END
C     ******************************************************************
      FUNCTION FJPS(N,LALF,LBET)
C     *
C     EVALUATE J-INTEGRALS BY USE OF A POWER SERIES.
C     SET UP THE PARAMETERS IN THE POWER SERIES DEPENDING ON THE
C     RELATIVE VALUES OF THE ALFA AND BETA VARIABLES.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (LIMAB=10)
C     SAVE !OpenMP
      COMMON
     ./FJCMN / ALF,BET,XI,XPLS,XMNS,XP
     ./FSICMN/ SI(26,4)
C$OMP THREADPRIVATE (/FJCMN/,/FSICMN/)
      DIMENSION DFCTRL(4)
      DATA DFCTRL/1.0D0,3.0D0,15.0D0,105.0D0/
C
      IF(ALF.GT.BET) THEN
         X   = BET*BET
         T   = BET**LBET/DFCTRL(LBET+1)
         L1  = LALF+1
         L2  = LBET+N+1
         L3  = LBET+LBET+1
      ELSE
         X   = ALF*ALF
         T   = ALF**LALF/DFCTRL(LALF+1)
         L1  = LBET+1
         L2  = LALF+N+1
         L3  = LALF+LALF+1
      ENDIF
C *** START THE POWER SERIES SUM.
      SUM    = T*SI(L2,L1)
      DO 10 K=1,LIMAB
      T      = T*X/(2*K*(K+K+L3))
      SUM    = SUM+T*SI(K+K+L2,L1)
   10 CONTINUE
      FJPS = SUM*(0.5D0*XI)**(N+1)
      RETURN
      END
C     ******************************************************************
      FUNCTION FM(L)
C     *
C     CALCULATE THE MODIFIED SPHERICAL BESSEL FUNCTION
C     MULTIPLIED BY A GAUSSIAN EXPONENTIAL.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     SAVE !OpenMP
C     SAVE T,EST,ECT ! OpenMP
      COMMON
     ./FJCMN / A,B,XI,XPLS,XMNS,XP
     ./FMCMN / T,EST,ECT
C$OMP THREADPRIVATE (/FJCMN/,/FMCMN/)
      LP2 = L+2
      GO TO (10,20,30,40,50),LP2
   10 T   = 2.0D0*A*B
      EST = 0.5D0*(XPLS-XMNS)
      ECT = 0.5D0*(XPLS+XMNS)
      FM  = T
      RETURN
   20 FM = EST/T
      RETURN
   30 FM = (-EST/T+ECT)/T
      RETURN
   40 T2 = T*T
      FM = (3.0D0/T2+1.0D0)*EST/T-3.0D0*ECT/T2
      RETURN
   50 T2 = T*T
      FM = -(15.0D0/T2+6.0D0)*EST/T2+(15.0D0/T2+1.0D0)*ECT/T
      RETURN
      END
C     ******************************************************************
      FUNCTION ERRF(Y)
C     *
C     EVALUATE THE ERROR FUNCTION (AN ODD-PARITY FUNCTION).
C     USE A PIECEWISE CHEBYSHEV POLYNOMIAL BY DEFAULT.
C     USE THE ASYMPTOTIC VALUE IF NEEDED.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     SAVE !OpenMP
      COMMON
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./ERRFCM/ C(142),IFIRST(20),ILAST(20),H
      X      = ABS(Y)
      IF(X.LT.4.86D0) THEN
         NX  = INT(X/H)
         NX  = NX+1
         IF  = IFIRST(NX)
         IL  = ILAST(NX)
         T   = C(IL)
         KL  = IL-IF
         DO 10 K=1,KL
         T   = C(IL-K)+X*T
   10    CONTINUE
         ERRF = T
      ELSE
         ERRF = ONE
      ENDIF
      IF(Y.LT.ZERO) ERRF=-ERRF
      RETURN
      END
C     ******************************************************************
      FUNCTION DAWF (Y)
C     *
C     EVALUATE THE DAWSON FUNCTION.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     SAVE !OpenMP
      COMMON
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./DAWFCM/ C(249),IFIRST(40),ILAST(40),H
      X      = ABS(Y)
      IF(X.LT.10.0D0) THEN
         NX  = INT(X/H)
         NX  = NX+1
         IF  = IFIRST(NX)
         IL  = ILAST(NX)
         T   = C(IL)
         KL  = IL-IF
         DO 10 K=1,KL
         T   = C(IL-K)+X*T
   10    CONTINUE
         DAWF= T
      ELSE
         TXT = PT5/(X*X)
         TX  = TXT*X
         DAWF= TX*(ONE+TXT*(ONE+TXT*(THREE+TXT*(15.0D0+105.0D0*TXT))))
      ENDIF
      IF(Y.LT.ZERO) DAWF=-DAWF
      RETURN
      END
C     ******************************************************************
      FUNCTION DAWERF (Y)
C     *
C     EVALUATE THE DAWSON-ERROR FUNCTION.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     SAVE !OpenMP
      COMMON
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./DERFCM/ C(246),IFIRST(40),ILAST(40),H
      X      = ABS(Y)
      IF(X.LT.10.0D0) THEN
         NX  = INT(X/H)
         NX  = NX+1
         IF  = IFIRST(NX)
         IL  = ILAST(NX)
         T   = C(IL)
         KL  = IL-IF
         DO 10 K=1,KL
         T   = C(IL-K)+X*T
   10    CONTINUE
         DAWERF = T
      ELSE
         TXT = PT5/(X*X)
         TX  = TXT*X
         DAWERF = TX*(ONE+TXT*(ONE+TXT*(THREE+TXT*(15.0D0+TXT*
     1               (105.0D0+945.0D0*TXT)))))
      ENDIF
      RETURN
      END
C     ******************************************************************
      SUBROUTINE RECUR1(L,R1,NIJ,LIJ,ALFI)
C     *
C     RECURSIONS FOR RADIAL INTEGRALS OF TYPE 1.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R1(NIJ,LIJ)
      IF(L.GT.0) GO TO 30
C *** CASE L=0.
      DO 20 LP=3,LIJ
      LP1    = LP-1
      LP2    = LP-2
      FP     = (LP1+LP2)*ALFI
CDIR$ NOVECTOR
C$DIR SCALAR
      DO 10 N=LP,NIJ,2
      R1(N,LP) = R1(N,LP2)-FP*R1(N-1,LP1)
   10 CONTINUE
   20 CONTINUE
      RETURN
C *** CASE L.GT.0.
   30 DO 50 LP=3,LIJ
      LP1    = LP-1
      LP2    = LP-2
      FP     = (LP1+LP2)*ALFI
      NF     = MAX(LP-L,2)
      NF     = NF+MOD(LP+NF,2)
C$DIR SCALAR
      DO 40 N=NF,NIJ,2
      R1(N,LP) = R1(N,LP2)-FP*R1(N-1,LP1)
   40 CONTINUE
   50 CONTINUE
      RETURN
      END
C     ******************************************************************
      SUBROUTINE RECUR2(L,R2,NIJ,LIJC,LIJB,ALFI,BETI)
C     *
C     RECURSIONS FOR RADIAL INTEGRALS OF TYPE 2.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R2(NIJ,LIJC,LIJB)
      LIJCB  = LIJC+LIJB
      IF(LIJCB.LT.4 .OR. NIJ.LT.2) RETURN
      IF(LIJC.GT.2) THEN
         FC  = 3.0D0*ALFI
CDIR$ NOVECTOR
C$DIR SCALAR
         DO 10 N=3,NIJ,2
         R2(N,3,1) = R2(N,1,1)-FC*R2(N-1,2,1)
   10    CONTINUE
      ENDIF
      IF(LIJB.GT.2) THEN
         FB  = 3.0D0*BETI
C$DIR SCALAR
         DO 20 N=3,NIJ,2
         R2(N,1,3) = R2(N,1,1)-FB*R2(N-1,1,2)
   20    CONTINUE
      ENDIF
      LL1    = L+L+1
      DO 60 LCB=5,LIJCB
      NF     = MAX(LCB-LL1,2)
      NF     = NF+MOD(LCB+NF-1,2)
      LCF    = MAX(LCB-LIJB,1)
      LCL    = MIN(LCB-1,LIJC)
      DO 50 LC=LCF,LCL
      LB     = LCB-LC
      IF(LC.GT.LB) THEN
         LC1 = LC-1
         LC2 = LC-2
         FC  = (LC1+LC2)*ALFI
C$DIR SCALAR
         DO 30 N=NF,NIJ,2
         R2(N,LC,LB) = R2(N,LC2,LB)-FC*R2(N-1,LC1,LB)
   30    CONTINUE
      ELSE
         LB1 = LB-1
         LB2 = LB-2
         FB  = (LB1+LB2)*BETI
C$DIR SCALAR
         DO 40 N=NF,NIJ,2
         R2(N,LC,LB) = R2(N,LC,LB2)-FB*R2(N-1,LC,LB1)
   40    CONTINUE
      ENDIF
   50 CONTINUE
   60 CONTINUE
      RETURN
      END
