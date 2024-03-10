      SUBROUTINE COMPHF (A,LM5,ICNTR,ESCAL)
C     *
C     COMPUTE ONE-ELECTRON CORE HAMILTONIAN MATRIX AND
C     TWO-ELECTRON CONTRIBUTIONS TO FOCK MATRIX FOR GIVEN DENSITY.
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     A(LM5)    AVAILABLE BUFFER (I,O,S).
C     ICNTR=0   FULL CALCULATION AT REFERENCE GEOMETRY (I).
C     ICNTR>0   PARTIAL CALCULATION FOR GRADIENT AT CENTER ICNTR (I).
C               FOR ICNTR.LE.NUMAT: ALL NONZERO CONTRIBUTIONS.
C               FOR ICNTR.GT.NUMAT: ONLY ONE-ELECTRON CONTRIBUTIONS.
C     ESCAL     CORE REPULSION ENERGY (O).
C     *
C     ARRAYS USED IN THIS ROUTINE (FIRST ADDRESS IN BUFFER).
C     A(LS7)    H(LM4) : CORE HAMILTONIAN MATRIX (O).
C     A(LS6)    FA(LM4): RHF OR UHF-ALPHA FOCK MATRIX (O).
C     A(LU6)    FB(LM4): UHF-BETA FOCK MATRIX (O).
C     A(LS8)    PA(LM4): RHF OR UHF-ALPHA DENSITY MATRIX (I).
C     A(LU8)    PB(LM4): UHF-BETA DENSITY MATRIX (I).
C     A(LS3)    Q(LM6) : SCRATCH ARRAY USED FOR VARIOUS PURPOSES (S).
C     A(LS9)    W(LM9) : TWO-ELECTRON INTEGRALS (I).
C     A(LCD2)   DIELECTRIC MATRIX FOR COSMO SOLVATION TREATMENT (S).
C     A(LC1)    INVERSE A-MATRIX FOR SURFACE CHARGE INTERACTIONS (S).
C     A(LCB1)   B-MATRIX FOR SOLUTE/SURFACE INTERACTIONS (S).
C     A(LCC1)   SCRATCH SPACE (S).
C     A(LCD1)   DIELECTRIC OPERATOR MATRIX (S).
C     A(LC2)    BUFFER 2 FOR SUBSEQUENT COSMO TREATMENT (S).
C     A(LG1)    S(LM4): OVERLAP MATRIX (S).
C     A(LG2)    B(LM4): RESONANCE MATRIX (S).
C     A(LG3)    COR() : LOCAL CORE MATRIX (S).
C     A(LG4)    HO(LM4): ORTHOGONALIZATION CORRECTIONS (S).
C     *
      USE LIMIT, ONLY: LM1, LMX, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./ENERGT/ EE,ENUCLR,EAT,ATHEAT
     ./FLAG4 / INTSUM
     ./INOPT2/ IN2(300)
     ./LIMITS/ LM2,LM3,LM4,LM6,LM7,LM8,LM9
     ./LMCSM1/ LC1,LCB1,LCC1,LCD1,LC2,LCA2,LCB2,LCD2
     ./LMCSM2/ LMC1,LMCA1,LMCB1,LMCC1,LMCD1,LMCD2,LMQ,NPS
     ./LMGRAD/ LG1,LG2,LG3,LG4,LG5,LG6,LG7,LG8
     ./LMSCF / LS1,LS2,LS3,LS4,LS5,LS6,LS7,LS8,LS9
     ./LMUHF / LU1,LU2,LU3,LU4,LU5,LU6,LU7,LU8
     ./NBFILE/ NBF(20)
     ./MMCOM1/ EMM
      DIMENSION A(LM5)
C     *
C *** INPUT OPTIONS.
      IOP    = IN2(2)
      MMINP  = IN2(12)
      NNHCO  = IN2(20)
      ICOSMO = IN2(28)
C     NPRINT = IN2(72)
      MMCOUP = IN2(121)
      IMODE  = IN2(211)
      IF(ICNTR.GT.0) GO TO 10
C *** PREPARATIONS FOR COSMO SOLVATION TREATMENT.
      IF(ICOSMO.GT.0) THEN
         CALL CONSTS (A(LC1),LMC1,NPS)
         CALL DYNCSM
         CALL DMAT   (A(LC1),A(LCB1),A(LCC1),A(LCD1),A(LC2),LMCA1,
     1                LMCB1,LMCC1,LMCD1,LC1-LC2,LCD1-LC2,LMQ,NPS)
      ENDIF
C *** INTEGRAL CALCULATION FOR SCF.
      CALL HCOREP (A(LS7),A(LS9),A(LG1),A(LG2),A(LG3),A(LG4),
     1             LM2,LM4,LM6,LM9)
C *** PEPTIDE BOND CORRECTION IN PM3.
      EMM    = 0.0D0
      IF(NNHCO.GT.0) CALL MMOKF (NNHCO)
C *** CONTRIBUTIONS FROM EXTERNAL POINT CHARGES.
      IF(MMINP.EQ.2 .AND. MMCOUP.GE.2) THEN
         CALL MMINT (A(LS7),A(LS3),LM4,LM6,ENUCLR)
      ENDIF
C *** CONTRIBUTIONS FROM COSMO SOLVATION MODEL.
      IF(ICOSMO.GT.0) THEN
         CALL ADDHCR (A(LS7),A(LCD2),A(LCD2),LM4,LMCD2,LMQ)
         CALL ADDNUC (ENUCLR,A(LCD2),A(LCD2),LMCD2,LMQ)
      ENDIF
C *** CORE REPULSION ENERGY.
      ESCAL = ENUCLR + EMM
C *** TWO-ELECTRON CONTRIBUTIONS TO FOCK MATRIX FOR GIVEN DENSITY.
      IF(ICNTR.LE.NUMAT) THEN
C        PREVENT ARRAY BOUNDARY VIOLATION WITHOUT COSMO.
         IF(LCD2.EQ.0) LCD2 = 1
         CALL FOCK2E (A(LS6),A(LU6),A(LS8),A(LU8),A(LS3),A(LS9),
     1                A(LCD2),LM2,LM4,LM6,LM9,LMCD2,LMQ)
      ENDIF
      RETURN
C
C *** SIMPLIFIED CODE FOR GRADIENT CALCULATIONS WITH DISPLACED ATOM.
   10 CONTINUE
C *** INTEGRAL CALCULATION FOR SCF.
      CALL HCORE2 (A(LS7),A(LS9),A(LG1),A(LG2),A(LG3),A(LG4),
     1             LM2,LM4,LM6,LM9,ICNTR)
C *** PEPTIDE BOND CORRECTION IN PM3.
      EMM    = 0.0D0
      IF(NNHCO.GT.0) CALL MMOKF (NNHCO)
C *** CONTRIBUTIONS FROM EXTERNAL POINT CHARGES.
      IF(MMINP.EQ.2 .AND. MMCOUP.GE.2) THEN
         CALL MMINT2 (A(LS7),A(LS3),LM4,LM6,ENUCLR,NUMAT,ICNTR)
      ENDIF
C *** CORE REPULSION ENERGY.
      ESCAL = ENUCLR + EMM
C *** TWO-ELECTRON CONTRIBUTIONS TO FOCK MATRIX FOR GIVEN DENSITY.
      IF(ICNTR.LE.NUMAT) THEN
         CALL FOCK2G (A(LS6),A(LU6),A(LS8),A(LU8),A(LS3),A(LS9),
     1                LM2,LM4,LM6,LM9,ICNTR)
      ENDIF
      RETURN
      END
C     ******************************************************************
      SUBROUTINE HCORE2 (H,W,S,B,COR,HO,LM2,LM4,LM6,LM9,ICNTR)
C     *
C     CORE HAMILTONIAN CONTRIBUTIONS FOR MNDO AND RELATED METHODS.
C     *
C     ONLY INTEGRALS INVOLVING A SELECTED ATOM ICNTR ARE INCLUDED
C     IN NUMERICAL DIFFERENTIATION FOR CARTESIAN COORDINATE AT ICNTR.
C     *
C     NOTATION. I=INPUT, O=OUTPUT.
C     H(LM4)    CORE HAMILTONIAN MATRIX (O).
C     W(LM9)    TWO-ELECTRON INTEGRALS (O).
C     S(LM4)    OVERLAP MATRIX (O).
C     B(LM4)    RESONANCE MATRIX (O).
C     COR()     LOCAL CORE MATRIX (O).
C     HO(LM4)   ORTHOGONALIZATION CORRECTIONS (O).
C     ICNTR     NUMBER OF SELECTED ATOM (I).
C     *
      USE LIMIT, ONLY: LM1, LMX, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (IZERO=0)
      LOGICAL NOTOMI,OM2,OM3,OM4,ODM2,ODM3
      COMMON
     ./ATOMC / COORD(3,LM1)
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./ENERGT/ EE,ENUCLR,EAT,ATHEAT
     ./EXTRA1/ RI(22),CORE(10,2),REPT(10,2),WW(2025)
     ./EXTRA2/ SIJ(14),T(14),YY(675)
      COMMON
     ./INDEX / INDX(LMX)
     ./INDEXW/ NW(LM1)
     ./INOPT2/ IN2(300)
     ./PARDER/ TORE(LMZ),EHEAT(LMZ),EISOL(LMZ)
     ./PAROPT/ USS(LMZ),UPP(LMZ),ZS(LMZ),ZP(LMZ),BETAS(LMZ),BETAP(LMZ),
     .         ALP(LMZ)
      DIMENSION H(LM4),W(LM9)
      DIMENSION S(LM4),B(LM4),HO(LM4)
      DIMENSION COR(LM2,LM2)
C *** INPUT OPTIONS.
      IOP    = IN2(2)
C *** INITIALIZE SOME VARIABLES.
      NOTOMI = IOP.NE.-5 .AND. IOP.NE.-6 .AND. IOP.NE.-8 .AND. IOP.NE.-9
     1                   .AND. IOP.NE.-22.AND. IOP.NE.-23
      OM2    = IOP.EQ.-6 .AND. NUMAT.GT.2
      OM3    = IOP.EQ.-8 .AND. NUMAT.GT.2
      OM4    = IOP.EQ.-9 
      ODM2   = IOP.EQ.-22 .AND. NUMAT.GT.2
      ODM3   = IOP.EQ.-23 .AND. NUMAT.GT.2
      IJPAIR = INDX(NUMAT)+NUMAT
      KR     = 0
      ENUCLR = ZERO
      H(1:LM4) = ZERO
C     W(1:LM9) = ZERO
      IF(OM2 .OR. OM3 .OR. OM4 .OR. ODM2 .OR. ODM3) THEN
C        NORBS     = NLAST(NUMAT)
C        S(1:LM4)  = ZERO
C        B(1:LM4)  = ZERO
C        HO(1:LM4) = ZERO
C        COR(1:NORBS,1:NUMAT) = ZERO
         GO TO 30
      ENDIF
C
C *** ONE-CENTER TERMS ARE NOT NEEDED FOR NUMERICAL GRADIENT.
C
C *** LOOP OVER ATOM PAIRS FOR OFFDIAGONAL TWO-CENTER TERMS.
C     ATOMS I AND J ARE IDENTIFIED AT THE BEGINNING OF THE LOOP.
C     CODE FOR METHODS WITHOUT THREE-CENTER CONTRIBUTIONS.
      DO 20 I=2,NUMAT
      DO 10 J=1,I-1
      IF(I.NE.ICNTR .AND. J.NE.ICNTR) GO TO 10
      NI     = NAT(I)
      NJ     = NAT(J)
      IA     = NFIRST(I)
      JA     = NFIRST(J)
      IS     = INDX(IA)+IA
      JS     = INDX(JA)+JA
      IORBS  = NLAST(I)-IA+1
      JORBS  = NLAST(J)-JA+1
      ISHELL = I
      JSHELL = J
C *** SPECIAL SECTION FOR H-H PAIR.
      IF(NI.EQ.1 .AND. NJ.EQ.1 .AND. IORBS.EQ.1 .AND. JORBS.EQ.1
     1           .AND. NOTOMI) THEN
         CALL HHPAIR (J,I,COORD,SIJSS,HIJ,WIJ,EN)
         H(INDX(IA)+JA) = HIJ
         H(IS)  = H(IS)-WIJ
         H(JS)  = H(JS)-WIJ
         ENUCLR = ENUCLR+EN
         IJP    = (NW(J)-1)*LM6+NW(I)
         KLP    = (NW(I)-1)*LM6+NW(J)
         W(IJP) = WIJ
         W(KLP) = WIJ
         GO TO 10
      ENDIF
C *** DISTANCE R (AU) AND ROTATION MATRIX.
      CALL ROTMAT (J,I,JORBS,IORBS,NUMAT,COORD,R,YY)
C *** TWO-ELECTRON INTEGRALS.
      IW     = INDX(IORBS)+IORBS
      JW     = INDX(JORBS)+JORBS
C     TWO-ELECTRON INTEGRALS IN LOCAL COORDINATES.
      IF(NOTOMI) THEN
C        COMPUTE AND STORE THE SEMIEMPIRICAL INTEGRALS.
         CALL REPP (NI,NJ,R,RI,CORE)
         IF(IORBS.GE.9 .OR. JORBS.GE.9) THEN
            CALL REPPD (NI,NJ,R,RI,CORE,WW,IW,JW,1)
         ENDIF
      ELSE
C        COMPUTE AND STORE THE ANALYTICAL INTEGRALS.
         CALL REPGAU (R,RI,ISHELL,JSHELL,NI,NJ,IZERO)
         CALL CORDEF (RI,REPT,NI,NJ)
C        SCALE THE ANALYTICAL INTEGRALS.
         CALL GSCALE (R,RI,FKO,NI,NJ)
C        SAVE INTEGRALS TO BE USED FOR CORE-ELECTRON ATTRACTIONS.
         CALL CORDEF (RI,CORE,NI,NJ)
      ENDIF
C     TRANSFORM TWO-ELECTRON INTEGRALS TO MOLECULAR COORDINATES.
      IP     = NW(I)
      JP     = NW(J)
      IF(IORBS.LE.4 .AND. JORBS.LE.4) THEN
         CALL ROTATE (IW,JW,IP,JP,KR,RI,YY,W,W,LM6,LM9,0)
      ELSE
         CALL ROTD  (WW,YY,IW,JW)
         CALL W2MAT (IP,JP,WW,W,LM6,IW,JW)
      ENDIF
      CALL W2COPY (IW,JW,IP,JP,W,LM6)
C *** RESONANCE INTEGRALS.
      IF(NOTOMI) THEN
         CALL BETAIJ (NI,NJ,R,SIJ,T)
         CALL ROTBET (IA,JA,IORBS,JORBS,T,YY,H,LM4)
      ELSE
         CALL BETOM  (IZERO,ISHELL,JSHELL,NI,NJ,R,SIJ,T)
         IF(OM2 .OR. OM3 .OR. OM4 .OR. ODM2 .OR. ODM3) THEN
            CALL ROTBET (IA,JA,IORBS,JORBS,SIJ,YY,S,LM4)
            CALL ROTBET (IA,JA,IORBS,JORBS,T,YY,B,LM4)
         ELSE
            CALL ROTBET (IA,JA,IORBS,JORBS,T,YY,H,LM4)
         ENDIF
      ENDIF
C *** CORE-ELECTRON ATTRACTIONS.
      IF(NOTOMI) THEN
         CALL ROTCOH (IA,JA,IORBS,JORBS,IS,JS,CORE,YY,H,LM4)
      ELSE
         CALL COROM  (ISHELL,JSHELL,NI,NJ,R,FKO,CORE,REPT,SIJ,T)
         CALL ROTCOH (IA,JA,IORBS,JORBS,IS,JS,CORE,YY,H,LM4)
      ENDIF
C *** CORE-CORE REPULSIONS.
      IF(NOTOMI) THEN
         WIJ = RI(1)
         CALL CORREP (I,J,NI,NJ,R,WIJ,EN)
      ELSE
         EN  = FKO*TORE(NI)*TORE(NJ)*EV/R
      ENDIF
      ENUCLR = ENUCLR+EN
   10 CONTINUE
   20 CONTINUE
      RETURN
C
C *** LOOP OVER ATOM PAIRS FOR OFFDIAGONAL TWO-CENTER TERMS.
C     ATOMS I AND J ARE IDENTIFIED AT THE BEGINNING OF THE LOOP.
C     CODE FOR METHODS WITH THREE-CENTER CONTRIBUTIONS (OM2,OM3,ODM2,ODM3).
   30 CONTINUE
      DO 50 I=2,NUMAT 
      DO 40 J=1,I-1
      NI     = NAT(I)
      NJ     = NAT(J)
      IA     = NFIRST(I)
      JA     = NFIRST(J)
      IS     = INDX(IA)+IA
      JS     = INDX(JA)+JA
      IORBS  = NLAST(I)-IA+1
      JORBS  = NLAST(J)-JA+1
      ISHELL = I
      JSHELL = J
C *** DISTANCE R (AU) AND ROTATION MATRIX.
      CALL ROTMAT (J,I,JORBS,IORBS,NUMAT,COORD,R,YY)
C *** TWO-ELECTRON INTEGRALS.
      IW     = INDX(IORBS)+IORBS
      JW     = INDX(JORBS)+JORBS
C     TWO-ELECTRON INTEGRALS IN LOCAL COORDINATES.
C     COMPUTE AND STORE THE ANALYTICAL INTEGRALS.
      IF(I.EQ.ICNTR .OR. J.EQ.ICNTR .OR. OM2 .OR. ODM2) THEN
         MODGAU = 0
         IF(I.NE.ICNTR .AND. J.NE.ICNTR) MODGAU=1
         CALL REPGAU (R,RI,ISHELL,JSHELL,NI,NJ,MODGAU)
         CALL CORDEF (RI,REPT,NI,NJ)
C        SCALE THE ANALYTICAL INTEGRALS.
         CALL GSCALE (R,RI,FKO,NI,NJ)
C        SAVE INTEGRALS TO BE USED FOR CORE-ELECTRON ATTRACTIONS.
         CALL CORDEF (RI,CORE,NI,NJ)
      ENDIF 
C     TRANSFORM TWO-ELECTRON INTEGRALS TO MOLECULAR COORDINATES.
      IF(I.EQ.ICNTR .OR. J.EQ.ICNTR) THEN
         IP  = NW(I)
         JP  = NW(J)
         IF(IORBS.LE.4 .AND. JORBS.LE.4) THEN
            CALL ROTATE (IW,JW,IP,JP,KR,RI,YY,W,W,LM6,LM9,0)
         ELSE
            CALL ROTD  (WW,YY,IW,JW)
            CALL W2MAT (IP,JP,WW,W,LM6,IW,JW)
         ENDIF
         CALL W2COPY (IW,JW,IP,JP,W,LM6)
      ENDIF
C *** RESONANCE INTEGRALS.
      CALL BETOM  (IZERO,ISHELL,JSHELL,NI,NJ,R,SIJ,T)
      CALL ROTBET (IA,JA,IORBS,JORBS,SIJ,YY,S,LM4)
      CALL ROTBET (IA,JA,IORBS,JORBS,T,YY,B,LM4)
C *** CORE-ELECTRON ATTRACTIONS.
C     SAVE LOCAL CORE INTEGRALS FOR USE IN BETORT (OM2 OR ODM2).
      IF(OM2 .OR. ODM2) THEN
         COR(IA,J) = USS(NI)+CORE(1,1)
         COR(JA,I) = USS(NJ)+CORE(1,2)
         IF(IORBS.GE.4) THEN
            TEMP = UPP(NI)+(CORE(3,1)+TWO*CORE(4,1))/THREE
            COR(IA+1,J) = TEMP
            COR(IA+2,J) = TEMP
            COR(IA+3,J) = TEMP
         ENDIF
         IF(JORBS.GE.4) THEN
            TEMP = UPP(NJ)+(CORE(3,2)+TWO*CORE(4,2))/THREE
            COR(JA+1,I) = TEMP
            COR(JA+2,I) = TEMP
            COR(JA+3,I) = TEMP
         ENDIF
      ENDIF
C     UPDATE CORE-ELECTRON ATTRACTIONS.
      IF(I.EQ.ICNTR .OR. J.EQ.ICNTR) THEN
         CALL COROM  (ISHELL,JSHELL,NI,NJ,R,FKO,CORE,REPT,SIJ,T)
         CALL ROTCOH (IA,JA,IORBS,JORBS,IS,JS,CORE,YY,H,LM4)
      ENDIF
C *** CORE-CORE REPULSIONS.
      IF(I.EQ.ICNTR .OR. J.EQ.ICNTR) THEN
         EN     = FKO*TORE(NI)*TORE(NJ)*EV/R
         ENUCLR = ENUCLR+EN
      ENDIF
   40 CONTINUE
   50 CONTINUE
C
C *** ORTHOGONALIZATION CORRECTIONS FOR RESONANCE INTEGRALS.
      IF(ICNTR.LE.NUMAT) THEN
         MODE = -ICNTR
         IF(OM2 .OR. ODM2) THEN
            H(1:LM4) = H(1:LM4)+B(1:LM4)
            CALL BETORT (HO,S,B,COR,B,LM2,LM4,MODE)
         ELSE IF(OM3 .OR. ODM3) THEN
            H(1:LM4) = H(1:LM4)+B(1:LM4)
            CALL BETOR3 (HO,S,B,B,LM2,LM4,MODE)
         ELSE IF(OM4) THEN
            CALL BETLOW (HO,S,B,LM2,LM4)
         ENDIF
         H(1:LM4) = H(1:LM4)+HO(1:LM4)
      ENDIF
      RETURN
      END
C     ******************************************************************
      SUBROUTINE W2COPY (IW,JW,IP,JP,W,LM6)
C     *
C     COPY TWO-CENTER TWO-ELECTRON INTEGRALS INSIDE A SQUARE MATRIX.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION W(LM6,LM6)
      DO 20 I=IP,IP+IW-1
      DO 10 J=JP,JP+JW-1
      W(J,I) = W(I,J)
   10 CONTINUE
   20 CONTINUE
      RETURN
      END
C     ******************************************************************
      SUBROUTINE MMINT2 (H,HM,LM4,LM6,ENUCLR,NUMAT,ICNTR)
C     *
C     ELECTROSTATIC CONTRIBUTIONS FROM EXTERNAL POINT CHARGES
C     TO THE CORE HAMILTONIAN AND THE CORE-CORE REPULSIONS.
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     H(LM4)    CORE HAMILTONIAN MATRIX (I,O).
C     HM(LM6)   CONTRIBUTIONS FROM GIVEN POINT CHARGE (S).
C     ENUCLR    CORE-CORE REPULSION ENERGY (I,O).
C     NUMAT     NUMBER OF ATOMS (I).
C     ICNTR     NUMBER OF SELECTED CENTER (I).
C     *
C     THE SEQUENTIAL VERSION DOES NOT NEED THE SCRATCH ARRAY HM(LM6).
C     IT CAN COLLECT THE CONTRIBUTIONS DIRECTLY IN H(LM4) BY USING: 
C     CALL QMINTC (H,LM4,ENUCQM,COORDM(1,M),PTCHG,2)
C     CALL QMINT  (H,LM4,ENUCQM,COORDM(1,M),PTCHG,2)
C     *
      USE LIMIT, ONLY: LM1M, LMI
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./FINDX1/ IP(LMI),IP1(LMI),IP2(LMI)
     ./INOPT2/ IN2(300)
     ./QMMM1 / COORDM(3,LM1M),CHARGM(LM1M)
      DIMENSION H(LM4),HM(LM6)
C *** INPUT OPTIONS.
      NUMATM  = IN2(120)
      MMPOT   = IN2(122)
C *** LOOP OVER EXTERNAL POINT CHARGES.
C.OMP PARALLEL DO SCHEDULE(STATIC)
C.OMP.PRIVATE(ENUCQM,HM,PTCHG)
      DO 20 M=1,NUMATM
      IF(ICNTR.GT.NUMAT .AND. ICNTR.NE.(M+NUMAT)) GO TO 20
      PTCHG  = CHARGM(M)
      IF(PTCHG.NE.ZERO) THEN
         HM(1:LM6) = ZERO
         IF(MMPOT.LE.2) THEN
            CALL QMINTC (HM,LM6,ENUCQM,COORDM(1,M),PTCHG,1,ICNTR)
         ELSE
            CALL QMINT  (HM,LM6,ENUCQM,COORDM(1,M),PTCHG,1,ICNTR)
         ENDIF
C.OMP CRITICAL
         ENUCLR = ENUCLR+ENUCQM
         H(IP(1:LM6)) = H(IP(1:LM6)) + HM(1:LM6)
C.OMP END CRITICAL
      ENDIF
   20 CONTINUE
C.OMP END PARALLEL DO
      RETURN
      END
C     ******************************************************************
      SUBROUTINE FOCK2E (FA,FB,PA,PB,Q,W,D2,LM2,LM4,LM6,LM9,LMD2,LMQ)
C     *
C     COMPUTE TWO-ELECTRON PART OF FOCK MATRIX FOR GIVEN DENSITY.
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     FA(LM4)   RHF OR UHF-ALPHA FOCK MATRIX (O).
C     FB(LM4)   UHF-BETA FOCK MATRIX (O).
C     PA(LM4)   RHF OR UHF-ALPHA DENSITY MATRIX (I).
C     PB(LM4)   UHF-BETA DENSITY MATRIX (I).
C     Q(LM6)    SCRATCH ARRAY USED FOR VARIOUS PURPOSES (S).
C     W(LM9)    TWO-ELECTRON INTEGRALS (I).
C     D2(LMD2)  DIELECTRIC MATRIX FOR COSMO SOLVATION TREATMENT (I).
C     LMQ       TOTAL NUMBER OF SOLUTE CHARGES IN COSMO (I).
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LFOCK
      LOGICAL UHF
      COMMON
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./INOPT2/ IN2(300)
     ./NBFILE/ NBF(20)
     ./ORBITS/ NUMB,NORBS,NMOS,NALPHA,NBETA
     ./UHF   / UHF
C     ARRAYS IN THE ARGUMENT LIST.
      DIMENSION FA(LM4),FB(LM4)
      DIMENSION PA(LM4),PB(LM4),Q(LM2*8),W(LM9),D2(LMD2)
C *** FILE NUMBERS.
      NB1    = NBF(1)
C *** INPUT OPTIONS.
      ICOSMO = IN2(28)
      IMODE  = IN2(211)
      INOUT  = IN2(212)
C *** INITIALIZATION.
      LFOCK  = IMODE.LT.-1 .AND. LM6.GE.(-IMODE)
C *** FILE HANDLING.
      IF(INOUT.GT.0) THEN
         REWIND NB1
         READ(NB1) PA
         IF(UHF) READ(NB1) PB
      ENDIF
C *** CONSTRUCT THE F-MATRIX FOR RHF OR THE FALPHA-MATRIX FOR UHF.
      FA(1:LM4) = ZERO
      IF(LFOCK) THEN
         CALL FOCK  (FA,PA,PB,Q,W,LM4,LM6)
      ELSE
         CALL FOCKX (FA,PA,PB,Q,W,LM4,LM6)
      ENDIF
C     INCLUDE TERMS FROM COSMO SOLVATION MODEL.
      IF(ICOSMO.GT.0) CALL ADDFCK (FA,PA,PB,Q,D2,D2,LM4,LMD2,LMQ)
C *** CONSTRUCT THE FBETA-MATRIX FOR UHF.
      IF(UHF) THEN
         FB(1:LM4) = ZERO
         IF(LFOCK) THEN
            CALL FOCK  (FB,PB,PA,Q,W,LM4,LM6)
         ELSE
            CALL FOCKX (FB,PB,PA,Q,W,LM4,LM6)
         ENDIF
C        INCLUDE TERMS FROM COSMO SOLVATION MODEL.
         IF(ICOSMO.GT.0) CALL ADDFCK (FB,PA,PB,Q,D2,D2,LM4,LMD2,LMQ)
      ENDIF
      RETURN
      END
C     ******************************************************************
      SUBROUTINE FOCK2G (FA,FB,PA,PB,Q,W,LM2,LM4,LM6,LM9,ICNTR)
C     *
C     COMPUTE TWO-ELECTRON PART OF FOCK MATRIX FOR GIVEN DENSITY.
C     INCLUDE ONLY CONTRIBUTIONS INVOLVING ATOM ICNTR.
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     FA(LM4)   RHF OR UHF-ALPHA FOCK MATRIX (O).
C     FB(LM4)   UHF-BETA FOCK MATRIX (O).
C     PA(LM4)   RHF OR UHF-ALPHA DENSITY MATRIX (I).
C     PB(LM4)   UHF-BETA DENSITY MATRIX (I).
C     Q(LM6)    SCRATCH ARRAY USED FOR VARIOUS PURPOSES (S).
C     W(LM9)    TWO-ELECTRON INTEGRALS (I).
C     ICNTR     NUMBER OF SELECTED ATOM (I).
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL UHF
      COMMON
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./INOPT2/ IN2(300)
     ./NBFILE/ NBF(20)
     ./ORBITS/ NUMB,NORBS,NMOS,NALPHA,NBETA
     ./UHF   / UHF
C     ARRAYS IN THE ARGUMENT LIST.
      DIMENSION FA(LM4),FB(LM4)
      DIMENSION PA(LM4),PB(LM4),Q(LM6),W(LM9)
C *** FILE NUMBERS.
      NB1    = NBF(1)
C *** INPUT OPTIONS.
      INOUT  = IN2(212)
C *** FILE HANDLING.
      IF(INOUT.GT.0) THEN
         REWIND NB1
         READ(NB1) PA
         IF(UHF) READ(NB1) PB
      ENDIF
C *** CONSTRUCT THE F-MATRIX FOR RHF OR THE FALPHA-MATRIX FOR UHF.
      FA(1:LM4) = ZERO
      CALL FOCK2X (FA,PA,PB,Q,W,LM4,LM6,ICNTR)
C *** CONSTRUCT THE FBETA-MATRIX FOR UHF.
      IF(UHF) THEN
         FB(1:LM4) = ZERO
         CALL FOCK2X (FB,PB,PA,Q,W,LM4,LM6,ICNTR)
      ENDIF
      RETURN
      END
C     ******************************************************************
      SUBROUTINE FOCK2X (F,PA,PB,Q,W,LM4,LM6,ICNTR)
C     *
C     TWO-ELECTRON CONTRIBUTIONS TO MNDO-TYPE FOCK MATRIX.
C     SCALAR CODE FOR TWO-CENTER EXCHANGE CONTRIBUTIONS.
C     INCLUDE ONLY TWO-CENTER CONTRIBUTIONS INVOLVING ATOM ICNTR.
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     F(LM4)    FOCK MATRIX (I,O).
C     PA(LM4)   RHF OR UHF-ALPHA DENSITY MATRIX (I).
C     PB(LM4)   UHF-BETA DENSITY MATRIX (I).
C     Q(LM6)    SCRATCH ARRAY FOR ONE-CENTER PAIR TERMS (S).
C     W(LM6,*)  TWO-ELECTRON INTEGRALS (I).
C     ICNTR     NUMBER OF SELECTED ATOM (I).
C     *
      USE LIMIT, ONLY: LM1, LMI, LMX
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL UHF
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./FINDX1/ IP(LMI),IP1(LMI),IP2(LMI)
     ./INDEX / INDX(LMX)
     ./INDEXW/ NW(LM1)
     ./UHF   / UHF
      DIMENSION F(LM4),PA(LM4),PB(LM4),Q(LM6),W(LM6,LM6)
      DIMENSION IWW(4,4)
      DATA IWW/1,2,4,7,2,3,5,8,4,5,6,9,7,8,9,10/
C
C *** TWO-CENTER COULOMB CONTRIBUTIONS.
C     PREPARE VECTOR OF ONE-CENTER DENSITY MATRIX ELEMENTS.
      DO 10 KL=1,LM6
      Q(KL)  = (PA(IP(KL))+PB(IP(KL)))
      IF(IP1(KL).NE.IP2(KL)) Q(KL)=Q(KL)*TWO
   10 CONTINUE
C     PROCESS NONZERO TWO-CENTER COULOMB CONTRIBUTIONS.
C     ALTERNATIVE CODE: REPLACE "CREP" AND ADD "CADD".
C     VERSION WITH DO 16 JJ=1,NUMAT SLIGHTLY FASTER.
      DO 21 II=1,NUMAT
      IORBS  = NLAST(II)-NFIRST(II)+1
      IW     = INDX(IORBS)+IORBS
      IJMIN  = NW(II)-1
      DO 20 IJ=IJMIN+1,IJMIN+IW
      SUMIJ  = ZERO
CREP  DO 16 JJ=1,II-1
      DO 16 JJ=1,NUMAT
CREP  IF(II.NE.ICNTR .AND. JJ.NE.ICNTR) GO TO 16
      IF((II.NE.ICNTR .AND. JJ.NE.ICNTR) .OR. II.EQ.JJ) GO TO 16
      JORBS  = NLAST(JJ)-NFIRST(JJ)+1
      JW     = INDX(JORBS)+JORBS
      KLMIN  = NW(JJ)-1
      DO 15 KL=KLMIN+1,KLMIN+JW
      SUMIJ     = SUMIJ+Q(KL)*W(KL,IJ)
CADD  F(IP(KL)) = F(IP(KL))+Q(IJ)*W(KL,IJ)
   15 CONTINUE
   16 CONTINUE
      F(IP(IJ)) = F(IP(IJ))+SUMIJ
   20 CONTINUE
   21 CONTINUE
C
C *** TWO-CENTER EXCHANGE CONTRIBUTIONS.
      DO 55 II=1,NUMAT
      IA     = NFIRST(II)
      IB     = NLAST(II)
      IC     = INDX(IA)
      IORBS  = IB-IA+1
      IW     = INDX(IORBS)+IORBS
      IJ     = NW(II)-1
      DO 50 JJ=1,II-1
      IF(II.NE.ICNTR .AND. JJ.NE.ICNTR) GO TO 50
      JA     = NFIRST(JJ)
      JB     = NLAST(JJ)
      JORBS  = JB-JA+1
      JW     = INDX(JORBS)+JORBS
      KL     = NW(JJ)-1
      IF(IW.EQ.1 .AND. JW.EQ.1) THEN
         IS  = IC+JA
         F(IS) = F(IS)-PA(IS)*W(IJ+1,KL+1)
      ELSE IF(IW.EQ.1 .AND. JW.EQ.10) THEN
         IS  = IC+JA
         IX  = IS+1
         IY  = IS+2
         IZ  = IS+3
         IJS = IJ+1
         F(IS) = F(IS)-PA(IS)*W(KL+1,IJS)-PA(IX)*W(KL+2,IJS)
     1                -PA(IY)*W(KL+4,IJS)-PA(IZ)*W(KL+7,IJS)
         F(IX) = F(IX)-PA(IS)*W(KL+2,IJS)-PA(IX)*W(KL+3,IJS)
     1                -PA(IY)*W(KL+5,IJS)-PA(IZ)*W(KL+8,IJS)
         F(IY) = F(IY)-PA(IS)*W(KL+4,IJS)-PA(IX)*W(KL+5,IJS)
     1                -PA(IY)*W(KL+6,IJS)-PA(IZ)*W(KL+9,IJS)
         F(IZ) = F(IZ)-PA(IS)*W(KL+7,IJS)-PA(IX)*W(KL+8,IJS)
     1                -PA(IY)*W(KL+9,IJS)-PA(IZ)*W(KL+10,IJS)
      ELSE IF(IW.EQ.10 .AND. JW.EQ.1) THEN
         IS  = IC+JA
         IX  = INDX(IA+1)+JA
         IY  = INDX(IA+2)+JA
         IZ  = INDX(IA+3)+JA
         KLS = KL+1
         F(IS) = F(IS)-PA(IS)*W(IJ+1,KLS)-PA(IX)*W(IJ+2,KLS)
     1                -PA(IY)*W(IJ+4,KLS)-PA(IZ)*W(IJ+7,KLS)
         F(IX) = F(IX)-PA(IS)*W(IJ+2,KLS)-PA(IX)*W(IJ+3,KLS)
     1                -PA(IY)*W(IJ+5,KLS)-PA(IZ)*W(IJ+8,KLS)
         F(IY) = F(IY)-PA(IS)*W(IJ+4,KLS)-PA(IX)*W(IJ+5,KLS)
     1                -PA(IY)*W(IJ+6,KLS)-PA(IZ)*W(IJ+9,KLS)
         F(IZ) = F(IZ)-PA(IS)*W(IJ+7,KLS)-PA(IX)*W(IJ+8,KLS)
     1                -PA(IY)*W(IJ+9,KLS)-PA(IZ)*W(IJ+10,KLS)
      ELSE IF(IW.EQ.10 .AND. JW.EQ.10) THEN
         DO 40 I=1,4
         IS  = INDX(IA+I-1)+JA
         IX  = IS+1
         IY  = IS+2
         IZ  = IS+3
         DO 30 K=1,4
         KS  = INDX(IA+K-1)+JA
         KX  = KS+1
         KY  = KS+2
         KZ  = KS+3
         IJK = IJ+IWW(I,K)
         F(IS) = F(IS)-PA(KS)*W(KL+1,IJK)-PA(KX)*W(KL+2,IJK)
     1                -PA(KY)*W(KL+4,IJK)-PA(KZ)*W(KL+7,IJK)
         F(IX) = F(IX)-PA(KS)*W(KL+2,IJK)-PA(KX)*W(KL+3,IJK)
     1                -PA(KY)*W(KL+5,IJK)-PA(KZ)*W(KL+8,IJK)
         F(IY) = F(IY)-PA(KS)*W(KL+4,IJK)-PA(KX)*W(KL+5,IJK)
     1                -PA(KY)*W(KL+6,IJK)-PA(KZ)*W(KL+9,IJK)
         F(IZ) = F(IZ)-PA(KS)*W(KL+7,IJK)-PA(KX)*W(KL+8,IJK)
     1                -PA(KY)*W(KL+9,IJK)-PA(KZ)*W(KL+10,IJK)
   30    CONTINUE
   40    CONTINUE
      ELSE
C     GENERAL CODE   - ALSO VALID FOR D ORBITALS.
C     FIRST  SECTION - CONTRIBUTIONS FROM (II,KL).
C     SECOND SECTION - CONTRIBUTIONS FROM (IJ,KL) WITH I.NE.J.
         DO 43 I=IA,IB
         KA     = INDX(I)
         IJW    = IJ+INDX(I-IA+2)
         KLW    = KL
         DO 42 K=JA,JB
         IK     = KA+K
         DO 41 L=JA,K
         IL     = KA+L
         KLW    = KLW+1
         A      = W(KLW,IJW)
         F(IK) =  F(IK) - A*PA(IL)
         IF(K.NE.L) F(IL) =  F(IL) - A*PA(IK)
   41    CONTINUE
   42    CONTINUE
   43    CONTINUE
         DO 47 I=IA+1,IB
         KA     = INDX(I)
         DO 46 J=IA,I-1
         KB     = INDX(J)
         IJW    = IJ+INDX(I-IA+1)+J-IA+1
         KLW    = KL
         DO 45 K=JA,JB
         IK     = KA+K
         JK     = KB+K
         DO 44 L=JA,K
         IL     = KA+L
         JL     = KB+L
         KLW    = KLW+1
         A      = W(KLW,IJW)
         F(IK)  = F(IK) - A*PA(JL)
         F(JK)  = F(JK) - A*PA(IL)
         IF(K.NE.L) THEN
            F(IL)  = F(IL) - A*PA(JK)
            F(JL)  = F(JL) - A*PA(IK)
         ENDIF
   44    CONTINUE
   45    CONTINUE
   46    CONTINUE
   47    CONTINUE
      ENDIF
   50 CONTINUE
   55 CONTINUE
      RETURN
      END
