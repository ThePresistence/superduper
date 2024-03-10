      SUBROUTINE GUESSP (CA,CB,PA,PB,LM2,LM3,LM4,MPDIAG)
C     *
C     DEFINITION OF INITIAL DENSITY MATRIX.
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     CA(LM2,*) RHF OR UHF-ALPHA MO EIGENVECTORS (O,S).
C     CB(LM2,*) UHF-BETA MO EIGENVECTORS (O,S).
C     PA(LM4)   RHF OR UHF-ALPHA DENSITY MATRIX (O).
C     PB(LM4)   UHF-BETA DENSITY MATRIX (O).
C     MPDIAG    FLAG TO INDICATE TYPE OF DENSITY MATRIX (O).
C               = 0  SYMMETRIC, WITH NONZERO OFF-DIAGONAL ELEMENTS.
C               = 1  DIAGONAL, NONZERO ELEMENTS DIFFERENT.
C               = 2  DIAGONAL, NONZERO ELEMENTS ALL EQUAL.
C               = 3  UNIT MATRIX (MULTIPLIED BY 0.5).
C     *
      USE LIMIT, ONLY: LM1, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (TINY=1.0D-10)
      LOGICAL UHF
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./DNBND / III(LMZ),IIID(LMZ)
     ./HALFE / IODD,JODD
     ./INOPT2/ IN2(300)
     ./NBFILE/ NBF(20)
     ./ORBITS/ NUMB,NORBS,NMOS,NALPHA,NBETA
     ./PARDER/ CORE(LMZ),EHEAT(LMZ),EISOL(LMZ)
     ./UHF   / UHF
      DIMENSION CA(LM2,LM3),CB(LM2,LM3),PA(LM4),PB(LM4)
C *** CHECK WHETHER INITIAL DENSITY MATRIX IS ALREADY AVAILABLE.
C     IN THE CASE OF A JOB CONTINUATION RUN (IN2(37)=MIDDLE.GT.0)
C     THE PROGRAM DEFINES IN2(67)=KTRIAL=-1 (SEE SUBROUTINE INPUT).
C     IN THIS CASE THE INITIAL DENSITY MATRIX IS READ VIA DENSAV.
      IF(IN2(67).LT.0) RETURN
C *** FILE NUMBERS.
C     NB5    = NBF(5)
      NB6    = NBF(6)
C *** INPUT OPTIONS.
      KHARGE = IN2(65)
      KTRIAL = IN2(67)
      NPRINT = IN2(72)
      INOUT  = IN2(212)
C *** READ DENSITY MATRIX FROM FILE NB11.
      MPDIAG = 0
      IF(KTRIAL.EQ.11) THEN
         NBIN = NBF(11)
         REWIND NBIN
         READ(NBIN,ERR=10,END=10) PA
         IF(UHF) READ(NBIN,ERR=10,END=10) PB
         GO TO 300
      ENDIF
C     READ EIGENVECTORS FROM FILE NB12.
      IF(KTRIAL.EQ.12) THEN
         NBIN = NBF(12)
         REWIND NBIN
         READ(NBIN,ERR=10,END=10) CA
         IF(UHF) READ(NBIN,ERR=10,END=10) CB
         GO TO 200
      ENDIF
C     READ RHF DENSITY MATRIX FOR UHF GUESS FROM FILE NB11.
      IF(KTRIAL.EQ.13) THEN
         NBIN = NBF(11)
         REWIND NBIN
         READ(NBIN,ERR=10,END=10) PA
         GO TO 100
      ENDIF
C     ERROR MESSAGE FOR FAULTY INPUT FROM FILE.
   10 IF(KTRIAL.GE.11 .AND. KTRIAL.LE.13) WRITE(NB6,400) NBIN
C *** INITIALIZATION.
      DO 20 I=1,LM4
      PA(I)  = ZERO
   20 CONTINUE
C *** COMPUTE SIMPLIFIED DIAGONAL TRIAL DENSITY MATRIX (KTRIAL=1).
C     THE ELECTRONS ARE DISTRIBUTED EVENLY OVER ALL ORBITALS.
      IF(KTRIAL.EQ.1) THEN
         XEL    = DBLE(NALPHA+NBETA)/DBLE(NORBS)
         XEL    = XEL*PT5
         DO 30 I=1,NORBS
         II     = (I*(I+1))/2
         PA(II) = XEL
   30    CONTINUE
         MPDIAG = 2
         IF((NALPHA+NBETA).EQ.NORBS) MPDIAG=3
         GO TO 100
      ENDIF
C *** COMPUTE SIMPLIFIED DIAGONAL TRIAL DENSITY MATRIX (KTRIAL=2).
C     THE ELECTRONS ARE DISTRIBUTED EVENLY OVER THE ORBITALS IN
C     EACH ATOM. ANY OVERALL CHARGE IS DISTRIBUTED EVENLY.
      IF(KTRIAL.EQ.2) THEN
         YY     = DBLE(KHARGE)/DBLE(NORBS)
         DO 50 I=1,NUMAT
         IA     = NFIRST(I)
         IB     = NLAST(I)
         NI     = NAT(I)
         DO 40 J=IA,IB
         II     = (J*(J+1))/2
         PA(II) = (CORE(NI)-YY)*PT5
   40    CONTINUE
   50    CONTINUE
         MPDIAG = 1
         GO TO 100
      ENDIF
C *** COMPUTE STANDARD DIAGONAL TRIAL DENSITY MATRIX (KTRIAL=0).
C     NSPORB IS THE NUMBER OF ATOMIC ORBITALS INITIALLY POPULATED.
C     D ORBITALS OF MAIN-GROUP ELEMENTS ARE NOT POPULATED.
C     P ORBITALS OF TRANSITION ELEMENTS ARE NOT POPULATED.
C     POLARIZATION FUNCTIONS IN AN EXTENDED BASIS ARE NOT POPULATED.
      NSPORB = NORBS
      DO 70 I=1,NUMAT
      NI     = NAT(I)
      IORBS  = NLAST(I)-NFIRST(I)+1
      IF(IORBS.EQ.9) THEN
         IF(III(NI).LE.IIID(NI)) NSPORB=NSPORB-5
         IF(III(NI).GT.IIID(NI)) NSPORB=NSPORB-3
      ENDIF
      IF(IORBS.EQ.4 .AND. NI.LE.2) NSPORB=NSPORB-3
   70 CONTINUE
      YY     = DBLE(KHARGE)/DBLE(NSPORB)
C     LOOP OVER ALL ATOMS.
      DO 80 I=1,NUMAT
      IA     = NFIRST(I)
      IORBS  = NLAST(I)-NFIRST(I)+1
      IS     = IA*(IA+1)/2
      NI     = NAT(I)
C     ATOMS WITH AN S BASIS (POSSIBLY WITH P POLARIZATION FUNCTIONS).
      IF(IORBS.EQ.1 .OR. (IORBS.EQ.4 .AND. NI.LE.2)) THEN
         PA(IS) = (CORE(NI)-YY)*PT5
C     ATOMS WITH AN SP-BASIS.
      ELSE IF(IORBS.EQ.4) THEN
         W   = (CORE(NI)*PT25-YY)*PT5
         PA(IS)        = W
         PA(IS+IA+1)   = W
         PA(IS+2*IA+3) = W
         PA(IS+3*IA+6) = W
C     MAIN-GROUP ELEMENTS WITH AN SPD-BASIS.
      ELSE IF(IORBS.EQ.9 .AND. (III(NI).LE.IIID(NI))) THEN
         W   = (CORE(NI)*PT25-YY)*PT5
         PA(IS)        = W
         PA(IS+IA+1)   = W
         PA(IS+2*IA+3) = W
         PA(IS+3*IA+6) = W
C     TRANSITION-METAL ELEMENTS WITH AN SPD-BASIS.
      ELSE IF(IORBS.EQ.9 .AND. (III(NI).GT.IIID(NI))) THEN
         TEMP = CORE(NI)-YY*6.0D0
C        UP TO 10 ELECTRONS, PUT INTO S AND D ORBITALS.
         IF(TEMP.LT.10.0D0) THEN
            W   = TEMP/12.0D0
            PA(IS)         = W
            PA(IS+4*IA+10) = W
            PA(IS+5*IA+15) = W
            PA(IS+6*IA+21) = W
            PA(IS+7*IA+28) = W
            PA(IS+8*IA+36) = W
C        MORE THAN 10 ELECTRONS, FILL D ORBITALS, PUT REST IN S ORBITAL.
         ELSE
            W   = (TEMP-10.0D0)*PT5
            PA(IS)         = W
            PA(IS+4*IA+10) = ONE
            PA(IS+5*IA+15) = ONE
            PA(IS+6*IA+21) = ONE
            PA(IS+7*IA+28) = ONE
            PA(IS+8*IA+36) = ONE
         ENDIF
      ENDIF
   80 CONTINUE
C *** ASSIGN TYPE OF DIAGONAL DENSITY MATRIX GENERATED.
      MPDIAG = 3
      DO 90 I=2,NORBS
      II = (I*(I+1))/2
      IF(ABS(PA(II)-PA(1)).GT.TINY) THEN
         MPDIAG = 1
         GO TO 100
      ENDIF
      IF(ABS(PA(II)-PT5).GT.TINY) THEN
         MPDIAG = 2
      ENDIF
   90 CONTINUE
C *** DETERMINE UHF TRIAL DENSITY MATRIX (KTRIAL=0,1,13).
C     PERTURB INITIAL DENSITY MATRICES FOR UHF SINGLETS.
  100 CONTINUE
      IF(UHF) THEN
         DA  = NALPHA
         DB  = NBETA
         FA  = TWO*DA/(DA+DB)
         FB  = TWO*DB/(DA+DB)
         DO 110 I=1,LM4
         PB(I) = PA(I)*FB
         PA(I) = PA(I)*FA
  110    CONTINUE
         MPDIAG = MIN(MPDIAG,2)
         IF(NALPHA.EQ.NBETA) THEN
            DA = 0.98D0
            DB = TWO-DA
            K  = 0
            DO 120 J=1,NORBS
            K  = K+J
            DC = DA
            DA = DB
            DB = DC
            PA(K) = PA(K)*DA
            PB(K) = PB(K)*DB
  120       CONTINUE
            MPDIAG = MIN(MPDIAG,1)
         ENDIF
      ENDIF
      GO TO 300
C *** INITIAL DENSITY MATRIX FROM EIGENVECTORS ON FILE NB12.
  200 KK     = 0
      DO 230 I=1,NORBS
      DO 220 J=1,I
      KK     = KK+1
      SUM    = ZERO
      DO 210 K=1,NUMB
      OCCNM  = ONE
      IF(K.EQ.IODD .OR. K.EQ.JODD) OCCNM=PT5
      SUM    = SUM+OCCNM*CA(I,K)*CB(J,K)
  210 CONTINUE
      PA(KK) = SUM
  220 CONTINUE
  230 CONTINUE
      IF(UHF) THEN
         KK  = 0
         DO 260 I=1,NORBS
         DO 250 J=1,I
         KK  = KK+1
         SUM = ZERO
         DO 240 K=1,NBETA
         SUM = SUM+CB(I,K)*CB(J,K)
  240    CONTINUE
         PB(KK) = SUM
  250    CONTINUE
  260    CONTINUE
      ENDIF
C *** PRINTING SECTION.
  300 IF(NPRINT.LT.0) GO TO 310
      IF(NPRINT.EQ.5 .AND. (KTRIAL.EQ.0 .OR. NORBS.GT.100)) THEN
         IF(UHF) THEN
            WRITE(NB6,610)
            WRITE(NB6,630) (PA(I*(I+1)/2),I=1,NORBS)
            WRITE(NB6,620)
            WRITE(NB6,630) (PB(I*(I+1)/2),I=1,NORBS)
         ELSE
            WRITE(NB6,600)
            WRITE(NB6,630) (PA(I*(I+1)/2)*2.0D0,I=1,NORBS)
         ENDIF
         WRITE(NB6,640) MPDIAG
      ELSE IF(NPRINT.GE.5) THEN
         IF(UHF) THEN
            WRITE(NB6,510)
            CALL VECPRT (PA,LM4,NORBS)
            WRITE(NB6,520)
            CALL VECPRT (PB,LM4,NORBS)
         ELSE
            WRITE(NB6,500)
            CALL VECPRT (PA,LM4,NORBS)
         ENDIF
         WRITE(NB6,640) MPDIAG
      ENDIF
C *** SAVE DENSITY MATRIX ON FILE NB1.
  310 IF(INOUT.GT.0) THEN
         NB1 = NBF(1)
         REWIND NB1
         WRITE(NB1) PA
         IF(UHF) WRITE(NB1) PB
      ENDIF
C *** SET FLAG THAT INITIAL DENSITY MATRIX HAS BEEN GENERATED.
      IN2(67) = -1
      RETURN
  400 FORMAT (///1X,'ERROR OR END-OF-FILE WHEN READING FILE',I3,
     1           1X,'IN SUBROUTINE GUESSP.',
     2        /  1X,'A STANDARD INITIAL DENSITY MATRIX IS USED.')
  500 FORMAT (///1X,'INITIAL DENSITY MATRIX (MULTIPLIED BY 0.5).'/)
  510 FORMAT (///1X,'INITIAL ALPHA DENSITY MATRIX.'/)
  520 FORMAT (///1X,'INITIAL BETA DENSITY MATRIX.'/)
  600 FORMAT (///1X,'INITIAL DIAGONAL DENSITY MATRIX ELEMENTS.'/)
  610 FORMAT (///1X,'INITIAL DIAGONAL ALPHA DENSITY MATRIX ELEMENTS.'/)
  620 FORMAT (///1X,'INITIAL DIAGONAL BETA DENSITY MATRIX ELEMENTS.'/)
  630 FORMAT (   1X,10F10.5)
  640 FORMAT (// 1X,'TYPE OF INITIAL DENSITY MATRIX: MPDIAG =',I2)
      END
