      SUBROUTINE MOOUT(C,EIG,ISUB,NSYM,NR,NC,LM2,LM3,NPRINT,ICASE)
C     *
C     PRINT ORBITAL ENERGIES EIG(NC) AND MO COEFFICIENTS C(NR,NC).
C     PRINT OCCUPATION NUMBERS IN THE CASE IMOCC.GT.1.
C     PRINT SYMMETRY LABELS IN THE CASE ISUB.GT.0.
C     *
      USE LIMIT, ONLY: LM1, LMX
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*2 ELEMNT
      CHARACTER*3 NGROUP
      CHARACTER*4 IRREP
      CHARACTER*5 ORB(9)
      LOGICAL UHF
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./ELEMTS/ ELEMNT(107)
     ./NBFILE/ NBF(20)
     ./OCCFL / IMOCC,NOCCA,NOCCB,MSUB,MOSUMA,MOSUMB,MOCCA(8),MOCCB(8)
     ./OCCNM / OCCA(LMX),OCCB(LMX)
     ./PARM1 / AA(3,LM1),NABC(LM1*3),NN(LM1),NATOMS
     ./SYMLAB/ NGROUP(7),IRREP(26)
     ./SYMLIM/ IRREPA(7),IRREPN(7)
     ./UHF   / UHF
      DIMENSION C(LM2,LM3),EIG(LM3),NSYM(LM3)
      DATA ORB/'  s  ','  x  ','  y  ','  z  ','x2-y2',' xz  ',
     1         ' z2  ',' yz  ',' xy  '/
C *** FILE NUMBERS.
      NB6    = NBF(6)
C *** CHECK PRINTING FLAG.
      IF(NPRINT.LT.-1) RETURN
      IF(NPRINT.LE. 2 .AND. ISUB.GT.0) WRITE(NB6,300) NGROUP(ISUB)
      IF(NPRINT.LT. 0) GO TO 100
C *** STANDARD PRINTOUT.
      IF(.NOT.UHF) THEN
         WRITE(NB6,400)
      ELSE IF(ICASE.EQ.1) THEN
         WRITE(NB6,410)
      ELSE
         WRITE(NB6,420)
      ENDIF
C     NUMBER OF BLOCKS (JBLOCK EIGENVALUES AND EIGENVECTORS PER BLOCK).
      JBLOCK = 9
      IF(NPRINT.GT.0) JBLOCK=5
      KBLOCK = 1+(NC-1)/JBLOCK
      IF(ISUB.GT.0) NLAB = IRREPA(ISUB)
C     LOOP OVER BLOCKS.
      DO 50 K=1,KBLOCK
      KA     = 1+(K-1)*JBLOCK
      KB     = MIN(NC,KA+JBLOCK-1)
      WRITE(NB6,500) (J,J=KA,KB)
      WRITE(NB6,510) (EIG(J),J=KA,KB)
      IF(IMOCC.GT.1) THEN
         IF(.NOT.UHF) THEN
            WRITE(NB6,510) (OCCA(J)*TWO,J=KA,KB)
         ELSE IF(ICASE.EQ.1) THEN
            WRITE(NB6,510) (OCCA(J),J=KA,KB)
         ELSE
            WRITE(NB6,510) (OCCB(J),J=KA,KB)
         ENDIF
      ENDIF
      IF(ISUB.GT.0) WRITE(NB6,520) (IRREP(NLAB+NSYM(J)),J=KA,KB)
      WRITE(NB6,530)
C     LOOP OVER ATOMS AND ATOMIC ORBITALS.
      IIDUM  = 0
      IORB   = 0
      DO 40 II=1,NATOMS
      NI     = NN(II)
      IF(NI.EQ.99) THEN
         IIDUM = IIDUM+1
         GO TO 40
      ENDIF
      IAT   = II-IIDUM
      IA    = NFIRST(IAT)
      IB    = NLAST(IAT)
      IF(NATOMS.LT.10) THEN
         DO 10 I=IA,IB
         IORB  = IORB+1
         WRITE(NB6,540) IORB,ELEMNT(NI),II,ORB(I-IA+1),(C(I,J),J=KA,KB)
   10    CONTINUE
      ELSE IF(NATOMS.LT.100) THEN
         DO 20 I=IA,IB
         IORB  = IORB+1
         WRITE(NB6,550) IORB,ELEMNT(NI),II,ORB(I-IA+1),(C(I,J),J=KA,KB)
   20    CONTINUE
      ELSE
         DO 30 I=IA,IB
         IORB  = IORB+1
         WRITE(NB6,560) IORB,ELEMNT(NI),II,ORB(I-IA+1),(C(I,J),J=KA,KB)
   30    CONTINUE
      ENDIF
      WRITE(NB6,530)
   40 CONTINUE
   50 CONTINUE
      RETURN
C *** MINIMUM PRINTOUT.
  100 CONTINUE
      IF(.NOT.UHF) THEN
         WRITE(NB6,430)
      ELSE IF(ICASE.EQ.1) THEN
         WRITE(NB6,440)
      ELSE
         WRITE(NB6,450)
      ENDIF
C     NUMBER OF LINES (JBLOCK EIGENVALUES PER LINE).
      JBLOCK = 5
      KBLOCK = 1+(NC-1)/JBLOCK
      IF(ISUB.GT.0) NLAB = IRREPA(ISUB)
      WRITE(NB6,600) (J,J=1,JBLOCK)
C     LOOP OVER LINES.
      DO 150 K=1,KBLOCK
      KK     = (K-1)*JBLOCK
      KA     = 1+KK
      KB     = MIN(NC,KA+JBLOCK-1)
      WRITE(NB6,610) KK,(EIG(J),J=KA,KB)
      IF(IMOCC.GT.1) THEN
         IF(.NOT.UHF) THEN
            WRITE(NB6,610) KK,(OCCA(J)*TWO,J=KA,KB)
         ELSE IF(ICASE.EQ.1) THEN
            WRITE(NB6,610) KK,(OCCA(J),J=KA,KB)
         ELSE
            WRITE(NB6,610) KK,(OCCB(J),J=KA,KB)
         ENDIF
      ENDIF
      IF(ISUB.GT.0) WRITE(NB6,620) (IRREP(NLAB+NSYM(J)),J=KA,KB)
  150 CONTINUE
      RETURN
  300 FORMAT(// 5X,'POINT GROUP ',A3,' ASSIGNED.')
  400 FORMAT(///5X,'EIGENVALUES AND EIGENVECTORS.')
  410 FORMAT(///5X,'EIGENVALUES AND EIGENVECTORS FOR ALPHA SPIN.')
  420 FORMAT(///5X,'EIGENVALUES AND EIGENVECTORS FOR BETA SPIN.')
  430 FORMAT(///5X,'EIGENVALUES (EV).')
  440 FORMAT(///5X,'EIGENVALUES (EV) FOR ALPHA SPIN.')
  450 FORMAT(///5X,'EIGENVALUES (EV) FOR BETA SPIN.')
  500 FORMAT(//21X,I5,8I12/)
  510 FORMAT(17X,9F12.5)
  520 FORMAT(16X,9(8X,A4))
  530 FORMAT(1X)
  540 FORMAT(1X,I2,2X,A2,'(',I1,')',1X,A5,F13.5,8F12.5)
  550 FORMAT(1X,I3,2X,A2,'(',I2,')',1X,A5,F11.5,8F12.5)
  560 FORMAT(1X,I4,2X,A2,'(',I3,')',1X,A5,F9.5,8F12.5)
  600 FORMAT(//9X,I5,9I12)
  610 FORMAT(1X,I4,10F12.5)
  620 FORMAT(4X,10(8X,A4))
      END
