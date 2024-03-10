      SUBROUTINE GRAPHM (C,E,LM2,LM3,MODE,NB13,NPRINT)
C     *
C     SAVE SCF EIGENVECTORS FOR POSTPROCESSING VIA MOLDEN.
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     C(LM2,*)  EIGENVECTORS (I).
C     E(LM3)    EIGENVALUES (I).
C     MODE      TYPE OF EIGENVECTORS (I).
C               = 0 RHF
C               = 1 UHF-ALPHA.
C               = 2 UHF-BETA.
C     NB13      FILE NUMBER (I).
C     NPRINT    PRINTING FLAG (I).
C     *
      USE LIMIT, ONLY: LM1, LMX
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*3 NGROUP
      CHARACTER*4 IRREP
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./DELEMT/ NELMD,NFOCK
     ./IJWORK/ NSYM(LMX),ISYM(3,LM1)
C    ./NBFILE/ NBF(20)
     ./OCCNM / OCCA(LMX),OCCB(LMX)
     ./ORBITS/ NUMB,NORBS,NMOS,NALPHA,NBETA
     ./SYMLAB/ NGROUP(7),IRREP(26)
     ./SYMLIM/ IRREPA(7),IRREPN(7)
      DIMENSION C(LM2,LM3),E(LM2)
C *** FILE NUMBERS.
C     NB6    = NBF(6)
C     NB13   = NBF(13)
C *** ASSIGN MO SYMMETRY.
      CALL MOSYM (C,E,LM2,LM3,LM3,NSYM,ISYM,ISUB,NPRINT)
      IF(ISUB.GT.0) NLAB = IRREPA(ISUB)
C *** SAVE MOS ON FILE NB13 USING MOLDEN CONVENTIONS.
      DO 30 I=1,NORBS
      IF(ISUB.GT.0) THEN
          WRITE(NB13,'(A,A)') ' Sym= ', IRREP(NLAB+NSYM(I))
      ELSE
          WRITE(NB13,'(A,A)') ' Sym= ', 'A'
      ENDIF
      WRITE(NB13,'(A,F20.10)') ' Ene= ', E(I)
      IF(MODE.EQ.0) THEN
         WRITE(NB13,'(A)') ' Spin= Alpha'
         WRITE(NB13,'(A,F5.1)') ' Occup= ', OCCA(I)*2.0D0
      ELSE IF(MODE.EQ.1) THEN
         WRITE(NB13,'(A)') ' Spin= Alpha'
         WRITE(NB13,'(A,F5.1)') ' Occup= ', OCCA(I)
      ELSE
         WRITE(NB13,'(A)') ' Spin= Beta'
         WRITE(NB13,'(A,F5.1)') ' Occup= ', OCCB(I)
      ENDIF
C     COPY EIGENVECTORS IF THERE ARE NO D ORBITALS.
      IF(NELMD.LE.0) THEN
         DO 10 J=1,NORBS
         WRITE(NB13,'(I5,5X,F20.10)') J, C(J,I)
   10    CONTINUE
      ELSE
C     IN CASE OF D ORBITALS, COPY COEFFICIENTS OF S AND P AOS,
C     AND ADAPT COEFFICIENTS OF D AOS TO MOLDEN CONVENTIONS.
         KA = 1
         DO 20 K=1,NUMAT
         JA     = NFIRST(K)
         KORBS  = NLAST(K)-NFIRST(K)+1
         WRITE(NB13,'(I5,5X,F20.10)') KA, C(JA,I)
         IF(KORBS.GE.4) THEN
            WRITE(NB13,'(I5,5X,F20.10)') KA+1, C(JA+1,I)
            WRITE(NB13,'(I5,5X,F20.10)') KA+2, C(JA+2,I)
            WRITE(NB13,'(I5,5X,F20.10)') KA+3, C(JA+3,I)
         ENDIF
         IF(KORBS.GE.9) THEN
            SQ3 = ONE / SQRT(THREE)
            CXX = PT5*C(JA+4,I) - PT5*SQ3*C(JA+6,I)
            CYY =-PT5*C(JA+4,I) - PT5*SQ3*C(JA+6,I)
            CZZ = SQ3*C(JA+6,I)
            WRITE(NB13,'(I5,5X,F20.10)') KA+4, CXX
            WRITE(NB13,'(I5,5X,F20.10)') KA+5, CYY
            WRITE(NB13,'(I5,5X,F20.10)') KA+6, CZZ
            WRITE(NB13,'(I5,5X,F20.10)') KA+7, C(JA+7,I)
            WRITE(NB13,'(I5,5X,F20.10)') KA+8, C(JA+8,I)
            WRITE(NB13,'(I5,5X,F20.10)') KA+9, C(JA+5,I)
            KA = KA+1
         ENDIF
         KA = KA+KORBS
   20    CONTINUE
      ENDIF
   30 CONTINUE
      RETURN
      END
