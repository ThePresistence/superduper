      SUBROUTINE PATHL (A,LM5,SCFCAL)
C     *
C     REACTION PATH CALCULATION.
C     LINEAR INTERPOLATION BETWEEN TWO INPUT GEOMETRIES.
C     *
      USE LIMIT, ONLY: LM1, LMV, LMR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL SCFCAL
      CHARACTER*80  KTITLE,KOMENT
      COMMON
     ./ATOML / AA1(3,LM1),AA2(3,LM1)
     ./DIPOL / DIP(7)
     ./ERG   / ENERGY,GRAD(LMV),GNORM,CNORM
     ./FLAG1 / KTITLE,KOMENT
     ./INOPT2/ IN2(300)
     ./NBFILE/ NBF(20)
     ./PARM1 / AA(3,LM1),NA(LM1,4),NATOMS
     ./PARM4 / RC(LMR),LREACT,LTOTAL
     ./PARM6 / LXYZ(3)
      DIMENSION A(LM5)
      DIMENSION HEATS(LMR),DIPS(LMR)
      DIMENSION CNORMS(LMR),GNORMS(LMR)
C *** FILE NUMBERS.
      NB6    = NBF(6)
      NB7    = NBF(7)
C *** INPUT OPTIONS.
      JOP    = IN2(3)
      NSAV7  = IN2(14)
      MAXEND = IN2(32)
C *** LOOP OVER ALL POINTS ON THE REACTION PATH.
      DO 20 L=1,LTOTAL
C     DETERMINE CURRENT GEOMETRY.
      DO 10 I=1,NATOMS
      DO 10 J=1,3
      AA(J,I)= (1.0D0-RC(L))*AA1(J,I) + RC(L)*LXYZ(J)*AA2(J,I)
   10 CONTINUE
      WRITE(NB6,500) KOMENT,KTITLE
      WRITE(NB6,510) L
      CALL SYMTRY(+1)
      CALL GMETRY(-1)
      CALL GEOSAV (NATOMS,NSAV7,NB7)
C *** DO THE CALCULATION FOR A GIVEN POINT.
C     SINGLE SCF CALCULATION.
      IF(JOP.EQ.-1 .OR. MAXEND.EQ.1) THEN
         ICALL = 0
         CALL SCF (A,LM5,ICALL,SCFCAL)
         IF(ICALL.EQ.-1) THEN
            LMAX = L-1
            GO TO 30
         ENDIF
C     SINGLE SCF AND GRADIENT CALCULATION.
      ELSE IF(JOP.EQ.-2) THEN
         ICALL = 1
         CALL SCF (A,LM5,ICALL,SCFCAL)
         IF(ICALL.EQ.-1) THEN
            LMAX = L-1
            GO TO 30
         ENDIF
      ENDIF
C *** STORE RESULTS.
      HEATS (L) = ENERGY
      DIPS  (L) = DIP(1)
      CNORMS(L) = CNORM
      GNORMS(L) = GNORM
   20 CONTINUE
C *** PRINTING SECTION.
   30 CONTINUE
      WRITE(NB6,520)
      DO 40 L=1,LTOTAL
      WRITE(NB6,530) L,HEATS(L),DIPS(L),GNORMS(L),CNORMS(L)
   40 CONTINUE
      RETURN
  500 FORMAT(///1X,A,/1X,A)
  510 FORMAT(/  1X,'POINT',I3,' ON INTERPOLATED PATHWAY',//)
  520 FORMAT(///1X,'SUMMARY OF RESULTS FROM REACTION PATH CALCULATION.',
     1       // 3X,'I      HEAT OF    DIPOLE   OPTIMIZED    CARTESIAN',
     2       /  3X,'      FORMATION   MOMENT    GRADIENT    GRADIENT',
     3       /  3X,'     (KCAL/MOL)  (DEBYE)      NORM        NORM'/)
  530 FORMAT(   1X,I3,F11.3,F13.5,F9.3,F11.3,F12.3)
      END
