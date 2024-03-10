      SUBROUTINE SYMTRY (JPRINT)
C     *
C     IMPOSE SYMMETRY CONDITIONS FOR MOLECULAR GEOMETRY.
C     OPTIONALLY PRINT MOLECULAR GEOMETRY (JPRINT.GE.0).
C     *
      USE LIMIT, ONLY: LM1, LMV, LMS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*1 BLSTAR(2),Q(3)
      COMMON
     ./CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./IJWORK/ IARRAY(3,LM1),IJKDUM(9*LM1)
     ./INOPT1/ IN1(300)
     ./NBFILE/ NBF(20)
     ./PARM1 / A(3,LM1),NC(LM1),NB(LM1),NA(LM1),NN(LM1),NATOMS
     ./PARM2 / NSYM,LPAR(LMS),LNUM(LMS),LDEP(LMS)
     ./PARM3 / LOC(LMV),NVAR
     ./PARM7 / DEPFAC
      DATA BLSTAR/' ','*'/
C *** FILE NUMBERS.
      NB6    = NBF(6)
C *** INPUT OPTIONS.
C     SYMTRY MAY BE CALLED FROM READMO/GETGEO, HENCE USE ARRAY IN1.
      IGEOM  = IN1(4)
      IF(NSYM.EQ.0) GO TO 20
C *** LOOP OVER SYMMETRY CONDITIONS.
      DO 10 L=1,NSYM
      I      = LPAR(L)
      J      = LDEP(L)
      K      = LNUM(L)
      IF(K.LE.3) THEN
         A(K,J) = A(K,I)
      ELSE IF(K.EQ.14) THEN
         A(3,J) =-A(3,I)
      ELSE IF(K.LE.13) THEN
         IF(K.EQ. 4) A(3,J)=PT5*PI-A(3,I)
         IF(K.EQ. 5) A(3,J)=PT5*PI+A(3,I)
         IF(K.EQ. 6) A(3,J)=TWO*PI/THREE-A(3,I)
         IF(K.EQ. 7) A(3,J)=TWO*PI/THREE+A(3,I)
         IF(K.EQ. 8) A(3,J)=PI-A(3,I)
         IF(K.EQ. 9) A(3,J)=PI+A(3,I)
         IF(K.EQ.10) A(3,J)=FOUR*PI/THREE-A(3,I)
         IF(K.EQ.11) A(3,J)=FOUR*PI/THREE+A(3,I)
         IF(K.EQ.12) A(3,J)=1.5D0*PI-A(3,I)
         IF(K.EQ.13) A(3,J)=1.5D0*PI+A(3,I)
      ELSE IF(K.LE.20) THEN
         IF(K.EQ.15) A(1,J)=A(1,I)*PT5
         IF(K.EQ.16) A(2,J)=A(2,I)*PT5
         IF(K.EQ.17) A(2,J)=PI-A(2,I)
         IF(K.EQ.18) A(1,J)=A(1,I)*DEPFAC
         IF(K.EQ.19) A(1,J)=A(1,I)*0.763932022D0
         IF(K.EQ.20) A(1,J)=A(1,I)*0.70710678118655D0
      ELSE IF(K.LE.23) THEN
         K      = K-20
         A(K,J) =-A(K,I)
      ELSE IF(K.LE.33) THEN
         IF(K.EQ.24) A(2,J)=A(1,I)
         IF(K.EQ.25) A(1,J)=A(2,I)
         IF(K.EQ.26) A(3,J)=A(1,I)
         IF(K.EQ.27) A(1,J)=A(3,I)
         IF(K.EQ.28) A(3,J)=A(2,I)
         IF(K.EQ.29) A(2,J)=A(3,I)
         IF(K.EQ.30) A(2,J)=PI+A(2,I)
         IF(K.EQ.31) A(3,J)=A(3,I)*PT5
         IF(K.EQ.32) A(3,J)=A(3,I)*TWO
         IF(K.EQ.33) A(3,J)=TWO*PI/THREE-TWO*A(3,I)
      ELSE
         WRITE(NB6,700) K
         STOP 'SYMTRY'
      ENDIF
   10 CONTINUE
C *** DETERMINE POSITION OF OPTIMIZED VARIABLES.
   20 IF(JPRINT.LT.0) RETURN
      DO 30 I=1,NATOMS
      IARRAY(1,I) = 1
      IARRAY(2,I) = 1
      IARRAY(3,I) = 1
   30 CONTINUE
      DO 40 I=1,NVAR
      J      = (LOC(I)+2)/3
      K      = LOC(I)-(J-1)*3
      IARRAY(K,J) = 2
   40 CONTINUE
C *** PRINT THE GEOMETRY IN INTERNAL COORDINATES.
      IF(IGEOM.GT.0) GO TO 60
      WRITE(NB6,600)
      WRITE(NB6,610) NN(1)
      IF(NATOMS.LT.2) GO TO 80
      Q(1)   = BLSTAR(IARRAY(1,2))
      WRITE(NB6,620) NN(2),A(1,2),Q(1),NA(2)
      IF(NATOMS.LT.3) GO TO 80
      Q(1)   = BLSTAR(IARRAY(1,3))
      Q(2)   = BLSTAR(IARRAY(2,3))
      W      = A(2,3)*AFACT
      WRITE(NB6,630) NN(3),A(1,3),Q(1),W,Q(2),NA(3),NB(3)
      IF(NATOMS.LT.4)  GO TO 80
      DO 50 I=4,NATOMS
      Q(1)   = BLSTAR(IARRAY(1,I))
      Q(2)   = BLSTAR(IARRAY(2,I))
      Q(3)   = BLSTAR(IARRAY(3,I))
      W      = A(2,I)*AFACT
      X      = A(3,I)*AFACT
      WRITE(NB6,640) I,NN(I),A(1,I),Q(1),W,Q(2),X,Q(3),NA(I),NB(I),NC(I)
   50 CONTINUE
      GO TO 80
C *** PRINT THE GEOMETRY IN CARTESIAN COORDINATES.
   60 WRITE(NB6,650)
      DO 70 I=1,NATOMS
      WRITE(NB6,640) I,NN(I),(A(J,I),BLSTAR(IARRAY(J,I)),J=1,3)
   70 CONTINUE
C *** CHECK FOR LINK ATOM.
   80 NLINK  = 0
      DO 90 I=1,NATOMS
      IF(NN(I).EQ.86) NLINK=NLINK+1
   90 CONTINUE
      IF(NLINK.EQ.1) THEN
         WRITE(NB6,660)
      ELSE IF(NLINK.GT.1) THEN
         WRITE(NB6,670) NLINK
      ENDIF
      RETURN
  600 FORMAT(/6X,'ATOM',6X,'ATOMIC',15X,'BOND LENGTH',10X,'BOND ANGLE',
     1 11X,'TWIST ANGLE' /5X,'NUMBER',5X,'NUMBER',15X,'(ANGSTROMS)',
     2 11X,'(DEGREES)',12X,'(DEGREES)' /8X, 'I' ,32X,'NA I',16X,
     3 'NB NA I',12X,'NC NB NA I',16X,'NA',5X,'NB',5X,'NC' /)
  610 FORMAT(8X,'1',9X,I2)
  620 FORMAT(8X,'2',9X,I2,F26.5,1X,A1,58X,I2)
  630 FORMAT(8X,'3',9X,I2,7X,2(F19.5,1X,A1),37X,2(I2,5X) )
  640 FORMAT(6X,I3,9X,I2,7X,3(F19.5,1X,A1),11X,3(4X,I3) )
  650 FORMAT(/6X,'ATOM',6X,'ATOMIC',14X,'X-COORDINATE',
     1        9X,'Y-COORDINATE',9X,'Z-COORDINATE'
     2       /5X,'NUMBER',5X,'NUMBER',15X,'(ANGSTROMS)',
     3       10X,'(ANGSTROMS)',10X,'(ANGSTROMS)'/)
  660 FORMAT(//5X,'ELEMENT 86 DENOTES A LINK ATOM FOR QM/MM STUDIES.',
     1       / 5X,'THERE IS ONE SUCH LINK ATOM.')
  670 FORMAT(//5X,'ELEMENT 86 DENOTES A LINK ATOM FOR QM/MM STUDIES.',
     1       / 5X,'NUMBER OF SUCH LINK ATOMS:',I5)
  700 FORMAT(//5X,'SYMMETRY RELATION',I3,' UNDEFINED.'//)
      END