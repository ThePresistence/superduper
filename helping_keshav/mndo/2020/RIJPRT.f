      SUBROUTINE RIJPRT (COORD,NN,NATOMS,KSTOP)
C     *
C     PRINT INTERATOMIC DISTANCES AND CHECK FOR ZERO DISTANCES.
C     *
C     THE FORMATS ARE THE SAME AS IN SUBROUTINE VECPRT.
C     THE DISTANCES ARE COMPUTED AS NEEDED FROM THE CURRENT CARTESIAN
C     COORDINATES SO THAT THERE IS NO LARGE ARRAY FOR DISTANCES.
C     *
C     NOTATION. I=INPUT, O=OUTPUT.
C     COORD     CARTESIAN COORDINATES (I).
C     NN        ATOMIC NUMBERS (I).
C     NATOMS    NUMBER OF ATOMS (I).
C     KSTOP     NUMBER OF ATOM PAIRS WITH ZERO DISTANCE (O).
C               DUMMY ATOMS ARE IGNORED IN THIS COUNT.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (SMALL=1.0D-06)
      CHARACTER*4 LINE(31)
      COMMON /NBFILE/ NBF(20)
      DIMENSION COORD(3,NATOMS),NN(NATOMS)
      DIMENSION RIJ(10)
      DATA LINE / 31*'----' /
C *** FILE NUMBERS.
      NB6    = NBF(6)
C *** LOOP OVER OUTPUT BLOCKS.
      KSTOP  = 0
      NPASS  = 1+(NATOMS-1)/10
      DO 40 N=1,NPASS
      JA     = 1+(N-1)*10
      JCOL   = MIN((NATOMS+1-JA),10)
      JB     = JA-1+JCOL
      WRITE(NB6,500) (J,J=JA,JB)
      WRITE(NB6,510) (LINE(K),K=1,3*JCOL+1)
C *** LOOP OVER COLUMNS OF ONE BLOCK.
      DO 30 I=JA,NATOMS
      JB     = JA-1+MIN(JCOL,I-JA+1)
      DO 10 J=JA,JB
      RIJ(J-JA+1) = SQRT( (COORD(1,I)-COORD(1,J))**2
     1                   +(COORD(2,I)-COORD(2,J))**2
     2                   +(COORD(3,I)-COORD(3,J))**2 )
   10 CONTINUE
      WRITE(NB6,520) I,(RIJ(J),J=1,JB-JA+1)
C *** CHECK FOR ZERO DISTANCES.
      IF(NN(I).NE.99) THEN
         DO 20 J=JA,JB
         IF(NN(J).NE.99 .AND. I.NE.J .AND. RIJ(J-JA+1).LT.SMALL) THEN
            KSTOP = KSTOP+1
         ENDIF
   20    CONTINUE
      ENDIF
   30 CONTINUE
   40 CONTINUE
      RETURN
  500 FORMAT(// 2X,10(7X,I5))
  510 FORMAT(   1X,31A4)
  520 FORMAT(   1X,I3,1X,10F12.5)
      END
