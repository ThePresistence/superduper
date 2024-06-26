      SUBROUTINE PRTGUG (CIVEC,EGEN,FGEN,IEGEN,IFGEN,
     1                   LMEGEN,LMFGEN,LMVEC,ICALL)
C     *
C     DEBUG PRINT OF DATA GENERATED IN SUBROUTINE GUGACI.
C     *
C     NOTATION. I=INPUT.
C     CIVEC(LMVEC)    CI VECTOR (I).
C     EGEN(LMEGEN)    ONE-ELECTRON GENERATORS (I).
C     FGEN(LMFGEN)    TWO-ELECTRON GENERATORS (I).
C     IEGEN(LMEGEN)   PACKED INDICES FOR ONE-ELECTRON GENERATORS (I).
C     IFGEN(LMFGEN)   PACKED INDICES FOR TWO-ELECTRON GENERATORS (I).
C     ICALL           CONTROLS THE NUMBER OF DATA PRINTED (I).
C                     = 1 PRINT ONLY DIMENSIONS OF ARRAYS.
C                     = 2 PRINT ALSO ONE-ELECTRON AND TWO-ELECTRON
C                         GENERATORS WITH THEIR ASSOCIATE INDICES.
C                     = 3 PRINT ALSO CI VECTOR.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./INOPT2/ IN2(300)
     ./NBFILE/ NBF(20)
      DIMENSION CIVEC(LMVEC),EGEN(LMEGEN),FGEN(LMFGEN)
      DIMENSION IEGEN(LMEGEN),IFGEN(LMFGEN)
C *** FILE NUMBERS.
      NB6    = NBF(6)
C *** INITIALIZATION.
      NCIO   = IN2(131)+IN2(132)
      NCIGAM = (NCIO*(NCIO+1)*(NCIO**2+NCIO+2))/8
C *** PRINT NUMBER OF CONFIGURATIONS AND GENERATORS.
      WRITE(NB6,500) LMVEC,LMEGEN,LMFGEN,NCIGAM
      IF(ICALL.LE.1) RETURN
C *** PRINT ONE-ELECTRON GENERATORS.
      WRITE(NB6,520)
      DO 10 I=1,LMEGEN-1
      CALL GUGAFU (IEGEN(I),NCIGAM,L,M,IP,IQ,IR,IS)
      WRITE(NB6,530) I,EGEN(I),IEGEN(I),L,M,IP,IQ,IR,IS
   10 CONTINUE
C *** PRINT TWO-ELECTRON GENERATORS.
      WRITE(NB6,540)
      DO 20 I=1,LMFGEN-1
      CALL GUGAFU (IFGEN(I),NCIGAM,L,M,IP,IQ,IR,IS)
      WRITE(NB6,530) I,FGEN(I),IFGEN(I),L,M,IP,IQ,IR,IS
   20 CONTINUE
      IF(ICALL.LE.2) RETURN
C *** PRINT CI VECTOR.
      WRITE(NB6,510) (CIVEC(I),I=1,LMVEC)
      RETURN
  500 FORMAT(///1X,'NUMBER OF GUGA-CI CONFIGURATIONS   :',I12,
     1       /  1X,'NUMBER OF ONE-ELECTRON GENERATORS  :',I12,
     2       /  1X,'NUMBER OF TWO-ELECTRON GENERATORS  :',I12,
     3       /  1X,'SIZE OF TWO-PARTICLE DENSITY MATRIX:',I12)
  510 FORMAT(///1X,'GUGA-CI VECTOR.',/(1X,10F10.5))
  520 FORMAT(// 1X,'ONE-ELECTRON GENERATORS.'
     1       //10X,'I     VALUE',
     2         15X,'INDEX    MU    NU    I   J   K   L'/)
  530 FORMAT(   1X,I10,F10.5,I20,2I6,1X,4I4)
  540 FORMAT(///1X,'TWO-ELECTRON GENERATORS.'
     1       //10X,'I     VALUE',
     2         15X,'INDEX    MU    NU    I   J   K   L'/)
      END
