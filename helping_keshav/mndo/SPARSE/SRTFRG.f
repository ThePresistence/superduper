      SUBROUTINE SRTFRG (PSA,JPSA,IPSA,QSA,JQSA,IQSA,NPSA,
     1                   NUMATS,NORBSS,NFS,NLS,MAP,INVMAP,JPRINT)
C     *
C     SORT SPARSE DENSITY MATRIX (PSA->QSA) BY REORDERING COLUMNS.
C     *
C     NOTATION. I=INPUT, O=OUTPUT.
C     PSA       ORIGINAL SPARSE DENSITY MATRIX OF MOLECULE (I).
C     JPSA      POINTERS FOR SPARSE MATRIX IN CSR FORMAT (I).
C     IPSA      COLUMN INDICES FOR SPARSE MATRIX IN CSR FORMAT (I).
C     PSA       ORIGINAL SPARSE DENSITY MATRIX OF MOLECULE (O).
C     JPSA      POINTERS FOR SPARSE MATRIX IN CSR FORMAT (O).
C     IPSA      COLUMN INDICES FOR SPARSE MATRIX IN CSR FORMAT (O).
C     NPSA      DIMENSION OF PSA,JPSA,QSA,JQSA (I).
C     NORBSS+1  DIMENSION OF IPSA,IQSA (I).
C     NUMATS    NUMBER OF ATOMS IN MOLECULE (I).
C     NORBSS    NUMBER OF ORBITALS IN MOLECULE (I).
C     NFS       NUMBER OF FIRST ORBITAL PER ATOM IN MOLECULE (I).
C     NLS       NUMBER OF LAST  ORBITAL PER ATOM IN MOLECULE (I).
C     MAP       FORWARD MAPPING FOR ATOMS IN MOLECULE AND FRAGMENTS (I).
C     INVMAP    INVERSE MAPPING FOR ATOMS IN MOLECULE AND FRAGMENTS (I).
C     JPRINT    PRINTING FLAG (I).
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION, POINTER :: PSA(:)
      INTEGER, POINTER :: IPSA(:),JPSA(:)
      DIMENSION QSA(NPSA),JQSA(NPSA),IQSA(NORBSS+1)
      DIMENSION NFS(NUMATS),NLS(NUMATS)
      DIMENSION MAP(NUMATS),INVMAP(NUMATS)
C     AUTOMATIC SCRATCH ARRAYS.
      DIMENSION JFS(NUMATS),JLS(NUMATS)
C *** DEFINE FIRST AND LAST ORBITAL PER ATOM FOR FRAGMENT LIST.
      KK     = 0
      DO 10 J=1,NUMATS
      KK     = KK+1
      JFS(J) = KK
      KK     = KK+NLS(INVMAP(J))-NFS(INVMAP(J))
      JLS(J) = KK
   10 CONTINUE
C *** OUTER LOOP OVER ATOMS IN MOLECULE (ORIGINAL LIST).
      I2     = 0
      DO 30 II=1,NUMATS
      JJ     = MAP(II)
      IA     = NFS(II)
      IB     = NLS(II)
      JA     = JFS(JJ)
      JB     = JLS(JJ)
      DO 20 I=IA,IB
      J      = JA+I-IA
      J1     = IPSA(J)
      J2     = IPSA(J+1)-1
      I1     = I2+1
      I2     = I1+J2-J1
      QSA(I1:I2)  = PSA(J1:J2)
      JQSA(I1:I2) = JPSA(J1:J2)
      IQSA(I)     = I1
   20 CONTINUE
   30 CONTINUE
      IQSA(NORBSS+1) = NPSA+1
C *** DEBUG PRINT.
      IF(JPRINT.GE.9) THEN
         WRITE(6,500)
         WRITE(6,510) NFS
         WRITE(6,511) NLS
         WRITE(6,512) JFS
         WRITE(6,513) JLS
         WRITE(6,514) MAP
         WRITE(6,515) INVMAP
         WRITE(6,520)
         WRITE(6,530) IPSA
         WRITE(6,540) IQSA
         WRITE(6,550) JPSA
         WRITE(6,560) JQSA
         WRITE(6,570) PSA
         WRITE(6,580) QSA
      ENDIF
      RETURN
  500 FORMAT(///1X,'SRTFRG: DEBUG PRINT'/)
  510 FORMAT(1X,'NFS :',10I5)
  511 FORMAT(1X,'NLS :',10I5)
  512 FORMAT(1X,'JFS :',10I5)
  513 FORMAT(1X,'JLS :',10I5)
  514 FORMAT(1X,'MAP :',10I5)
  515 FORMAT(1X,'INV :',10I5)
  520 FORMAT(1X)
  530 FORMAT(1X,'IPSA:',10I5)
  540 FORMAT(1X,'IQSA:',10I5)
  550 FORMAT(1X,'JPSA:',10I5)
  560 FORMAT(1X,'JQSA:',10I5)
  570 FORMAT(1X,'PSA :',10F8.4)
  580 FORMAT(1X,'QSA :',10F8.4)
      END