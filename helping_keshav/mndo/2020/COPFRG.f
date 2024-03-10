      SUBROUTINE COPFRG (IFRAG,NUMATS,NFS,NFRAGS,P,LM2F,PA,LM4)
C     *
C     COPY FRAGMENT DENSITY INTO MOLECULAR DENSITY ARRAY.
C     *
C     NOTATION. I=INPUT, O=OUTPUT.
C     IFRAG     NUMBER OF CURRENT FRAGMENT (I).
C     NUMATS    NUMBER OF ATOMS IN MOLECULE (I).
C     NFS       FIRST ORBITAL OF ATOMS IN MOLECULE (I).
C     NFRAGS    FRAGMENT ASSIGNMENT OF ATOMS IN MOLECULE (I).
C     P(LM2F,*) SQUARE DENSITY MATRIX OF FRAGMENT (I).
C     PA(LM4)   INITIAL DENSITY MATRIX OF MOLECULE (O).
C     *
      USE LIMIT, ONLY: LM1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
      DIMENSION NFS(NUMATS),NFRAGS(NUMATS)
      DIMENSION P(LM2F,LM2F)
      DIMENSION PA(LM4)
C     AUTOMATIC SCRATCH ARRAY.
      DIMENSION MAP(NUMAT)
C *** ASSIGN ATOMS IN CURRENT FRAGMENT TO ATOMS IN MOLECULE.
      II     = 0
      DO 10 K=1,NUMATS
      IF(IFRAG.EQ.NFRAGS(K)) THEN
         II  = II+1
         MAP(II) = K
      ENDIF
   10 CONTINUE
C *** OUTER LOOP OVER ATOMS IN FRAGMENT.
      DO 50 II=1,NUMAT
      IA     = NFIRST(II)
      IB     = NLAST(II)
      IAS    = NFS(MAP(II))
      DO 40 I=IA,IB
      IS     = IAS+I-IA
      ISBAS  = (IS*(IS-1))/2
C *** INNER LOOP OVER ATOMS IN FRAGMENT.
      DO 30 JJ=1,II
      JA     = NFIRST(JJ)
      JB     = NLAST(JJ)
      JAS    = NFS(MAP(JJ))
      JMAX   = MIN(JB,I)
      DO 20 J=JA,JMAX
      PA(ISBAS+JAS+J-JA) = P(I,J)
   20 CONTINUE
   30 CONTINUE
   40 CONTINUE
   50 CONTINUE
      RETURN
      END
