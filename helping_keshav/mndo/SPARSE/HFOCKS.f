      SUBROUTINE HFOCKS (F,FS,IFS,JFS,NFS,N,IA,IB,MODE,CUTF)
C     *
C     TRANSFER NONZERO MATRIX ELEMENTS BETWEEN LOCAL ARRAY F(N,*)
C     AND SPARSE STORAGE IN CSR FORMAT (FS,IFS,JFS).
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     F(N,*)    PARTIAL LOCAL MATRIX (O/I).
C     FS(*)     SPARSE GLOBAL MATRIX (I/O).
C     IFS(N+1)  POINTER ARRAY FOR SPARSE MATRIX IN CSR FORMAT (I/O).
C     JFS(*)    COLUMN INDICES FOR SPARSE MATRIX IN CSR FORMAT (I/O).
C     NFS       TOTAL NUMBER OF NONZERO ELEMENTS (O/I,O).
C     N         DIMENSION OF LOCAL MATRIX (I).
C     IA        LOWER INDEX IN LOCAL MATRIX (I).
C     IB        UPPER INDEX IN LOCAL MATRIX (I).
C     MODE      DIRECTION OF TRANSFER (I).
C               = 1 FROM SPARSE STORAGE TO LOCAL ARRAY.
C               = 2 FROM LOCAL ARRAY TO SPARSE STORAGE.
C     CUTF      CUTOFF (I).
C     *
C     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT NONE
      INTEGER :: N,IA,IB,K,J,NFSNEW,NFS,L,I,MODE
      DOUBLE PRECISION :: CUTF,F(N,*)
      DOUBLE PRECISION, DIMENSION (:), POINTER :: FS
      INTEGER, DIMENSION (:), POINTER :: IFS,JFS
C
C *** TRANSFER FROM SPARSE STORAGE TO LOCAL ARRAY.
      IF(MODE.LE.1) THEN
         DO 20 J=1,IB-IA+1
         DO 10 I=1,N
         F(I,J) = 0.0D0
   10    CONTINUE
   20    CONTINUE
C
         DO 40 K=IA,IB
         J      = K-IA+1
         DO 30 L=IFS(K),IFS(K+1)-1
         I      = JFS(L)
         IF(ABS(FS(L)).GT.CUTF) THEN
            F(I,J) = FS(L)
         ENDIF
   30    CONTINUE
   40    CONTINUE
C
         NFS = IFS(N+1)-1
      ELSE
C *** TRANSFER FROM LOCAL ARRAY TO SPARSE STORAGE.
         IF(MODE.EQ.2) THEN
            DO 60 K=IA,IB
            J      = K-IA+1
            NFSNEW = NFS+1
            DO 50 I=1,N
            IF(ABS(F(I,J)).GT.CUTF) THEN
               NFS = NFS+1
               FS(NFS) = F(I,J)
               JFS(NFS) = I
               IF(NFS.EQ.NFSNEW) THEN
                  IFS(K) = NFS
               ENDIF
            ENDIF
   50       CONTINUE
   60       CONTINUE
            IFS(IB+1) = NFS+1
         ELSE
            DO 61 K=IA,IB
            J      = K-IA+1
            DO 51 I=1,N
            IF(ABS(F(I,J)).GT.CUTF) THEN
               NFS = NFS+1
            ENDIF
   51       CONTINUE
   61       CONTINUE
         ENDIF
      ENDIF
      RETURN
      END
