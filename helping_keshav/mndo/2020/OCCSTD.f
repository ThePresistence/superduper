      SUBROUTINE OCCSTD
C     *
C     DEFINE STANDARD OCCUPATION NUMBERS.
C     *
C     NOTATION IN COMMON BLOCK OCCFL.
C     IMOCC    EQUAL TO ABS(IUHF) (INPUT VARIABLE IUHF).
C              STANDARD OCCUPATION NUMBERS ARE DEFINED FOR IMOCC.LE.2.
C     NOCCA    NUMBER OF HIGHEST OCCUPIED RHF OR UHF-ALPHA ORBITAL.
C     NOCCB    NUMBER OF HIGHEST OCCUPIED UHF-BETA ORBITAL.
C     *
      USE LIMIT, ONLY: LMX
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL UHF
      COMMON
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./HALFE / IOD,JOD
     ./OCCFL / IMOCC,NOCCA,NOCCB,MSUB,MOSUMA,MOSUMB,MOCCA(8),MOCCB(8)
     ./OCCNM / OCCA(LMX),OCCB(LMX)
     ./ORBITS/ NUMB,NORBS,NMOS,NALPHA,NBETA
     ./UHF   / UHF
      IF(IMOCC.GT.2 .AND. IMOCC.LT.5) RETURN
      IF(UHF) THEN
         NOCCA  = NALPHA
         DO 10 I=1,NOCCA
         OCCA(I) = ONE
   10    CONTINUE
         DO 20 I=NOCCA+1,NORBS
         OCCA(I) = ZERO
   20    CONTINUE
         NOCCB  = NBETA
         DO 30 I=1,NOCCB
         OCCB(I) = ONE
   30    CONTINUE
         DO 40 I=NOCCB+1,NORBS
         OCCB(I) = ZERO
   40    CONTINUE
      ELSE
         NOCCA  = NUMB
         DO 50 I=1,NUMB
         OCCA(I) = ONE
   50    CONTINUE
         IF(IOD.GT.0) THEN
            OCCA(IOD) = PT5
            IF(JOD.GT.0) OCCA(JOD)=PT5
         ENDIF
         DO 60 I=NUMB+1,NORBS
         OCCA(I) = ZERO
   60    CONTINUE
      ENDIF
      RETURN
      END
