      SUBROUTINE DENSAV (NB,A,LM5,MODE)
C     *
C     SAVE THE DENSITY MATRIX ON   FILE NF (MODE=0).
C     READ THE DENSITY MATRIX FROM FILE NF (MODE=1).
C     *
C     CONVENTION FOR THE FILE NUMBER NF:
C     NF=NB      FOR NB.GT.0, THE FILE IS NOT REWOUND.
C     NF=NBF(11) FOR NB.EQ.0, THE FILE IS REWOUND.
C     NF=NBF(10) FOR NB.LT.0, THE FILE IS REWOUND.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL UHF
      COMMON
     ./INOPT2/ IN2(300)
     ./LIMITS/ LM2,LM3,LM4,LM6,LM7,LM8,LM9
     ./LMSCF / LS1(7),LS8,LS9
     ./LMUHF / LU1(7),LU8
     ./NBFILE/ NBF(20)
     ./UHF   / UHF
      DIMENSION A(LM5)
C *** INPUT OPTIONS.
      INTDIR = IN2(55)
      LINDMS = IN2(56)
      INOUT  = IN2(212)
C *** CHECK FOR AVAILABILITY OF DENSITY MATRIX.
      IF(INTDIR.GT.3 .AND. LINDMS.GT.3) RETURN
C *** DEFINE FILE NUMBER.
      IF(NB.GT.0) THEN
         NF  = NB
      ELSE IF(NB.EQ.0) THEN
         NF  = NBF(11)
         REWIND NF
      ELSE
         NF  = NBF(10)
         REWIND NF
      ENDIF
C *** SAVE OR READ DENSITY MATRIX.
      IF(MODE.EQ.0) THEN
         WRITE(NF) (A(I),I=LS8,LS8-1+LM4)
         IF(UHF) WRITE(NF) (A(I),I=LU8,LU8-1+LM4)
      ELSE IF(MODE.EQ.1) THEN
         READ(NF) (A(I),I=LS8,LS8-1+LM4)
         IF(UHF) READ(NF) (A(I),I=LU8,LU8-1+LM4)
         IF(INOUT.GT.0) THEN
            NB1 = NBF(1)
            REWIND NB1
            WRITE(NB1) (A(I),I=LS8,LS8-1+LM4)
            IF(UHF) WRITE(NB1) (A(I),I=LU8,LU8-1+LM4)
         ENDIF
      ENDIF
      RETURN
      END