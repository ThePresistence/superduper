      SUBROUTINE DENEXT (NB,A,LM5)
C     *
C     EXTRAPOLATE THE DENSITY MATRIX IN FINITE-DIFFERENCE PROCEDURES.
C     *
C     THERE ARE THREE RELEVANT GEOMETRIES:
C     A.    REFERENCE GEOMETRY.
C     B.    A GIVEN VARIABLE WITH A POSITIVE DISPLACEMENT.
C     C.    A GIVEN VARIABLE WITH A NEGATIVE DISPLACEMENT.
C     CORRESPONDING DENSITY MATRICES:
C     A.    AVAILABLE FROM FILE NF.
C     B.    AVAILABLE AS CURRENT DENSITY MATRIX AT A(LS8).
C     C.    ESTIMATED AS P(C)=P(A)-(P(B)-P(A))=TWO*P(A)-P(B).
C     *
C     CONVENTION FOR THE FILE NUMBER NF:
C     NF=NB      FOR NB.GT.0, THE FILE IS NOT REWOUND.
C     NF=NBF(11) FOR NB.EQ.0, THE FILE IS REWOUND.
C     NF=NBF(10) FOR NB.LT.0, THE FILE IS REWOUND.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (TWO=2.0D0)
      LOGICAL UHF
      COMMON
     ./INOPT2/ IN2(300)
     ./LIMITS/ LM2,LM3,LM4,LM6,LM7,LM8,LM9
     ./LMSCF / LS1(6),LS7,LS8,LS9
     ./LMUHF / LU1(7),LU8
     ./NBFILE/ NBF(20)
     ./UHF   / UHF
      DIMENSION A(LM5)
C *** INPUT OPTION.
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
C *** RHF (OR UHF ALPHA) DENSITY MATRIX.
C     COPY P(B) FROM A(LS8) TO A(LS7).
      DO 10 I=1,LM4
      A(LS7-1+I) = A(LS8-1+I)
   10 CONTINUE
C     READ P(A) INTO A(LS8).
      READ(NF) (A(I),I=LS8,LS8-1+LM4)
C     FORM P(C) AT A(LS8).
      DO 20 I=1,LM4
      A(LS8-1+I) = TWO*A(LS8-1+I)-A(LS7-1+I)
   20 CONTINUE
C *** UHF BETA DENSITY MATRIX.
      IF(UHF) THEN
         DO 30 I=1,LM4
         A(LS7-1+I) = A(LU8-1+I)
   30    CONTINUE
         READ(NF) (A(I),I=LU8,LU8-1+LM4)
         DO 40 I=1,LM4
         A(LU8-1+I) = TWO*A(LU8-1+I)-A(LS7-1+I)
   40    CONTINUE
      ENDIF
C *** SAVE EXTRAPOLATED DENSITY ON FILE NB1 (IF NEEDED).
      IF(INOUT.GT.0) THEN
         NB1 = NBF(1)
         REWIND NB1
         WRITE(NB1) (A(I),I=LS8,LS8-1+LM4)
         IF(UHF) WRITE(NB1) (A(I),I=LU8,LU8-1+LM4)
      ENDIF
      RETURN
      END
