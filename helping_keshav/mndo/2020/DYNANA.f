      SUBROUTINE DYNANA (LEN)
C     *
C     DYNAMIC MEMORY ALLOCATION FOR ANALYTIC DERIVATIVE CODE.
C     *
C     NOTATION. I=INPUT.
C     LEN       MAXIMUM LENGTH OF ALLOCATABLE BUFFER (I).
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./INOPT1/ IN1(300)
     ./INOPT2/ IN2(300)
     ./LIMITS/ LM2,LM3,LM4,LM6,LM7,LM8,LM9
     ./NBFILE/ NBF(20)
     ./PSDOPT/ RPSOPT(10),IPSOPT(28)
     ./SKIPA / ISKPA
C *** FILE NUMBERS.
      NB6    = NBF(6)
C *** INPUT VARIABLES.
      NMR    = IN2(13)
      IPSANA = IN2(29)
      JPRINT = IN2(42)
C     *
C *** SECTION FOR ANALYTICAL DERIVATIVES.
C     *
C     THE RELEVANT ROUTINES HAVE THEIR OWN DYNAMIC MEMORY ALLOCATION.
C     ONLY THE MAXIMUM LENGTH OF THE ALLOCATABLE BUFFER IS DEFINED HERE.
C     THIS VALUE (IN MEGABYTE) CAN BE CHOSEN BY INPUT FOR ICORE=IN1(92).
C     IT IS CONVERTED TO WORDS AND STORED IN IN2(92) AND IPSOPT(7),
C     THE DEFAULT VALUE BEING GIVEN BY LEN (WORDS).
C     *
      IF(IPSANA.GT.0 .OR. NMR.NE.0 .OR. IN2(121).EQ.4) THEN
C        THE INPUT VARIABLE ICORE=IN1(92) DEFINES THE REQUESTED LENGTH
C        OF THE IN-CORE BUFFER USED IN THE ANALYTICAL DERIVATIVE CODE.
C        CHECK FOR INPUT ERRORS DISABLED.
C        IF(IN1(92)*131072.GT.LEN) THEN
C           WRITE(NB6,500) IN1(92),LEN
C           STOP 'DYNANA'
C        ENDIF
C        INITIALIZE OR REDEFINE THE BUFFER LENGTH (IF NECESSARY).
C        INPUT OF IN1(92) IN MEGABYTE, IN2(92) IN WORDS.
         IN2(92) = IN1(92)*131072
         IF(IN2(92).LE.0) THEN
            IN2(92) = LEN
         ENDIF
C        SAVE THE BUFFER LENGTH.
         IPSOPT(7)  = IN2(92)
C        PRINTING SECTION.
         IF(JPRINT.GE.2) THEN
            WRITE(NB6,510) IPSOPT(7)
            IF(IN1(92).GT.LEN) WRITE(NB6,520)
         ENDIF
C        CHECK DEFAULT SELECTIONS.
         IF(IN1(29).EQ.0) THEN
            IF(ISKPA.EQ.-1) THEN
               MEMREQ = 20*LM2**2
            ELSE IF(ISKPA.EQ.-2) THEN
               MEMREQ = 30*LM2**2
            ELSE
               MEMREQ = 0
            ENDIF
            IF(IN2(92).LT.MEMREQ) THEN
               WRITE(NB6,530) IN2(92),MEMREQ
               WRITE(NB6,540)
               IF(ISKPA.EQ.-1) WRITE(NB6,550)
               IF(ISKPA.EQ.-2) WRITE(NB6,560)
               STOP 'DYNANA'
            ENDIF
         ENDIF
      ENDIF
      RETURN
C 500 FORMAT(///1X,'MEMORY REQUESTED FOR ANALYTICAL DERIVATIVES',I10,
C    1       /  1X,'MAXIMUM AVAILABLE MEMORY (WORDS)           ',I10,
C    2       /  1X,'PROGRAM WILL STOP.'//)
  510 FORMAT(///1X,'MEMORY ALLOCATION FOR ANALYTICAL DERIVATIVES',
     1       // 1X,'MAXIMUM LENGTH OF ALLOCATABLE BUFFER       ',I10)
  520 FORMAT(   1X,'WARNING: REQUESTED MEMORY LARGER THAN DEFAULT')
  530 FORMAT(///1X,'ANALYTICAL DERIVATIVES WERE SELECTED BY DEFAULT.',
     1       /  1X,'AVAILABLE LENGTH OF BLANK COMMON           ',I10,
     2       /  1X,'REQUIRED  LENGTH FOR EFFICIENT COMPUTATION ',I10)
  540 FORMAT(///1X,'PLEASE CONSIDER THE FOLLOWING INPUT OPTIONS.',
     1       /  1X,'IPSANA=-1 WILL SELECT NUMERICAL DERIVATIVES.',
     2       /  1X,'IPSANA= 1 WILL ATTEMPT A LESS EFFICIENT ANALYTICAL',
     3          1X,'COMPUTATION USING LESS MEMORY.')
  550 FORMAT(   1X,'BOTH OPTIONS WILL REQUIRE LARGE EXECUTION TIMES.')
  560 FORMAT(   1X,'IPSANA= 1 SHOULD STILL BE FASTER IF THE JOB FITS',
     1          1X,'INTO MEMORY.',
     2          1X,'IPSANA=-1 SHOULD BE FEASIBLE OTHERWISE.')
      END
