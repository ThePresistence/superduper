      SUBROUTINE BNDPOP(PA,PB,LM4,LM2,BOMAT,DPRINT)
C     CALCULATES THE BOND ORDER MATRIX
C     FOR BOTH UHF AND RHF WAVEFUNCTIONS
C     NOTE:
C     SETS THE DIAGONAL ELEMENTS TO UNITY
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (LM1=600)
      PARAMETER (LMX=9*LM1)
      PARAMETER (LMZ=86)
      INTEGER LM4,LM2
      DIMENSION PA(LM4),PB(LM4),F(LM2,LM2)
      LOGICAL UHF
      LOGICAL DPRINT
C      PARAMETER (LMPAIR=LM1*(LM1-1)/2)
      COMMON
C     ./ATOMC / COORD(3,LM1)
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
C.##IF CHARMM
C.##ELSE
     ./INDEX / INDX(LMX)
C.##ENDIF
     ./ORBITS/ NUMB,NORBS,NMOS,NALPHA,NBETA
C     ./CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./UHF   / UHF
C     ./PARDER/ TORE(LMZ),EHEAT(LMZ),EISOL(LMZ)
      DIMENSION BOMAT(NUMAT,NUMAT)
      DOUBLE PRECISION X
      INTEGER IL,IH,IA,IB,JA,JB
c      DPRINT = .FALSE.
C      WRITE(6,*)"A",DPRINT
      IF ( DPRINT ) WRITE(6,'(A)')"_____ENTERED_BNDPOP____"
      DO 1 I=1,NUMAT
         DO 2 J=1,NUMAT
            BOMAT(I,J) = ZERO
 2       CONTINUE
 1    CONTINUE
      DO 3 I=1,LM2
         DO 4 J=1,LM2
            F(I,J) = ZERO
 4       CONTINUE
 3    CONTINUE
C      WRITE(6,*)"B",DPRINT
      IF ( UHF ) THEN
         FACT = TWO
      ELSE
         FACT = FOUR
      END IF
C      IF ( DPRINT ) THEN
C         WRITE(6,'(A)')"i     PA(i)      PB(i)"
C         DO 7 I=1,LM4
C            WRITE(6,'(I3,2F10.5)')I,PA(I),PB(I)
C 7       CONTINUE
C      END IF
C      WRITE(6,*)NORBS
      CALL SQUARE (PA,F,NORBS,LM2,LM4)
C      DO 1000 I=1,LM2
C         DO 1001 J=1,LM2
C            WRITE(6,*)I,J,F(I,J)
C 1001    CONTINUE
C 1000 CONTINUE
C      WRITE(6,*)"C",DPRINT
C      IF ( DPRINT ) THEN
C         WRITE(6,'(A)')"F MATRIX (PA)"
C         DO 8 I=1,LM2
C            WRITE(6,'(I3,100F10.5:)')I,(F(I,J),J=1,LM2)
C 8       CONTINUE
C      END IF
      DO 40 I=1,NUMAT
         IA     = NFIRST(I)
         IB     = NLAST(I)
         DO 20 J=1,I
            JA     = NFIRST(J)
            JB     = NLAST(J)
            X      = ZERO
            DO 11 IL=IA,IB
               DO 10 IH=JA,JB
C                  WRITE(6,*)LM2,IL,IH,DPRINT,IA,IB,JA,JB
                  X      = X+F(IL,IH)*F(IL,IH)
C                  WRITE(6,*)DPRINT
 10            CONTINUE
 11         CONTINUE
            BOMAT(I,J)  = X*FACT
            BOMAT(J,I)  = BOMAT(I,J)
C            IF ( DPRINT ) THEN
C               WRITE(6,'(A,2I3,2F10.5)')"I,J,BO,X",I,J,BOMAT(I,J),X
C            END IF
  20      CONTINUE
         X      =-BOMAT(I,I)
C         IF ( DPRINT ) THEN
C            WRITE(6,'(A,I3,F10.5)')"I,X",I,X
C         END IF
         DO 30 J=IA,IB
            X      = X+FACT*F(J,J)
 30      CONTINUE
C         IF ( DPRINT ) THEN
C            WRITE(6,'(A,2I3,F10.5)')"I,J,X",I,J,X
C         END IF
C         WRITE(6,*)"BNDPOP",NUMAT,I,J
         BOMAT(I,I)  = X
C         BOMAT(J,I)  = BOMAT(I,J)
 40   CONTINUE
C      WRITE(6,*)"D",DPRINT
      IF ( .NOT. UHF ) GOTO 999
      DO 5 I=1,LM2
         DO 6 J=1,LM2
            F(I,J) = ZERO
 6       CONTINUE
 5    CONTINUE
C      WRITE(6,*)"A",DPRINT
C      WRITE(6,*)NORBS
      CALL SQUARE (PB,F,NORBS,LM2,LM4)
C      WRITE(6,*)"A",DPRINT
      DO 400 I=1,NUMAT
         IA     = NFIRST(I)
         IB     = NLAST(I)
         DO 200 J=1,I
            JA     = NFIRST(J)
            JB     = NLAST(J)
            X      = ZERO
            DO 110 IL=IA,IB
               DO 100 IH=JA,JB
                  X      = X+F(IL,IH)*F(IL,IH)
 100           CONTINUE
 110        CONTINUE
            BOMAT(I,J)  = BOMAT(I,J)+X*FACT
            BOMAT(J,I)  = BOMAT(I,J)
C            TEMP  = X*FACT
 200     CONTINUE
         X      =-BOMAT(I,I)
         DO 300 J=IA,IB
            X      = X+FACT*F(J,J)
 300     CONTINUE
         BOMAT(I,I)  = X
C         BOMAT(J,I)  = BOMAT(I,J)
 400  CONTINUE
 999  CONTINUE
      DO 450 I=1,NUMAT
         BOMAT(I,I) = ONE
 450  CONTINUE
C      WRITE(6,*)DPRINT
      IF ( DPRINT ) THEN
         WRITE(6,'(A)')"BOND ORDER MATRIX"
         DO 500 I=1,NUMAT
            DO 501 J=1,NUMAT
               WRITE(6,'(2I3,20F11.5:)')I,J,BOMAT(I,J)
 501        CONTINUE
 500     CONTINUE
         WRITE(6,'(A)')"_____EXITING_BNDPOP____"
      END IF
      END SUBROUTINE
