      SUBROUTINE QPOP(QAT,PA,PB,LM4)
C     SUBROUTINE CALCULATES THE MULLIKEN PARTIAL CHARGE
C     OF EACH ATOM FROM THE DENSITY MATRIX
C     THIS ROUTINE WORKS FOR BOTH UHF AND RHF WAVEFUNCTIONS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (LM1=600)
      PARAMETER (LMX=9*LM1)
      PARAMETER (LMZ=86)
      DIMENSION PA(LM4),PB(LM4), QAT(LM1)
      LOGICAL UHF
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
C.##IF CHARMM
C.##ELSE
     ./INDEX / INDX(LMX)
C.##ENDIF
     ./CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./UHF   / UHF
     ./PARDER/ TORE(LMZ),EHEAT(LMZ),EISOL(LMZ)
C     Calculate mulliken point charges
      IF(UHF) THEN
      DO 30 I=1,NUMAT
         QAT(I)= ZERO
         IA      = NFIRST(I)
         LL      = INDX(IA)+IA
         QAT(I)= QAT(I)-PA(LL)-PB(LL)
         QAT(I)= QAT(I)+TORE(NAT(I))
         IORBS   = NLAST(I)-IA+1
         IF(IORBS.GE.4) THEN
            DO 21 J= 1,3
            LL      = INDX(J+IA)+J+IA
   21       QAT(I)= QAT(I)-PA(LL)-PB(LL)
            IF(IORBS.GE.9) THEN
               DO 22 J= 4,8
               LL      = INDX(J+IA)+J+IA
   22          QAT(I)= QAT(I)-PA(LL)-PB(LL)
            ENDIF
         ENDIF
   30 CONTINUE
      ELSE
      DO 50 I=1,NUMAT
         QAT(I)= ZERO
         IA      = NFIRST(I)
         LL      = INDX(IA)+IA
         QAT(I)= QAT(I)-TWO*PA(LL)
         QAT(I)= QAT(I)+TORE(NAT(I))
         IORBS   = NLAST(I)-IA+1
         IF(IORBS.GE.4) THEN
            DO 41 J= 1,3
            LL      = INDX(J+IA)+J+IA
   41       QAT(I)= QAT(I)-TWO*PA(LL)
            IF(IORBS.GE.9) THEN
              DO 42 J= 4,8
              LL      = INDX(J+IA)+J+IA
   42         QAT(I)= QAT(I)-TWO*PA(LL)
            ENDIF
         ENDIF
   50 CONTINUE
      END IF
      END SUBROUTINE
