      SUBROUTINE PQDISP(EE,PA,PB,LM4,LM2)
C     SUBROUTINE CALCULATES THE PELLENQ
C     DISPERSION ENERGY CONTRIBUTION TO
C     THE ELECTRONIC ENERGY
C     AS A POST-SCF CORRECTION
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (LM1=600)
      PARAMETER (LMX=9*LM1)
      PARAMETER (LMZ=86)
      PARAMETER (LPQ=LMZ+LMZ*(LMZ-1)/2)
      DIMENSION PA(LM4),PB(LM4)
      DIMENSION EFFN(LM1)
      DIMENSION SL(3),ETAL(3,LM1), TTDAMP(6), PQC(6), POLL(3)
      DIMENSION THKROM(3)
      LOGICAL DPRINT
      COMMON
     ./ATOMC / COORD(3,LM1)
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./PARDER/ TORE(LMZ),EHEAT(LMZ),EISOL(LMZ)
C     PELLENQ DISPERSION PARAMETERS
      COMMON
     ./PQPAR / PQNEF(LMZ),PQDP(LMZ),PQQP(LMZ),PQOP(LMZ),PQB(LPQ),
     .         PQBO0(LMZ),PQBOZ(LMZ),PQRL(LMZ),PQRH(LMZ),PQDS,PQQS,PQOS
     ./PQOPT / NFLPQD,NPTPQD
     ./PQOUT / PQENE
      DIMENSION BOMAT(NUMAT,NUMAT)
C     DEBUG PRINTING FLAG
      IF ( NPTPQD .EQ. 1 ) THEN
         DPRINT = .TRUE.
      ELSE
         DPRINT = .FALSE.
      END IF
C     THAKKER OMEGA COEFFICIENTS TO CALCULATE
C     C12, C14, and C16
C     See Thakker_JChemPhys_1988_v89_p2092 equation 29
      THKROM(1) = 1.0281430356399119D0
      THKROM(2) = 0.9749145042676423D0
      THKROM(3) = 0.9735032057454021D0
      IF ( DPRINT ) THEN
         WRITE(6,'(1X)')
         WRITE(6,'(A)')"*******************************"
         WRITE(6,'(A)')"        ENTERED PQDISP"
         WRITE(6,'(A)')" CALCULATES PELLENQ DISPERSION"
         WRITE(6,'(A)')" ENERGY CORRECTION"
         WRITE(6,'(A)')"*******************************"
      END IF
      PQENE = ZERO
C     GET THE EFFECTIVE NUMBER OF ELECTRONS
      CALL MULPOP(PA,PB,LM4,EFFN)
C     GET THE BOND ORDER MATRIX
C      WRITE(6,'(A)')"CALL BNDPOP"
      IF ( NFLPQD .EQ. 2 ) THEN
         CALL BNDPOP(PA,PB,LM4,LM2,BOMAT,.FALSE.)
         IF ( DPRINT ) THEN
            WRITE(6,'(A)')"BOND ORDER MATRIX"
            DO 3 I=1,NUMAT
               WRITE(6,'(100(F8.4:))')(BOMAT(I,J),J=1,NUMAT)
 3          CONTINUE
         END IF
C     CREATE A DAMPING MATRIX BASED ON THE BOND ORDER
         CALL BNDDMP(BOMAT)
         IF ( DPRINT ) THEN
            WRITE(6,'(A)')"BOND ORDER DAMPING MATRIX"
            DO 4 I=1,NUMAT
               WRITE(6,'(100F8.4:)')(BOMAT(I,J),J=1,NUMAT)
 4          CONTINUE
         END IF
      ELSE
         DO 5 I=1,NUMAT
            DO 6 J=1,NUMAT
               BOMAT(I,J) = ONE
 6          CONTINUE
 5       CONTINUE
      END IF
      SVNPT5 = THREE*(THREE+TWO)/TWO
      DO 20 I=1,NUMAT
         POLL(1) = PQDP(NAT(I))*PQDS
         POLL(2) = PQQP(NAT(I))*PQQS
         POLL(3) = PQOP(NAT(I))*PQOS
         PQLAM = PQNEF(NAT(I)) / TORE(NAT(I))
         VNE = EFFN(I) * PQLAM
         IF ( VNE .LT. ZERO ) VNE = ZERO
         IF ( DPRINT ) THEN
            WRITE(6,'(A,I3,A,F11.5,A,F11.5,A,F11.5)')" ATOM= ",I
     A           ," LAM= ",PQLAM," NVAL= ",
     B           EFFN(I)," NEFF= ",VNE
         END IF
         SL(1) = VNE
         SL(2) = THREE * SQRT( VNE*POLL(1) )
         SL(3) = SVNPT5 * SQRT( SL(2)*POLL(2) )
         DO 25 J=1,3
            IF ( POLL(J).GT.ZERO ) THEN
               ETAL(J,I) = SQRT( SL(J) / POLL(J) )
            ELSE
               ETAL(J,I) = SQRT( SL(J) )
            END IF
            IF ( DPRINT ) THEN
               WRITE(6,'(A,2I4,3F12.4)')"I,ZI,POL(J),SL(J),ETAL(J,I)",
     A              I,NAT(I),POLL(J),SL(J),ETAL(J,I)
C               WRITE(6,'(A,2F12.4)')"SVNPT5,SQRT",
C     A              SVNPT5,SQRT( SL(2)*POLL(2) )
            END IF
 25      CONTINUE
 20   CONTINUE
      PREFA = THREE/TWO
      PREFB = THREE*(THREE+TWO)/FOUR
      PREFC = THREE+FOUR
      PREFD = (THREE+TWO)/TWO
      ICNT = 0
      DO 50 I=1,NUMAT-1
         DPI = PQDP(NAT(I))*PQDS
         QPI = PQQP(NAT(I))*PQQS
         OPI = PQOP(NAT(I))*PQOS
         DO 55 J=I+1,NUMAT
            ICNT=ICNT+1
            DPJ = PQDP(NAT(J))*PQDS
            QPJ = PQQP(NAT(J))*PQQS
            OPJ = PQOP(NAT(J))*PQOS
            PQC(1) = PREFA*ABOAPB(ETAL(1,I),ETAL(1,J))*DPI*DPJ
C            WRITE(6,*)"C6 RAW",PQC(1),DPI,DPJ,ETAL(1,I),ETAL(1,J),PREFA
            PQC(2) = PREFB * (
     A           ABOAPB(ETAL(1,I),ETAL(2,J))*DPI*QPJ +
     B           ABOAPB(ETAL(1,J),ETAL(2,I))*DPJ*QPI )
            PQC(3) = PREFC * (
     A           ABOAPB(ETAL(1,I),ETAL(3,J))*DPI*OPJ +
     B           ABOAPB(ETAL(1,J),ETAL(3,I))*DPJ*OPI +
     C           PREFD * ABOAPB(ETAL(2,J),ETAL(2,I))*QPJ*QPI )
C     CALCULATE C12,14,16 FROM C6,8,10 USING THE THAKKER COEFFICIENTS
            DO 27 K=4,6
               PQC(K)=THKROM(K-3)*PQC(K-3)*((PQC(K-1)/PQC(K-2))**THREE)
 27         CONTINUE
            RIJ = SQRT(
     A           ( COORD(1,I) - COORD(1,J) )**2 +
     B           ( COORD(2,I) - COORD(2,J) )**2 +
     C           ( COORD(3,I) - COORD(3,J) )**2 )/A0
C---------------------------------------------------
C     SCALE THE EFFECTIVE DISTANCE
            RIJ0 = RIJ
            IF ( NFLPQD .EQ. 3 ) THEN
               RL = (PQRL(NAT(I))+PQRL(NAT(J)))/2.0D0
               RH = (PQRH(NAT(I))+PQRH(NAT(J)))/2.0D0
C     CONVERT FROM ANGSTROM TO BOHR
               RL = RL/A0
               RH = RH/A0
               IF ( RIJ .GE. RH ) THEN
                  RIJ = RIJ
               ELSE IF ( RIJ < RL ) THEN
                  RIJ = (RL+RH)/2.0D0
               ELSE
                  RS = (RL-RIJ)/(RH-RL)
                  RS2 = RS*RS
                  RS3 = RS2*RS
                  RIJ = (RL-RIJ)*RS3*(2.5D0+THREE*RS+RS2)+PT5*(RL+RH)
               END IF
            END IF
C-----------------------------------
C COMPUTE THE INDEX FOR THE IJ PAIR AND CALL IT LIJ
            LATI = NAT(I)
            LATJ = NAT(J)
            IF ( LATI.LT.LATJ ) THEN
               LIJ = LATJ + (LATI-1)*(2*LMZ-LATI)/2
            ELSE
               LIJ = LATI + (LATJ-1)*(2*LMZ-LATJ)/2
            END IF
C            IF ( DPRINT ) WRITE(6,'(A,1X,I5)')"LIJ = ",LIJ
C     CONVERT FROM 1/A to 1/Bohr
C     A0 has units of A/Bohr
            BIJ = PQB(LIJ)*A0
C     CALCULATE THE TANG-TOENNIES DAMPING
            CALL TTDMPF( RIJ,BIJ,TTDAMP )
C            TTDAMP(4) = ZERO
C            TTDAMP(5) = ZERO
C            TTDAMP(6) = ZERO
            IF ( DPRINT ) THEN
               WRITE(6,'(A,9F11.5)')"B,R,DAMP=",BIJ,RIJ0,RIJ,
     A              (TTDAMP(K),K=1,6)
            END IF
            IF ( DPRINT ) THEN
               WRITE(6,'(A,2I3,6(1X,E12.4))')"C6-16,I,J:",
     A              I,J,(PQC(K),K=1,6)
            END IF
C            PQOLD = PQENE
            DENE = 0.0D0
            DO 60 K=1,6
C               IF ( DPRINT ) THEN
C                  WRITE(6,'(A,I2,A,I3,I3,E12.4)')"I,J, C(",
C     A                 2*K+4," )",LATI,LATJ,PQC(K)
C               END IF
               DENE=DENE-BOMAT(I,J)*TTDAMP(K)*PQC(K)/RIJ**(TWO*(K+2))
C               PQENE = PQENE + DENE(K)
C               PQENE = PQENE - BOMAT(I,J)*TTDAMP(K)*
C     A              PQC(K)/RIJ**(TWO*(K+2))
 60         CONTINUE
            PQENE = PQENE+DENE
            IF ( NPTPQD .GT. 0 ) THEN
               WRITE(6,'(A,2I3,2E12.4)')"I,J,RIJ,Disp:",
     A              I,J,RIJ,DENE
            END IF
 55      CONTINUE
 50   CONTINUE
      IF ( DPRINT ) WRITE(6,'(A,E20.12)')"   PQENE (AU) = ",PQENE
      PQENE = PQENE * EV
C      WRITE(6,*)"FOO"
      IF ( NPTPQD .GT. 0 ) THEN
         WRITE(6,'(A,E20.12)')"   PQENE (EV) = ",PQENE
      END IF
C      IF ( DPRINT ) WRITE(6,'(A,E20.12)')"EE       (EV) = ",EE
      EE = EE + PQENE
C      IF ( DPRINT ) WRITE(6,'(A,E20.12)')"EE+PQENE (EV) = ",EE
C      WRITE(6,*)"EXITING PELLENQ"
      END SUBROUTINE
      SUBROUTINE TTDMPF(R,B,TT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION TT(6)
      COMMON
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
      DO 10 I=1,6
         TT(I) = ZERO
 10   CONTINUE
      BR = B*R
      EXPBR = EXP( -BR )
      KFACT = ONE
      DFACT = ONE
      DO 20 I=1,6,1
         KFACT = KFACT*I
         DFACT = DFACT*BR
         TT(1) = TT(1) + DFACT/KFACT
 20   CONTINUE
C     THIS IS THE K=0 PART
      TT(1) = TT(1) + ONE
      TT(2)=TT(1)
      DO 30 I=7,8,1
         KFACT = KFACT*I
         DFACT = DFACT*BR
         TT(2) = TT(2) + DFACT/KFACT
 30   CONTINUE
      TT(3) = TT(2)
      DO 40 I=9,10,1
         KFACT = KFACT*I
         DFACT = DFACT*BR
         TT(3) = TT(3) + DFACT/KFACT
 40   CONTINUE
      TT(4) = TT(3)
      DO 50 I=11,12,1
         KFACT = KFACT*I
         DFACT = DFACT*BR
         TT(4) = TT(4) + DFACT/KFACT
 50   CONTINUE
      TT(5) = TT(4)
      DO 60 I=13,14,1
         KFACT = KFACT*I
         DFACT = DFACT*BR
         TT(5) = TT(5) + DFACT/KFACT
 60   CONTINUE
      TT(6) = TT(5)
      DO 70 I=15,16,1
         KFACT = KFACT*I
         DFACT = DFACT*BR
         TT(6) = TT(6) + DFACT/KFACT
 70   CONTINUE
      DO 80 I=1,6
C         WRITE(6,*)"TTDAMP 2n=",2*I+4,TT(I),EXPBR,TT(I)*EXPBR
         TT(I) = ONE - TT(I)*EXPBR
         IF ( TT(I) .LT. ZERO ) TT(I) = ZERO
 80   CONTINUE
      END SUBROUTINE
C**************************************************************
C**************************************************************
C**************************************************************
      FUNCTION ABOAPB(A,B)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      C = A+B
      IF (C.NE.0.0D0) THEN
         ABOAPB = A*B/C
      ELSE
         ABOAPB = 1.D10
      END IF
C      WRITE(6,*)"ABOAPB = ",ABOAPB
      RETURN
      END FUNCTION
C*************************************************************
C*************************************************************
C*************************************************************
      SUBROUTINE BNDDMP(BOMAT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (LM1=600)
      PARAMETER (LMX=9*LM1)
      PARAMETER (LMZ=86)
      PARAMETER (LPQ=LMZ+LMZ*(LMZ-1)/2)
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
      COMMON
     ./PQPAR / PQNEF(LMZ),PQDP(LMZ),PQQP(LMZ),PQOP(LMZ),PQB(LPQ),
     .     PQBO0(LMZ),PQBOZ(LMZ),PQRL(LMZ),PQRH(LMZ),PQDS,PQQS,PQOS
      DIMENSION BOMAT(NUMAT,NUMAT),BNDSCL(NUMAT,NUMAT)
      BLO = ZERO
      BHI = ONE/THREE
C     LET'S TAKE IT TO THE NEXT LEVEL
C     GOD THAT SOUNDS CORNY
      DO 100 I=1,NUMAT-1
         DO 110 J=I+1,NUMAT
C------------------------------------------------------------------
C     1-2 interactions
            B = BOMAT(I,J)
            IF ( B .GE. BHI ) GOTO 999
C     1-3 interactions
            DO 120 K=1,NUMAT
               IF ( K .NE. I .AND. K .NE. J ) THEN
                  Z = BOMAT(I,K)*BOMAT(K,J)
                  IF ( Z .GT. B ) THEN
                     B = Z
                     IF ( B .GT. BHI ) GOTO 999
                  END IF
               END IF
 120        CONTINUE
 999        CONTINUE
            BNDSCL(I,J) = B
            BNDSCL(J,I) = B
 110     CONTINUE
 100  CONTINUE
      DO 10 I=1,NUMAT-1
         DO 20 J=I+1,NUMAT
            B = BNDSCL(I,J)
            CALL GENSW( B,BLO,BHI,SW )
            B = ONE-SW
            BOMAT(I,J)=B
            BOMAT(J,I)=B
 20      CONTINUE
 10   CONTINUE
      DO 30 I=1,NUMAT
         BOMAT(I,I) = ZERO
 30   CONTINUE
      END SUBROUTINE
c$$$      SUBROUTINE BNDDMP(BOMAT)
c$$$      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c$$$      PARAMETER (LM1=600)
c$$$      PARAMETER (LMX=9*LM1)
c$$$      PARAMETER (LMZ=86)
c$$$      PARAMETER (LPQ=LMZ+LMZ*(LMZ-1)/2)
c$$$      COMMON
c$$$     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
c$$$     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
c$$$      COMMON
c$$$     ./PQPAR / PQNEF(LMZ),PQDP(LMZ),PQQP(LMZ),PQOP(LMZ),PQB(LPQ),
c$$$     .         PQBO0(LMZ),PQBOZ(LMZ)
c$$$      DIMENSION BOMAT(NUMAT,NUMAT)
c$$$
c$$$C CALCULATE PELLENQ DAMPING BASED ON BOND ORDER
c$$$
c$$$      DO 10 I=1,NUMAT
c$$$         RI = PQBO0(NAT(I))
c$$$         ZI = PQBOZ(NAT(I))
c$$$         DO 20 J=1,I-1
c$$$            RJ = PQBO0(NAT(J))
c$$$            ZJ = PQBOZ(NAT(J))
c$$$            R = (RI+RJ)/TWO
c$$$            IF ( R .LT. ZERO ) R = ZERO
c$$$            Z = SQRT(ZI*ZJ)
c$$$            IF ( BOMAT(I,J) .LT. ZERO ) THEN
c$$$C     THIS SHOULDN'T HAPPEN, BUT JUST IN CASE
c$$$               BOMAT(I,J) = ZERO
c$$$            END IF
c$$$            IF ( BOMAT(I,J) .LE. R ) THEN
c$$$               BOMAT(I,J) = ONE
c$$$            ELSE
c$$$               RIJ = BOMAT(I,J)-R
c$$$               RIJ2 = RIJ*RIJ
c$$$               BOMAT(I,J) = EXP(-Z*RIJ2)
c$$$            END IF
c$$$C     KEEP IT SYMMETRIC
c$$$            BOMAT(J,I) = BOMAT(I,J)
c$$$
c$$$ 20      CONTINUE
c$$$         BOMAT(I,I) = ZERO
c$$$ 10   CONTINUE
c$$$
c$$$      END SUBROUTINE
C*************************************************************
C*************************************************************
C*************************************************************
      SUBROUTINE MULPOP(PA,PB,LM4,VALPOP)
C     SUBROUTINE CALCULATES THE NUMBER OF VALENCE ELECTRONS
C     FROM A MULLIKEN POPULATION ANALYSIS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (LM1=600)
      PARAMETER (LMX=9*LM1)
      PARAMETER (LMZ=86)
      DIMENSION PA(LM4),PB(LM4),VALPOP(LM1)
      LOGICAL UHF
      PARAMETER (LMPAIR=LM1*(LM1-1)/2)
      COMMON
C     .     /ATOMC / COORD(3,LM1)
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
C.##IF CHARMM
C.##ELSE
     ./INDEX / INDX(LMX)
C.##ENDIF
C     .     /CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./UHF   / UHF
     ./PARDER/ TORE(LMZ),EHEAT(LMZ),EISOL(LMZ)
      IF(UHF) THEN
      DO 30 I=1,NUMAT
         VALPOP(I)= ZERO
         IA      = NFIRST(I)
         LL      = INDX(IA)+IA
         VALPOP(I)= VALPOP(I)+PA(LL)+PB(LL)
C         VALPOP(I)= VALPOP(I)+TORE(NAT(I))
         IORBS   = NLAST(I)-IA+1
         IF(IORBS.GE.4) THEN
            DO 21 J= 1,3
            LL      = INDX(J+IA)+J+IA
   21       VALPOP(I)= VALPOP(I)+PA(LL)+PB(LL)
            IF(IORBS.GE.9) THEN
               DO 22 J= 4,8
               LL      = INDX(J+IA)+J+IA
   22          VALPOP(I)= VALPOP(I)+PA(LL)+PB(LL)
            ENDIF
         ENDIF
   30 CONTINUE
      ELSE
      DO 50 I=1,NUMAT
         VALPOP(I)= ZERO
         IA      = NFIRST(I)
         LL      = INDX(IA)+IA
         VALPOP(I)= VALPOP(I)+TWO*PA(LL)
C         VALPOP(I)= VALPOP(I)+TORE(NAT(I))
         IORBS   = NLAST(I)-IA+1
         IF(IORBS.GE.4) THEN
            DO 41 J= 1,3
            LL      = INDX(J+IA)+J+IA
   41       VALPOP(I)= VALPOP(I)+TWO*PA(LL)
            IF(IORBS.GE.9) THEN
              DO 42 J= 4,8
              LL      = INDX(J+IA)+J+IA
   42         VALPOP(I)= VALPOP(I)+TWO*PA(LL)
            ENDIF
         ENDIF
   50 CONTINUE
      END IF
      END SUBROUTINE
