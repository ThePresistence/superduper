      SUBROUTINE SDISP(EE,PA,PB,LM4,LM2)
C     SUBROUTINE CALCULATES THE SIMPLE
C     DISPERSION ENERGY CONTRIBUTION TO
C     THE ELECTRONIC ENERGY
C     AS A POST-SCF CORRECTION
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (LM1=600)
      PARAMETER (LMX=9*LM1)
      PARAMETER (LMZ=86)
C      PARAMETER (LPQ=LMZ+LMZ*(LMZ-1)/2)
      DIMENSION PA(LM4),PB(LM4)
      DIMENSION EFFN(LM1),QCHRG(LM1)
      DIMENSION SL(2),ETAL(2,LM1), SDC(2), POLL(2),RIJV(2)
C      DIMENSION THKROM(3)
      LOGICAL DPRINT
      COMMON
     ./ATOMC / COORD(3,LM1)
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./PARDER/ TORE(LMZ),EHEAT(LMZ),EISOL(LMZ)
     ./SDISPB/ SDRLD(LMZ),SDRLQ(LMZ),SDRHD(LMZ),SDRHQ(LMZ),
     . SDDP(LMZ),SDQP(LMZ),SDDPQ(LMZ),SDQPQ(LMZ),SDNEF(LMZ),NFLSD
C     DEBUG PRINTING FLAG
C      DPRINT = .TRUE.
      DPRINT = .FALSE.
C     THAKKER OMEGA COEFFICIENTS TO CALCULATE
C     C12, C14, and C16
C     See Thakker_JChemPhys_1988_v89_p2092 equation 29
C      THKROM(1) = 1.0281430356399119D0
C      THKROM(2) = 0.9749145042676423D0
C      THKROM(3) = 0.9735032057454021D0
      IF ( DPRINT ) THEN
         WRITE(6,'(1X)')
         WRITE(6,'(A)')"*******************************"
         WRITE(6,'(A)')"        ENTERED SDISP"
         WRITE(6,'(A)')" CALCULATES SIMPLE DISPERSION"
         WRITE(6,'(A)')"      ENERGY CORRECTION"
         WRITE(6,'(A)')"*******************************"
      END IF
      SDENE = ZERO
C     GET THE EFFECTIVE NUMBER OF ELECTRONS
C     MULPOP APPEARS IN PELLENQ.src
C      CALL MULPOP(PA,PB,LM4,EFFN)
      CALL QPOP(PA,PB,LM4,QCHRG)
      DO 10 I=1,NUMAT
         EFFN(I) = TORE(NAT(I)) - QCHRG(I)
 10   CONTINUE
      DO 20 I=1,NUMAT
         POLL(1) = SDDP(NAT(I))*EXP(-SDDPQ(NAT(I))*QCHRG(I))
         POLL(2) = SDQP(NAT(I))*EXP(-SDQPQ(NAT(I))*QCHRG(I))
         SDLAM = SDNEF(NAT(I)) / TORE(NAT(I))
         VNE = EFFN(I) * SDLAM
C RARE GAS DIMER HACK
C         VNE = SDNEF(NAT(I))
C----HACK------------
         IF ( VNE .LT. ZERO ) VNE = ZERO
         IF ( DPRINT ) THEN
            WRITE(6,'(A,I3,A,F11.5,A,F11.5,A,F11.5)')" ATOM= ",I
     A           ," LAM= ",SDLAM," NVAL= ",
     B           EFFN(I)," NEFF= ",VNE
         END IF
         SL(1) = VNE
         SL(2) = THREE * SQRT( VNE*POLL(1) )
         DO 25 J=1,2
            IF ( POLL(J).GT.ZERO ) THEN
               ETAL(J,I) = SQRT( SL(J) / POLL(J) )
            ELSE
               ETAL(J,I) = 0.0D0
CSQRT( SL(J) )
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
         DPI = SDDP(NAT(I))
         QPI = SDQP(NAT(I))
         DO 55 J=I+1,NUMAT
            ICNT=ICNT+1
C            WRITE(6,*)"RL D ",SDRLD(NAT(I)),SDRLD(NAT(J))
C            WRITE(6,*)"RL Q ",SDRLQ(NAT(I)),SDRLQ(NAT(J))
C            WRITE(6,*)"RH D ",SDRHD(NAT(I)),SDRHD(NAT(J))
C            WRITE(6,*)"RH Q ",SDRHQ(NAT(I)),SDRHQ(NAT(J))
            DPJ = SDDP(NAT(J))*EXP(-SDDPQ(NAT(J))*QCHRG(J))
            QPJ = SDQP(NAT(J))*EXP(-SDQPQ(NAT(J))*QCHRG(J))
C----------------------------------------------------
C COMPUTE THE DISPERISON COEFFICIENTS
            SDC(1) = PREFA*ABOAPB(ETAL(1,I),ETAL(1,J))*DPI*DPJ
C            WRITE(6,*)"C6 RAW",SDC(1),DPI,DPJ,ETAL(1,I),ETAL(1,J),PREFA
            SDC(2) = PREFB * (
     A           ABOAPB(ETAL(1,I),ETAL(2,J))*DPI*QPJ +
     B           ABOAPB(ETAL(1,J),ETAL(2,I))*DPJ*QPI )
            IF ( DPRINT ) THEN
               WRITE(6,'(A,2(1PE14.5))')"DISPCOEF = ",SDC(1),SDC(2)
            END IF
C---------------------------------------------------
C GET THE RIJ SCALAR
            RIJ = SQRT(
     A           ( COORD(1,I) - COORD(1,J) )**2 +
     B           ( COORD(2,I) - COORD(2,J) )**2 +
     C           ( COORD(3,I) - COORD(3,J) )**2 )/A0
C---------------------------------------------------
C     SCALE THE EFFECTIVE DISTANCE
            RIJ0 = RIJ
            RL = (SDRLD(NAT(I))+SDRLD(NAT(J)))
C/2.0D0
            RH = (SDRHD(NAT(I))+SDRHD(NAT(J)))
C/2.0D0
C     CONVERT FROM ANGSTROM TO BOHR
            RL = RL/A0
            RH = RH/A0
            IF ( RIJ0 .GE. RH ) THEN
               RIJ = RIJ0
            ELSE IF ( RIJ0 < RL ) THEN
               RIJ = (RL+RH)/2.0D0
            ELSE
               RS = (RL-RIJ0)/(RH-RL)
               RS2 = RS*RS
               RS3 = RS2*RS
               RIJ = (RL-RIJ0)*RS3*(2.5D0+THREE*RS+RS2)+PT5*(RL+RH)
            END IF
            RIJV(1) = RIJ
CCC (NOW DO THE SAME THING FOR THE C8 VERSION
            RL = (SDRLQ(NAT(I))+SDRLQ(NAT(J)))
C/2.0D0
            RH = (SDRHQ(NAT(I))+SDRHQ(NAT(J)))
C/2.0D0
C     CONVERT FROM ANGSTROM TO BOHR
            RL = RL/A0
            RH = RH/A0
            IF ( RIJ0 .GE. RH ) THEN
               RIJ = RIJ0
            ELSE IF ( RIJ0 < RL ) THEN
               RIJ = (RL+RH)/2.0D0
            ELSE
               RS = (RL-RIJ0)/(RH-RL)
               RS2 = RS*RS
               RS3 = RS2*RS
               RIJ = (RL-RIJ0)*RS3*(2.5D0+THREE*RS+RS2)+PT5*(RL+RH)
            END IF
            RIJV(2) = RIJ
C-----------------------------------
C     COMPUTE THE DISPERSION ENERGY FOR THIS IJ INTERACTION
            DENE = 0.0D0
            DO 60 K=1,2
               DENE=DENE-SDC(K)/RIJV(K)**(TWO*(K+2))
 60         CONTINUE
            IF ( DPRINT ) THEN
               WRITE(6,'(A,2I3,2E12.4)')"I,J,RIJ,Disp:",
     A              I,J,RIJ0,DENE
            END IF
C-----------------------------------
C     ADD THIS IJ ENERGY TO THE TOTAL ENERGY
            IF ( DPRINT ) THEN
               IF ( I.EQ.1.AND.J.EQ.2 ) THEN
                  WRITE(6,'(A,3F11.5,3E16.8)')"RGREP",
     a                 RIJ0,RIJV(1),RIJV(2),
     a                 -SDC(1)/RIJV(1)**(TWO*(1+2)),
     a                 -SDC(2)/RIJV(2)**(TWO*(2+2)),
     a                 DENE
C     a              EE,
C     a              EE+DENE*EV
               END IF
            END IF
            SDENE = SDENE+DENE
 55      CONTINUE
 50   CONTINUE
      IF ( DPRINT ) WRITE(6,'(A,E20.12)')"   SDENE (AU) = ",SDENE
      SDENE = SDENE * EV
      IF ( DPRINT ) WRITE(6,'(A,E20.12)')"   SDENE (EV) = ",SDENE
C--------------------------------------------
C     ADD THE TOTAL DISPERSION ENERGY TO THE ELECTRONIC ENERGY
      EE = EE + SDENE
      IF ( DPRINT ) WRITE(6,'(A,E20.12)')"EE+SDENE (EV) = ",EE
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     THIS APPEARS IN PELLENQ.src,
C     SO IT SHOULDN'T BE LISTED A SECOND TIME
c$$$      FUNCTION ABOAPB(A,B)
c$$$      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c$$$      C = A+B
c$$$      IF (C.NE.0.0D0) THEN
c$$$         ABOAPB = A*B/C
c$$$      ELSE
c$$$         ABOAPB = 1.D10
c$$$      END IF
c$$$      RETURN
c$$$      END FUNCTION
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE INISD()
C     INITIALIZES THE SIMPLE DISPERSION PARAMETERS
C     BASED ON THE PELLENQ PARAMETERS HARDWIRED
C     INTO BLOCK1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (LMZ=86)
      PARAMETER (LPQ=LMZ+LMZ*(LMZ-1)/2)
      COMMON
     ./SDISPB/ SDRLD(LMZ),SDRLQ(LMZ),SDRHD(LMZ),SDRHQ(LMZ),
     . SDDP(LMZ),SDQP(LMZ),SDDPQ(LMZ),SDQPQ(LMZ),SDNEF(LMZ),NFLSD
     ./PQPAR / PQNEF(LMZ),PQDP(LMZ),PQQP(LMZ),PQOP(LMZ),PQB(LPQ),
     .     PQBO0(LMZ),PQBOZ(LMZ),PQRL(LMZ),PQRH(LMZ),PQDS,PQQS,PQOS
      DO 10 I=1,LMZ
C     EFFECTIVE NUMBER OF ELECTRONS (NEUTRAL ATOM)
         SDNEF(I) = PQNEF(I)
C     DIPOLE POLARIZABILITY (NEUTRAL ATOM)
         SDDP(I)  = PQDP(I)
C     QUADRUPOLE POLARIZABILITY (NEUTRAL ATOM)
         SDQP(I)  = PQQP(I)
C     SCLALED R-LOW USED FOR C6
         SDRLD(I) = PQRL(I)
C     SCALED R-LOW USED FOR C8
         SDRLQ(I) = PQRL(I)
C     SCALED R-HIGH USED FOR C6
         SDRHD(I) = PQRH(I)
C     SCALED R-HIGH USED FOR C8
         SDRHQ(I) = PQRH(I)
C     EXPONENTIAL CHARGE DEPENDENCE OF THE DIPOLE POLARIZABILITY
         SDDPQ(I) = 0.0D0
C     EXPONENTIAL CHARGE DEPENDENCE OF THE QUADRUPOLE POLARIZABILITY
         SDQPQ(I) = 0.0D0
 10   CONTINUE
      END SUBROUTINE
      SUBROUTINE SDGRAD(CG,PA,PB,LM4,LM2)
C     SUBROUTINE CALCULATES THE SIMPLE
C     DISPERSION CONTRIBUTION TO THE GRADIENT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (LM1=600)
      PARAMETER (LMX=9*LM1)
      PARAMETER (LMZ=86)
      DIMENSION CG(3,*)
      DIMENSION PA(LM4),PB(LM4)
      DIMENSION EFFN(LM1),QCHRG(LM1)
      DIMENSION SL(2),ETAL(2,LM1), SDC(2), POLL(2),RIJV(2)
      DIMENSION GRADIJ(3)
      DIMENSION URIJ(3)
C      DIMENSION THKROM(3)
      LOGICAL DPRINT
      COMMON
     ./ATOMC / COORD(3,LM1)
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./PARDER/ TORE(LMZ),EHEAT(LMZ),EISOL(LMZ)
     ./SDISPB/ SDRLD(LMZ),SDRLQ(LMZ),SDRHD(LMZ),SDRHQ(LMZ),
     . SDDP(LMZ),SDQP(LMZ),SDDPQ(LMZ),SDQPQ(LMZ),SDNEF(LMZ),NFLSD
C     DEBUG PRINTING FLAG
C      DPRINT = .TRUE.
      DPRINT = .FALSE.
C     THAKKER OMEGA COEFFICIENTS TO CALCULATE
C     C12, C14, and C16
C     See Thakker_JChemPhys_1988_v89_p2092 equation 29
C      THKROM(1) = 1.0281430356399119D0
C      THKROM(2) = 0.9749145042676423D0
C      THKROM(3) = 0.9735032057454021D0
      IF ( DPRINT ) THEN
         WRITE(6,'(1X)')
         WRITE(6,'(A)')"*******************************"
         WRITE(6,'(A)')"        ENTERED SDGRAD"
         WRITE(6,'(A)')"*******************************"
      END IF
      SDENE = ZERO
C     GET THE EFFECTIVE NUMBER OF ELECTRONS
C     MULPOP APPEARS IN PELLENQ.src
C      CALL MULPOP(PA,PB,LM4,EFFN)
      CALL QPOP(PA,PB,LM4,QCHRG)
      DO 10 I=1,NUMAT
         EFFN(I) = TORE(NAT(I)) - QCHRG(I)
 10   CONTINUE
C      SVNPT5 = THREE*(THREE+TWO)/TWO
      DO 20 I=1,NUMAT
         POLL(1) = SDDP(NAT(I))*EXP(-SDDPQ(NAT(I))*QCHRG(I))
         POLL(2) = SDQP(NAT(I))*EXP(-SDQPQ(NAT(I))*QCHRG(I))
         SDLAM = SDNEF(NAT(I)) / TORE(NAT(I))
         VNE = EFFN(I) * SDLAM
C RARE GAS DIMER HACK
C         VNE = SDNEF(NAT(I))
C----HACK------------
         IF ( VNE .LT. ZERO ) VNE = ZERO
         IF ( DPRINT ) THEN
            WRITE(6,'(A,I3,A,F11.5,A,F11.5,A,F11.5)')" ATOM= ",I
     A           ," LAM= ",SDLAM," NVAL= ",
     B           EFFN(I)," NEFF= ",VNE
         END IF
         SL(1) = VNE
         SL(2) = THREE * SQRT( VNE*POLL(1) )
         DO 25 J=1,2
            IF ( POLL(J).GT.ZERO ) THEN
               ETAL(J,I) = SQRT( SL(J) / POLL(J) )
            ELSE
               ETAL(J,I) = 0.0D0
CSQRT( SL(J) )
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
         DPI = SDDP(NAT(I))
         QPI = SDQP(NAT(I))
         DO 55 J=I+1,NUMAT
            ICNT=ICNT+1
C            WRITE(6,*)"RL D ",SDRLD(NAT(I)),SDRLD(NAT(J))
C            WRITE(6,*)"RL Q ",SDRLQ(NAT(I)),SDRLQ(NAT(J))
C            WRITE(6,*)"RH D ",SDRHD(NAT(I)),SDRHD(NAT(J))
C            WRITE(6,*)"RH Q ",SDRHQ(NAT(I)),SDRHQ(NAT(J))
            DPJ = SDDP(NAT(J))*EXP(-SDDPQ(NAT(J))*QCHRG(J))
            QPJ = SDQP(NAT(J))*EXP(-SDQPQ(NAT(J))*QCHRG(J))
C----------------------------------------------------
C COMPUTE THE DISPERISON COEFFICIENTS
            SDC(1) = PREFA*ABOAPB(ETAL(1,I),ETAL(1,J))*DPI*DPJ
C            WRITE(6,*)"C6 RAW",SDC(1),DPI,DPJ,ETAL(1,I),ETAL(1,J),PREFA
            SDC(2) = PREFB * (
     A           ABOAPB(ETAL(1,I),ETAL(2,J))*DPI*QPJ +
     B           ABOAPB(ETAL(1,J),ETAL(2,I))*DPJ*QPI )
            IF ( DPRINT ) THEN
               WRITE(6,'(A,2(1PE14.5))')"DISPCOEF = ",SDC(1),SDC(2)
            END IF
C---------------------------------------------------
C GET THE RIJ SCALAR
            RIJ = SQRT(
     A           ( COORD(1,I) - COORD(1,J) )**2 +
     B           ( COORD(2,I) - COORD(2,J) )**2 +
     C           ( COORD(3,I) - COORD(3,J) )**2 )/A0
            URIJ(1) = ((COORD(1,I) - COORD(1,J))/A0)/RIJ
            URIJ(2) = ((COORD(2,I) - COORD(2,J))/A0)/RIJ
            URIJ(3) = ((COORD(3,I) - COORD(3,J))/A0)/RIJ
C---------------------------------------------------
C     SCALE THE EFFECTIVE DISTANCE
            RIJ0 = RIJ
            RL = (SDRLD(NAT(I))+SDRLD(NAT(J)))
C/2.0D0
            RH = (SDRHD(NAT(I))+SDRHD(NAT(J)))
C/2.0D0
C     CONVERT FROM ANGSTROM TO BOHR
            RL = RL/A0
            RH = RH/A0
            RL6 = RL
            RH6 = RH
            IF ( RIJ0 .GE. RH ) THEN
               RIJ = RIJ0
            ELSE IF ( RIJ0 < RL ) THEN
               RIJ = (RL+RH)/2.0D0
            ELSE
               RS = (RL-RIJ0)/(RH-RL)
               RS6 = RS
               RS2 = RS*RS
               RS3 = RS2*RS
               RIJ = (RL-RIJ0)*RS3*(2.5D0+THREE*RS+RS2)+PT5*(RL+RH)
            END IF
            RIJV(1) = RIJ
CCC (NOW DO THE SAME THING FOR THE C8 VERSION
            RL = (SDRLQ(NAT(I))+SDRLQ(NAT(J)))
C/2.0D0
            RH = (SDRHQ(NAT(I))+SDRHQ(NAT(J)))
C/2.0D0
C     CONVERT FROM ANGSTROM TO BOHR
            RL = RL/A0
            RH = RH/A0
            RL8 = RL
            RH8 = RH
            IF ( RIJ0 .GE. RH ) THEN
               RIJ = RIJ0
            ELSE IF ( RIJ0 < RL ) THEN
               RIJ = (RL+RH)/2.0D0
            ELSE
               RS = (RL-RIJ0)/(RH-RL)
               RS8=RS
               RS2 = RS*RS
               RS3 = RS2*RS
               RIJ = (RL-RIJ0)*RS3*(2.5D0+THREE*RS+RS2)+PT5*(RL+RH)
            END IF
            RIJV(2) = RIJ
C-----------------------------------
C     COMPUTE THE DISPERSION GRADIENT FOR THIS IJ INTERACTION
C     C6 SECTION -----
            IF ( RIJ0 .LE. RL6 ) THEN
               GRAD6 = 0.D0
            ELSE
               DEDS6 = -6.D0*SDC(1)/(RIJV(1)**7.D0)
               IF ( RIJ0 .GE. RH6 ) THEN
                  DSDR6 = 1.D0
               ELSE
                  DSDR6 = -(RS6**3.D0) *
     A                 (10.D0+3.D0*RS6*(5.D0+2.D0*RS6))
               END IF
               GRAD6 = DEDS6*DSDR6
            END IF
C     C8 SECTION -----
            IF ( RIJ0 .LE. RL8 ) THEN
               GRAD8 = 0.D0
            ELSE
               DEDS8 = -8.D0*SDC(2)/(RIJV(2)**9.D0)
               IF ( RIJ0 .GE. RH8 ) THEN
                  DSDR8 = 1.D0
               ELSE
                  DSDR8 = -(RS8**3.D0) *
     A                 (10.D0+3.D0*RS8*(5.D0+2.D0*RS8))
               END IF
               GRAD8 = DEDS8*DSDR8
            END IF
            GRADIJ(1) = (GRAD6+GRAD8)*URIJ(1)
            GRADIJ(2) = (GRAD6+GRAD8)*URIJ(2)
            GRADIJ(3) = (GRAD6+GRAD8)*URIJ(3)
            DO 1000 K=1,3
C               WRITE(6,'(A,3F15.8)')"CG(K,I)",CG(K,I),
C     A              -GRADIJ(K)*EV*EVCAL/A0,
C     A              CG(K,I) - GRADIJ(K)*EV*EVCAL/A0
               CG(K,I) = CG(K,I) - GRADIJ(K)*EV*EVCAL/A0
               CG(K,J) = CG(K,J) + GRADIJ(K)*EV*EVCAL/A0
 1000       CONTINUE
            IF ( DPRINT ) THEN
               WRITE(6,'(A,2I3,5E15.7)')"I,J,RIJ,S6,S8",I,J,RIJ0,
     A              RIJV(1),RIJV(2),RS6,RS8
               WRITE(6,'(A,5E15.7)')"Grad",
     A              GRADIJ(1)*EV/A0,GRADIJ(2)*EV/A0,GRADIJ(3)*EV/A0,
     A              DSDR6,DSDR8
            END IF
C-----------------------------------
C     ADD THIS IJ ENERGY TO THE TOTAL ENERGY
C            IF ( I.EQ.1.AND.J.EQ.2 ) THEN
c               WRITE(6,'(A,3F11.5,3E16.8)')"RGREP",RIJ0,RIJV(1),RIJV(2),
c     a              -SDC(1)/RIJV(1)**(TWO*(1+2)),
c     a              -SDC(2)/RIJV(2)**(TWO*(2+2)),
c     a              DENE
C             END IF
 55      CONTINUE
 50   CONTINUE
      END SUBROUTINE
