      SUBROUTINE CPEFCN(EE,PA,PB,LM4,LM2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (LM1=600)
      PARAMETER (LMX=9*LM1)
      PARAMETER (LMZ=86)
      INTEGER LM4,LM2
      DIMENSION PA(LM4),PB(LM4)
      LOGICAL UHF
      COMMON
     .     /ATOMC / COORD(3,LM1)
     .     /ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
C.##IF CHARMM
C.##ELSE
     ./INDEX / INDX(LMX)
C.##ENDIF
     .     /CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
     .     /CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     .     /UHF   / UHF
     .     /PARDER/ TORE(LMZ),EHEAT(LMZ),EISOL(LMZ)
     .     /FIELD2/ FIFI(3)
     .     /CPEOPT/ NFLCPE,NPTCPE,NCPEZ
     .     /CPEPAR/ CPEZ(LMZ),CPEED0(LMZ),CPEQ0(LMZ),ATPS(LMZ),
     .     ANPS(LMZ),ATPT(LMZ), ANPT(LMZ),CPESL(LMZ),CPEEM0(LMZ),
     .     CPEM0(LMZ),CPERS(LMZ),CPERR(LMZ),CMXCAT(LMZ),CMXAN(LMZ),
     .     CPEFT(LMZ),CPERC(LMZ),CPEQM(LMZ),CPEPS,CPEZA(LMZ),CPEZB(LMZ)
     .     /CPEOUT/ CPEENE,CPEDIP(3)
c      PARAMETER (LMPAIR=LM1*(LM1-1)/2)
      DIMENSION CPEETA(NUMAT*4,NUMAT*4)
      DIMENSION CPEINV(NUMAT*4,NUMAT*4)
      DIMENSION CPENU(NUMAT*4)
C,CPEMNU(NUMAT*4)
      DIMENSION CPERES(NUMAT*4)
      DIMENSION QREF(NUMAT)
      LOGICAL NOFLD, PTFILE, DOSS
      DIMENSION ETAID(NUMAT*4)
      DIMENSION ETAIN(NUMAT*4)
      DIMENSION CPED(NUMAT*4)
      DIMENSION GLBCON(NUMAT*4,1),GLBC(1)
      DIMENSION DIPCON(6)
      DIMENSION EIGVAL(NUMAT),EIGVEC(NUMAT,NUMAT)
      DIMENSION PMOL(NUMAT,NUMAT)
      DIMENSION CONVEC(NUMAT),WRKVEC(NUMAT)
      DIMENSION WRKMAT(NUMAT,NUMAT)
      DIMENSION ZEROV(NUMAT)
      DIMENSION HDIP(1,NUMAT)
C      DOUBLE PRECISION SCLBND
      DIMENSION SCLBND(NUMAT,NUMAT)
c      DIMENSION CONMT(NUMAT*4,NUMAT*(NUMAT-1)/2+3*NUMAT)
c      INTEGER NP,NC
      NCPE = NUMAT*4
c      NP = NUMAT*(NUMAT-1)/2
c      NC = NP + 3*NUMAT
c      CALL PTV(PA,LM4,1,"PA(IN) ",6,.FALSE.)
c      CALL PTV(PB,LM4,1,"PB(IN) ",6,.FALSE.)
C-----------------------------------------------
C     ARE WE GOING TO DO A SSCPE SOLUTION?
C-----------------------------------------------
      IF ( NFLCPE.EQ.1.OR.NFLCPE.EQ.4 ) THEN
         DOSS = .FALSE.
      ELSE
         DOSS = .TRUE.
      END IF
C-----------------------------------------------
C     SHOULD WE EVEN BE HERE ???
C-----------------------------------------------
      IF(NFLCPE.LT.1) RETURN
      IF(NPTCPE.GT.0) THEN
         WRITE(6,*)
         WRITE(6,'(A)')"----------------------------"
         WRITE(6,'(A)')"       ENTERED CPEFCN"
         WRITE(6,'(A/)')"----------------------------"
         WRITE(6,'(A,3F12.5,A)')"APPLIED FIELD = ",
     1        (FIFI(I),I=1,3)," AU"
      END IF
C-----------------------------------------------
C     DO WE PRINT DETAILED INFO TO A FILE?
C-----------------------------------------------
      ENTRPY = ZERO
      SUM = ZERO
      DO 2 I = 1,3
         SUM = SUM + FIFI(I)**2
 2    CONTINUE
      SUM = SQRT(SUM)
      IF ( SUM .LT. 0.00001D0 .AND. NPTCPE .EQ. -1 ) THEN
         PTFILE = .TRUE.
      ELSE
         PTFILE = .FALSE.
      END IF
      IF ( PTFILE ) OPEN(FILE="MNDO_CPE.DAT",UNIT=58)
C-----------------------------------------------
C     ZERO OUT ARRAYS
C-----------------------------------------------
      CPEEN = ZERO
      CPEENE= ZERO
      CALL CLRV(CPEDIP,3)
      CALL CLRV(DIPCON,6)
      CALL CLRM(3,HDIP,NUMAT)
      CALL CLRM(NCPE,CPEETA,NCPE)
      CALL CLRM(NCPE,CPEINV,NCPE)
      CALL CLRV(CPENU,NCPE)
      CALL CLRV(CPERES,NCPE)
      CALL CLRV(QREF,NUMAT)
c      CALL CLRM(NUMAT,GAMMA,NUMAT)
      CALL CLRM(NCPE,GLBCON,1)
      CALL CLRV(ZEROV,NUMAT)
      CALL CLRM(NUMAT,PMOL,NUMAT)
c      CALL CLRM(NCPE,CONMT,NC)
      GLBC(1) = ZERO
      DO 10 I=1,NCPE,4
         GLBCON(I,1)=ONE
 10   CONTINUE
C-----------------------------------------------
C     CALCULATE THE ARRAYS
C-----------------------------------------------
C     CALCULATE THE MULLIKEN PARTIAL CHARGES
      CALL QPOP(QREF,PA,PB,LM4)
C     CALCULATE BOND SCALING
C      CALL BOSCL(PA,PB,LM4,LM2,NUMAT,SCLBND)
C     CALCULATE THE HYBRID DIPOLE MOMENTS
      CALL HYBDIP(HDIP,PA,PB,LM4,.FALSE.)
C     CALCULATE THE ETA MATRIX AND NU (AKA "M") VECTOR
      CALL CPEDRV(QREF,HDIP,NCPE,PTFILE,CPEETA,CPENU,SCLBND)
C     OUR PROBLEM IS ETA.C = -NU WHICH IS ANALOGOUS TO A.X = B
C     WHERE B = -NU, SO WE NEED TO PASS -NU (CPEMNU)
C      CALL CLRP(CPEETA,CPENU,NUMAT,NCPE)
C      CALL SCLV(-ONE,CPENU,NCPE,CPEMNU)
C-----------------------------------------------
C     SOLVE FOR THE RESPONSE COEFFICIENTS      |
C-----------------------------------------------
      IF ( NPTCPE.GT.0 ) THEN
         WRITE(6,'(/A)')"------------------------------------------"
      END IF
 300  CONTINUE
      IF ( .NOT. DOSS ) THEN
C     *************************************************************
C     GLOBAL NORMALIZATION
C     *************************************************************
C     d
C     -- ( C.NU + 1/2 C.ETA.C - MU(C.D - 0) ) = 0
C     dc
C
C     WHERE D = VOLUME OF BASIS FUNCTION
         IF ( NPTCPE.GT.0 ) THEN
            WRITE(6,'(A)')"SOLVING RESPONSE UNDER GLOBAL NORMALIZATION"
         END IF
C     THE EASY WAY TO DO THIS IS
c     CALL SYMLIN(CPEETA,NCPE,CPEMNU,CPERES,GLBCON,1,GLBC)
C     BUT, I WANT TO REUSE SOME OF THE WORK ARRAYS, SO
C     I WRITE IT OUT BY HAND TO STORE THEM
         CALL CLRV(CPED,NCPE)
         DO 100 I=1,NCPE,4
            CPED(I)=ONE
 100     CONTINUE
         CALL SYMINV(CPEETA,NCPE,CPEINV)
         IF ( NFLCPE.EQ.4 ) THEN
            IF ( NPTCPE.GT.0 ) THEN
               WRITE(6,'(A)')"RESPONSE BASED ON FERMI POPULATIONS"
            END IF
            IF(NPTCPE.GT.0) THEN
               CALL CPEFRM(QREF,CPEINV,CPENU,
     a              CPERES,NUMAT,NCPE,ENTRPY,.TRUE.)
            ELSE
               CALL CPEFRM(QREF,CPEINV,CPENU,
     a              CPERES,NUMAT,NCPE,ENTRPY,.FALSE.)
            END IF
         ELSE
            CALL YSV(NCPE,CPEINV,CPED,ETAID)
            CALL YSV(NCPE,CPEINV,CPENU,ETAIN)
            CALL DOTPRD(CPED,NCPE,ETAID,CONDEN)
            CALL DOTPRD(CPED,NCPE,ETAIN,CONNUM)
            CONSMU = CONNUM/CONDEN
            IF ( NPTCPE.GT.0 ) THEN
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C PRINT ETA INV
c$$$         WRITE(6,*)
c$$$         WRITE(6,*)"----------------------------"
c$$$         WRITE(6,*)"     CPE ETAINV MATRIX"
c$$$         WRITE(6,*)"----------------------------"
c$$$         WRITE(6,'(A/)')"**** SELF ENERGY BLOCKS **** "
c$$$         DO 270 I=1,NUMAT
c$$$            WRITE(6,'(/A,2I3/)')"SELF ENERGY BLOCK: ATOM I/Z",I,NAT(I)
c$$$            IBLK = (I-1)*4
c$$$            WRITE(6,'(15X,4(2(A,I3),A,1X))')" [",I,"|",NAT(I),"]"
c$$$     A           ," [",I,"|",NAT(I),"]"," [",I,"|",NAT(I),"]"
c$$$     B           ," [",I,"|",NAT(I),"]"
c$$$            WRITE(6,'(12X,4A11)')"S   ","P(X)","P(Y)","P(Z)"
c$$$            WRITE(6,'(2(A,I3),A,4(F11.5))')"[",I,"|",NAT(I),"] S   "
c$$$     A           ,(CPEINV(IBLK+1,IBLK+K),K=1,4)
c$$$            WRITE(6,'(2(A,I3),A,4(F11.5))')"[",I,"|",NAT(I),"] P(X)"
c$$$     A           ,(CPEINV(IBLK+2,IBLK+K),K=1,4)
c$$$            WRITE(6,'(2(A,I3),A,4(F11.5))')"[",I,"|",NAT(I),"] P(Y)"
c$$$     A           ,(CPEINV(IBLK+3,IBLK+K),K=1,4)
c$$$            WRITE(6,'(2(A,I3),A,4(F11.5))')"[",I,"|",NAT(I),"] P(Z)"
c$$$     A           ,(CPEINV(IBLK+4,IBLK+K),K=1,4)
c$$$ 270     CONTINUE
c$$$
c$$$         WRITE(6,'(//A/)')"**** COUPLING BLOCKS **** "
c$$$         DO 271 I=1,NUMAT-1
c$$$            IBLK = (I-1)*4
c$$$            DO 272 J=I+1,NUMAT
c$$$               JBLK = (J-1)*4
c$$$               WRITE(6,'(/A,4I3/)')"COUPLING BLOCK: ATOM I/Z,J/Z",
c$$$     A              I,NAT(I),J,NAT(J)
c$$$               WRITE(6,'(15X,4(2(A,I3),A,1X))')" [",J,"|",NAT(J),"]"
c$$$     A              ," [",J,"|",NAT(J),"]"," [",J,"|",NAT(J),"]"
c$$$     B              ," [",J,"|",NAT(J),"]"
c$$$               WRITE(6,'(12X,4A11)')"S   ","P(X)","P(Y)","P(Z)"
c$$$               WRITE(6,'(2(A,I3),A,4(F11.5))')"[",I,"|",NAT(I),"] S   "
c$$$     A              ,(CPEINV(IBLK+1,JBLK+K),K=1,4)
c$$$               WRITE(6,'(2(A,I3),A,4(F11.5))')"[",I,"|",NAT(I),"] P(X)"
c$$$     A              ,(CPEINV(IBLK+2,JBLK+K),K=1,4)
c$$$               WRITE(6,'(2(A,I3),A,4(F11.5))')"[",I,"|",NAT(I),"] P(Y)"
c$$$     A              ,(CPEINV(IBLK+3,JBLK+K),K=1,4)
c$$$               WRITE(6,'(2(A,I3),A,4(F11.5))')"[",I,"|",NAT(I),"] P(Z)"
c$$$     A              ,(CPEINV(IBLK+4,JBLK+K),K=1,4)
c$$$ 272        CONTINUE
c$$$ 271     CONTINUE
C----------------------------------------------------------------------
C----------------------------------------------------------------------
               WRITE(6,*)"LAGRANGE MULTIPLIER = ",CONSMU
            END IF
            CALL ADDVV(CONSMU,ETAID,NCPE,-ONE,ETAIN,CPERES)
         END IF
      ELSE IF ( DOSS ) THEN
C*****************************************************
C     SUBSYSTEM NORMALIZATION
C*****************************************************
C     CHOOSE THE TYPE OF SUBSYSTEMS
         IF ( NPTCPE.GT.0 ) THEN
           WRITE(6,'(A,A)')"SOLVING FOR RESPONSE UNDER ",
     A           "SUBSYSTEM CONSTRAINTS"
        END IF
        IF ( NFLCPE.LE.3 ) THEN
c           CALL CPEMNM(PA,PB,LM4,LM2,EIGVEC,EIGVAL,PMOL)
           CALL CPEMNM(PA,PB,LM4,LM2,PMOL)
         ELSE
            WRITE(6,*)"**********************"
            WRITE(6,*)"FATAL ERROR"
            WRITE(6,*)"FLCPE = ",NFLCPE
            WRITE(6,*)"NO SUBSYSTEM NORMALIZATION IS DEFINED"
            WRITE(6,*)"FOR THIS INPUT OPTION"
            WRITE(6,*)"THIS IS A USER ERROR, NOT A PROGRAMMING ERROR"
            WRITE(6,*)"ABNORMAL TERMINATION"
            STOP
         END IF
c------------------- PRINT SECTION
         IF ( PTFILE ) THEN
            DO 200 I=1,NUMAT
               DO 201 J=1,NUMAT
                  WRITE(58,'(E28.20)')PMOL(I,J)
 201           CONTINUE
 200        CONTINUE
         END IF
c         IF ( NPTCPE.GT.0 ) THEN
c            WRITE(6,'(/A)')"CHARGE TRANSFER PROBABILITY MATRIX"
c            WRITE(6,'(3X,30(I4,3X))')(NAT(I),I=1,NUMAT)
c            DO 202 I=1,NUMAT
c               WRITE(6,'(I3,200F7.3)')NAT(I),(PMOL(I,J),J=1,NUMAT)
c 202        CONTINUE
c         END IF
         CALL PPCPE(CPEETA,CPENU,CPERES,PMOL,NUMAT,NCPE)
C         STOP
C     I NEED ETAID AND ETAIN FOR THE ENERGY CALCULATION
         CALL SYMINV(CPEETA,NCPE,CPEINV)
         CALL YSV(NCPE,CPEINV,CPED,ETAID)
         CALL YSV(NCPE,CPEINV,CPENU,ETAIN)
C         CALL SYMLIN(CPEETA,NCPE,CPEMNU,CPERES,CONMT,NUMAT,ZEROV)
C         CALL SYMINV(CPEETA,NCPE,CPEINV)
C     SOLVE FOR THE RESPONSE COEFFICIENTS
c         CALL SSCPE(CPERES,CPEINV,CPENU,EIGVEC,EIGVAL,ICON,LCON)
C         CALL SSCPE(CPERES,CPEINV,CPENU,PMOL)
      END IF
      IF(NPTCPE.GT.0) THEN
         WRITE(6,*)
         WRITE(6,*)"----------------------------"
         WRITE(6,*)" CPE RESPONSE COEFFICIENTS"
         WRITE(6,*)"----------------------------"
         WRITE(6,'(2A3,4A12)')"I","Z","S","P(X)","P(Y)","P(Z)"
         DO 358 I=1,NUMAT,1
            J = (I-1)*4
            WRITE(6,'(2I3,4F12.5:)')I,NAT(I),(CPERES(J+K),K=1,4)
 358     CONTINUE
      END IF
C     ------------------------------------------------------------
C     CALCULATE THE ENERGY
C     ------------------------------------------------------------
C     NOW WE HAVE ALL THE INFORMATION WE NEED TO CALCULATE THE ENERGY
C     CPEE1 IS THE FIRST TERM INVOLVING A DOT PRODUCT OF NU AND C
      CPEE1 = ZERO
      DO 400 I=1,NCPE
         CPEE1 = CPEE1 + CPERES(I)*CPENU(I)
 400  CONTINUE
C     CPEE2 IS THE SECOND TERM INVOLVING C.ETA.C
      CPEE2 = ZERO
C     I START WITH PERFORMING ETA.C AND STORING IT IN CPENU
      CALL DSYMV("U",NCPE,ONE,CPEETA,NCPE,CPERES,1,ZERO,CPENU,1)
C     I THEN TAKE THE DOT PRODUCT OF C.(ETA.C), BUT (ETA.C) IS
C     STORED IN NU, SO IT LOOKS LIKE C.NU
      DO 410 I=1,NCPE
         CPEE2 = CPEE2 + CPERES(I)*CPENU(I)
 410  CONTINUE
C     THE TOTAL ENERGY IS BELOW.  DON'T FORGET THE 0.5 PREFACTOR
C                       (WHICH SHOWS UP AS PT5)
      CPEEN = CPEE1 + PT5*CPEE2
C      CPEEN = CPEE1
      CPEEN = CPEEN - ENTRPY
      IF(NPTCPE.GT.0) THEN
         WRITE(6,*)
         WRITE(6,*)"**************************************************",
     A     "*********************"
         WRITE(6,'(A,E12.5,A)')"  CPE CONTRIBUTION TO THE ENERGY = "
     A        ,CPEEN," HARTREE"
         WRITE(6,'(2A,E12.5,A)')"  ENTROPY CONTRIBUTION TO THE ",
     A        "CPE ENERGY (-T*S) = "
     A        ,-ENTRPY," HARTREE"
         WRITE(6,'(A,E12.5,A)')"  TOTAL ENERGY       = ",EE," EV"
      END IF
      CPEENE = CPEEN*EV
      EE = EE + CPEENE
      IF(NPTCPE.GT.0) THEN
         WRITE(6,'(A,E12.5,A)')"  TOTAL ENERGY + CPE = ",EE," EV"
      END IF
C     ----------------------------------------------------------------
C     CALCULATE THE DIPOLE MOMENT CAUSED BY THE CPE RESPONSE
C     ----------------------------------------------------------------
      DO 500 ICRD=1,3
         CPEDIP(ICRD) = ZERO
 500  CONTINUE
      IF ( NPTCPE.GT.0 ) THEN
         WRITE(6,*)"**************************************************",
     A     "*********************"
         WRITE(6,*)"--------------------------------------------------",
     A        "---------------------"
         WRITE(6,'(A)')" CPE CONTRIBUTION TO THE DIPOLE MOMENT (D)"
         WRITE(6,'(A,3F11.4)')" FIELD STRENGTH = ",(FIFI(ICRD),ICRD=1,3)
         WRITE(6,*)"--------------------------------------------------",
     A        "---------------------"
         WRITE(6,'(2A3,6A11)')"I","Z","S(X)","S(Y)","S(Z)",
     A        "P(X)","P(Y)","P(Z)"
      END IF
      DO 510 IL=1,NUMAT
         IBLK = (IL-1)*4+1
C     THE MINUS SIGNS ARE BECAUSE OF THE DIFFERENCE BETWEEN
C     RESPONSE COEFFICIENTS OF DENSITY VERSUS RESPONSE COEFFICIENTS
C     OF CHARGE!!!
         IF ( NPTCPE.GT.0 ) THEN
            WRITE(6,'(2I3,6F11.4)')IL,NAT(IL),
     A           CPERES(IBLK)*COORD(1,IL)*4.80324D0,
     B           CPERES(IBLK)*COORD(2,IL)*4.80324D0,
     C           CPERES(IBLK)*COORD(3,IL)*4.80324D0,
     D           (CPERES(IBLK+ICRD)*A0*4.80324D0,ICRD=1,3)
            DO 509 ICRD=1,3
               DIPCON(ICRD) = DIPCON(ICRD) +
     A              CPERES(IBLK)*COORD(ICRD,IL)*4.80324D0
               DIPCON(ICRD+3) = DIPCON(ICRD+3) +
     A              CPERES(IBLK+ICRD)*A0*4.80324D0
 509        CONTINUE
         END IF
         DO 520 ICRD=1,3
C     MONOPOLE CONTRIBUTION
            CPEDIP(ICRD) = CPEDIP(ICRD)-CPERES(IBLK)*COORD(ICRD,IL)/A0
C     DIPOLE CONTRIBUTION
            CPEDIP(ICRD) = CPEDIP(ICRD)-CPERES(IBLK+ICRD)
 520     CONTINUE
 510  CONTINUE
C     CONVERT FROM ELECTRON-BOHR TO FR (ESU)
      DO 530 ICRD=1,3
         CPEDIP(ICRD) = -CPEDIP(ICRD)*A0*4.80324D0
 530  CONTINUE
      IF(NPTCPE.GT.0) THEN
         WRITE(6,*)"--------------------------------------------------",
     A        "---------------------"
         WRITE(6,'(A6,6F11.4)')" TOTAL",(DIPCON(J),J=1,6)
         WRITE(6,'(A6,3F11.4)')" FINAL",(CPEDIP(J),J=1,3)
         WRITE(6,*)"--------------------------------------------------",
     A        "---------------------"
         FMAG = ZERO
         DO 610 ICRD=1,3
            FAMG = FMAG + FIFI(ICRD)**2
 610     CONTINUE
C         IF ( NPTCPE .GT. 0 ) THEN
C         IF ( FMAG .EQ. ZERO ) THEN
C            DO 620 I=1,NUMAT
C               ICPE = (I-1)*4+1
C               WRITE(6,'(A,2I3,F12.5)')"CPE MONOPOLE RESPONSE ",
C     A              I,NAT(I),CPERES(ICPE)
C 620        CONTINUE
C         END IF
C     END IF
         WRITE(6,*)"NET MULLIKEN CHARGES"
         WRITE(6,*)"--------------------------------------------------",
     A        "---------------------"
         WRITE(6,'(2A3,3A11)')"I","Z","TOTAL","CPE","AB INITIO"
         DO 600 I=1,NUMAT
            WRITE(6,'(2I3,3F11.4)')I,NAT(I),QREF(I)+CPERES((I-1)*4+1),
     A           CPERES((I-1)*4+1),QREF(I)
 600     CONTINUE
         WRITE(6,*)"**************************************************",
     A     "*********************"
         WRITE(6,*)
      END IF
      IF ( PTFILE ) CLOSE(UNIT=58)
      RETURN
      END SUBROUTINE CPEFCN
      SUBROUTINE CPEDRV(QREF,HDIP,NCPE,PTFILE,CPEETA,CPENU,SCLBND)
C*************************************************
C     CALCULATES THE ETA MATRIX AND NU VECTOR
C*************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (LM1=600)
      PARAMETER (LMZ=86)
      COMMON
     .     /ATOMC / COORD(3,LM1)
     .     /ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     .     /CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
     .     /CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     .     /FIELD2/ FIFI(3)
     .     /CPEOPT/ NFLCPE,NPTCPE,NCPEZ
     .     /CPEPAR/ CPEZ(LMZ),CPEED0(LMZ),CPEQ0(LMZ),ATPS(LMZ),
     .     ANPS(LMZ),ATPT(LMZ), ANPT(LMZ),CPESL(LMZ),CPEEM0(LMZ),
     .     CPEM0(LMZ),CPERS(LMZ),CPERR(LMZ),CMXCAT(LMZ),CMXAN(LMZ),
     .     CPEFT(LMZ),CPERC(LMZ),CPEQM(LMZ),CPEPS,CPEZA(LMZ),CPEZB(LMZ)
     .     /CPEOUT/ CPEENE,CPEDIP(3)
      PARAMETER (LMPAIR=LM1*(LM1-1)/2)
      DIMENSION CPEETA(NUMAT*4,NUMAT*4)
      DIMENSION CPENU(NUMAT*4)
      DIMENSION CPERES(NUMAT*4)
      DIMENSION CRDS(3)
      DIMENSION QREF(NUMAT)
      DIMENSION HDIP(3,NUMAT)
      LOGICAL PTFILE
C      DOUBLE PRECISION SCLBND
      DIMENSION SCLBND(NUMAT,NUMAT)
      IF ( NPTCPE.GT.0 ) THEN
         WRITE(6,'(2A3,5A15)')"I","Z","Qmull","Zeta","Dx","Dy","Dz"
      END IF
C     -----------------------------------------------------------
C     CONSTRUCT THE COULOMB (ETA) CPE MATRIX
C     AS WE DO THIS, WE ARE ALSO GOING TO CALCULATE
C     THE NU VECTOR - WHICH IS THE COULOMB INTERACTION
C     OF THE REFERENCE POTENTIAL WITH THE RESPONSE BASIS
C     THE GAUSSIAN IS (ZI*ZI/PI)**1.5 * EXP(ZI**2 * ABS(R-RI)**2)
C     -----------------------------------------------------------
C     HERE ARE SOME CONSTANTS THAT ARE NEEDED FOR THE SELF
C     GAUSSIAN COULOMB ENERGY
      PT333 = ONE/THREE
      SQ2OPI = SQRT(TWO/PI)
      LOGSL  = LOG(ONE/(FOUR*(TWO+THREE)))
      TSQOP  = (TWO/SQRT(PI))
      FTSQPI = FOUR/(THREE*SQRT(PI))
C     -----------------------------------------------------------
C     ILOOP------------------------------------------------------
      DO 100 IL=1,NUMAT,1
         NATIL = NAT(IL)
C     THIS IS THE (INDEX-1) OF THE FIRST CARTESIAN DIPOLE OF ATOM I
         IBLK = (IL-1)*4+1
C     GET THE ZETA EXPONENT OF ATOM I
C         CALL CPEZET(IL,ZIOLD,QREF(IL))
         CALL NEWZET(IL,ZI,QREF(IL))
C         WRITE(6,*)"NEW/OLD ZETAS",ZI,ZIOLD
         ZI3 = ZI**3
C     THIS IS THE SELF ENERGY OF THE GAUSSIAN DIPOLE
         SELFED = PT333*ZI3*SQ2OPI
C     THIS IS THE SELF ENERGY OF THE GAUSSIAN MONOPOLE
         SELFEM = ZI*SQ2OPI
         IF ( NPTCPE.GT.0 ) THEN
            WRITE(6,'(2I3,5F15.9)')IL,NAT(IL),QREF(IL),ZI,
     1           HDIP(1,IL),HDIP(2,IL),HDIP(3,IL)
         END IF
C     -----------------------------------------
C     JLOOP------------------------------------
         DO 200 JL=1,NUMAT,1
            NATJL = NAT(JL)
            JBLK = (JL-1)*4+1
C     GET THE ZETA EXPONENT OF ATOM J
C            CALL CPEZET(JL,ZJ,QREF(JL))
            CALL NEWZET(JL,ZJ,QREF(JL))
            IF ( IL .EQ. JL ) THEN
C----------------------------------------------------
C      DIAGONAL BLOCK OF THE HARDNESS MATRIX
C----------------------------------------------------
C     DIAGONAL ELEMENT, SELF ENERGY CALCULATION
C     MONOPOLE SELF ENERGY
               CPEETA(IBLK,JBLK) = SELFEM
               DO 105 ICRD = 1,3
C     ZERO OUT THE MONOPOLE-DIPOLE CROSS TERMS
                  CPEETA(IBLK+ICRD,JBLK) = ZERO
                  CPEETA(IBLK,JBLK+ICRD) = ZERO
                  DO 106 JCRD = 1,3
                     CPEETA(IBLK+ICRD,JBLK+JCRD) = ZERO
 106              CONTINUE
C     (1/3)*SQRT(2/PI)*ZI**3
C     DIPOLE-DIPOLE SELF ENERGY
                  CPEETA(IBLK+ICRD,JBLK+ICRD) = SELFED
 105           CONTINUE
C     THE REFERENCE POTENTIAL OF THE JL'TH ATOM
C     MAKES NO CONTRIBUTION TO THE NU VECTOR DUE
C     TO SYMMETRY... THE DIPOLE RESPONSE BASIS
C     LOCATED ON THE IL'TH ATOM IS LOCATED ON
C     THE SAME CENTER AS THE JL'TH REFERENCE POTENTIAL
C     NOTE THAT THE REFERENCE POTENTIAL IS SIMPLY
C     TREATED AS THE SUM OF POINT CHARGES LOCATED
C     ON EACH OF THE NUCELEI
C     ***HOWEVER**** THE MONOPOLE RESPONSE FUNCTION
C     INTERACTS WITH THE POINT CHARGE ON TOP OF IT
C     IF OUR ZETA IS TOO BIG, THEN THIS IS INFINITY, SO
C     WE SET IT EQUAL TO ZERO IF IT IS A SELF POINT CHARGE INTERACTION
               IF ( ZI .LT. 10.0D0 ) THEN
C                  CPENU(IBLK)=CPENU(IBLK) + CPERS(NAT(IL)) *
C     A                 QREF(IL)*CPEQ0(NAT(IL))*TSQOP*ZI
                  CPENU(IBLK)=CPENU(IBLK)+QREF(IL)*TSQOP*ZI
               END IF
C     ADD IN THE INTERACTION OF REFERENCE DIPOLE JL WITH
C     RESPONSE DIPOLE IL BOTH LOCATED ON THE SAME CENTER
c$$$               IF ( ZI .LT. 10.0D0 ) THEN
c$$$c$$$                  DO 108 ICRD=1,3
c$$$c$$$                     CPENU(IBLK+ICRD)=CPENU(IBLK+ICRD)+
c$$$c$$$     A                    CPERS(NAT(IL)) *
c$$$c$$$     B                    HDIP(ICRD,JL)*CPEQ0(NAT(IL))*
c$$$c$$$     C                    FTSQPI*ZI**3
c$$$c$$$ 108              CONTINUE
c$$$                  DO 108 ICRD=1,3
c$$$                     CPENU(IBLK+ICRD)=CPENU(IBLK+ICRD)+
c$$$     A                    CPERS(NAT(IL)) *
c$$$     B                    HDIP(ICRD,JL)*CPEQ0(NAT(IL))*
c$$$     C                    SELFED
c$$$ 108              CONTINUE
c$$$               END IF
            ELSE
C --------------------------------------------------
C     OFF DIAGONAL BLOCK OF THE HARDNESS MATRIX
C       WE CONSTRUCT THE NU VECTOR HERE TOO
C --------------------------------------------------
C     THESE DISTANCES ARE NOW BOHR DUE TO THE A0 CONVERSTION FACTOR
               CRDS(1) = (COORD(1,IL)-COORD(1,JL))/A0
               CRDS(2) = (COORD(2,IL)-COORD(2,JL))/A0
               CRDS(3) = (COORD(3,IL)-COORD(3,JL))/A0
               RIJ  = SQRT(CRDS(1)**2+CRDS(2)**2+CRDS(3)**2)
C     THIS DAMPS THE INTERACTION BY SOME NUMBER BETWEEN
C     0 AND 1
               DAMPFN = ZERO
               CALL CPEDMP(DAMPFN,RIJ,NATIL,NATJL)
C********************************************************
C     ---------------------------------------------
C     ---- ETA MATRIX OFF DIAGONAL SECTION ----
C     ---------------------------------------------
C********************************************************
C********************************************************
C     GAUSSIAN DIPOLE - GAUSSIAN DIPOLE SUBBLOCK
C********************************************************
               ZIJ    = ZI*ZJ/SQRT(ZI**2+ZJ**2)
               RIJ5   = RIJ**5
               RIJ3   = RIJ**3
               ZR     = ZIJ * RIJ
               ZR2    = ZR*ZR
               ERFZR  = ERF(ZR)
               ZEXR   = TSQOP*ZR*EXP(-ZR2)
               TZR    = THREE+TWO*ZR2
               C11OFF = (TZR*ZEXR-THREE*ERFZR)/RIJ5
               C11ON  = (ERFZR-ZEXR)/RIJ3
C PLACE IN THE EFFECTIVE DIELECTRIC PREFACTOR
               C11OFF = C11OFF*CPEPS
               C11ON  = C11ON*CPEPS
C     ICRD AND JCRD LOOP OVER THE X Y Z CARTESIAN INDEX
C     THESE LOOPS PERFORM A CROSS PRODUCT
               DO 115 ICRD = 1,3
                  DO 120 JCRD = 1,3
                     CPEETA(IBLK+ICRD,JBLK+JCRD) =
     A                    C11OFF*CRDS(ICRD)*CRDS(JCRD)
 120              CONTINUE
C     THE DIPOLE-DIPOLE INTERACTION HAS A KRONECKER DELTA FUNCTION
C     WHICH APPEARS AND IS 1 WHEN THE CARTESIAN INDICIES ARE THE SAME
C     SO THIS STUFF IS OUTSIDE THE JCRD LOOP AND ICRD IS USED
C     IN BOTH ARRAY DIMENSIONS
                  CPEETA(IBLK+ICRD,JBLK+ICRD) =
     A                 CPEETA(IBLK+ICRD,JBLK+ICRD) + C11ON
 115           CONTINUE
C***********************************************************
C   END GAUSSIAN DIPOLE - GAUSSIAN DIPOLE SUBBLOCK
C***********************************************************
C********************************************************
C     GAUSSIAN DIPOLE - GAUSSIAN MONOPOLE SUBBLOCK
C                         AND
C     GAUSSIAN MONOPOLE - GAUSSIAN DIPOLE SUBBLOCK
C********************************************************
               DO 123 JCRD = 1,3
                  CPEETA(IBLK,JBLK+JCRD)= -CRDS(JCRD)*C11ON
                  CPEETA(IBLK+JCRD,JBLK)=  CRDS(JCRD)*C11ON
 123           CONTINUE
C********************************************************
C  END DIPOLE-MONOPOLE SUBBLOCK
C********************************************************
C********************************************************
C  RESPONSE MONOPOLE- RESPONSE MONOPOLE SUBBLOCK
C********************************************************
               CPEETA(IBLK,JBLK) = CPEPS*DAMPFN* ERFZR/RIJ
C********************************************************
C  END MONOPOLE-MONOPOLE SUBBLOCK
C********************************************************
C     --- END ETA MATRIX SECTION ---
C     --- NOW WE CALCULATE THE NU VECTOR ---
C     THE ILOOP GOES OVER RESPONSE
C     THE JLOOP GOES OVER ATOMIC POINT CHARGES
C********************************************************
C  RESPONSE DIPOLE- REFERENCE POINT SUBBLOCK
C********************************************************
C     THIS IS THE COULOMB INTERACTION OF THE RESPONSE DIPOLE ICRD
C     WITH THE PARTIAL CHARGES OF THE JL ATOMS
               ZIR    = ZI*RIJ
               ERFZIR = ERF(ZIR)
               CPENUT = TSQOP*ZIR*EXP(-ZIR*ZIR)-ERFZIR
               CPENUT = CPENUT/RIJ3
C               CALL CPEDMP(DAMPFN,RIJ,CPESL(NAT(IL)),CPEQ0(NAT(JL)),
C     a              CPERC(NAT(IL)))
               QEFF   =
C     THE NEXT LINE WORKS, BUT CHANGES THE POLARIZABILITY
c     a              QREF(JL) + RIJ*(CPEQM(NAT(JL))-CPEQM(NAT(IL)))/2
     a              CPEPS * QREF(JL) * DAMPFN
C     b              + DAMPFN*(CPEQM(NAT(JL))-CPEQM(NAT(IL)))
C               WRITE(6,'(A,2I3,F12.5)')"DAMP",IL,JL,DAMPFN
               CPENUT = CPENUT*QEFF
               DO 150 ICRD = 1,3
                  CPENU(IBLK+ICRD)=CPENU(IBLK+ICRD)+
C SCLBND(IL,JL)*
     a                 CPENUT*CRDS(ICRD)
 150           CONTINUE
C**********************************************************
C  RESPONSE DIPOLE-CHEMICAL POTENTIAL FIELD SUBBLOCK
C**********************************************************
C               ZZ = ZI*ZI*ZJ*ZJ/(ZI*ZI+ZJ*ZJ)
c               ZZ = CPERS(NATIL)**2 * CPERS(NATJL)**2 /
c     a              ( CPERS(NATIL)**2 + CPERS(NATJL)**2 )
c               FLDMAG =
c     a              (ONE-SCLBND(IL,JL))*
c     a              (CPEQM(NATJL)-CPEQM(NATIL)) *
c     b              (ZZ/PI)**1.5D0 * EXP(-ZZ*RIJ)
c               DO 149 ICRD=1,3
c                  CPENU(IBLK+ICRD)=CPENU(IBLK+ICRD)+FLDMAG*
c     a                 CRDS(ICRD)/RIJ
c 149           CONTINUE
c$$$               DO 149 ICRD=1,3
c$$$                  CPENU(IBLK+ICRD)=CPENU(IBLK+ICRD)+
c$$$     c                 CRDS(ICRD)*(CPEQM(NAT(IL))-CPEQM(NAT(JL)))
c$$$     d                 / (2*RIJ*RIJ)
c$$$ 149           CONTINUE
C********************************************************
C  END DIPOLE-POINT SUBBLOCK
C********************************************************
C********************************************************
C  RESPONSE MONOPOLE- REFERENCE POINT SUBBLOCK
C********************************************************
               CPENU(IBLK) = CPENU(IBLK)+
cSCLBND(IL,JL)*
     a              QEFF*ERFZIR/RIJ
C********************************************************
C  END MONOPOLE-POINT SUBBLOCK
C********************************************************
C********************************************************
C  RESPONSE DIPOLE-REFERENCE POINT DIPOLE SUBBLOCK
C********************************************************
               ZIJ    = ZI
               RIJ5   = RIJ**5
               RIJ3   = RIJ**3
               ZR     = ZIJ * RIJ
               ZR2    = ZR*ZR
               ERFZR  = ERF(ZR)
               ZEXR   = TSQOP*ZR*EXP(-ZR2)
               TZR    = THREE+TWO*ZR2
               C11OFF = CPEPS*DAMPFN*(TZR*ZEXR-THREE*ERFZR)/RIJ5
               C11ON  = CPEPS*DAMPFN*(ERFZR-ZEXR)/RIJ3
C     ICRD AND JCRD LOOP OVER THE X Y Z CARTESIAN INDEX
C     THESE LOOPS PERFORM A CROSS PRODUCT
               DO 151 ICRD = 1,3
                  DO 152 JCRD = 1,3
                     CPENU(IBLK+ICRD) = CPENU(IBLK+ICRD) +
     B                    HDIP(JCRD,JL) *
     A                    C11OFF*CRDS(ICRD)*CRDS(JCRD)
c*
c     C                    SCLBND(IL,JL)
 152              CONTINUE
C     THE DIPOLE-DIPOLE INTERACTION HAS A KRONECKER DELTA FUNCTION
C     WHICH APPEARS AND IS 1 WHEN THE CARTESIAN INDICIES ARE THE SAME
C     SO THIS STUFF IS OUTSIDE THE JCRD LOOP AND ICRD IS USED
C     IN BOTH ARRAY DIMENSIONS
                  CPENU(IBLK+ICRD) =
     A                 CPENU(IBLK+ICRD) +
     B                 HDIP(ICRD,JL) *C11ON
c*
c     C                 SCLBND(IL,JL)
 151           CONTINUE
C********************************************************
C  END DIPOLE-POINT DIPOLE SUBBLOCK
C********************************************************
C********************************************************
C  RESPONSE MONOPOLE-REFERENCE POINT DIPOLE SUBBLOCK
C********************************************************
               DO 153 JCRD = 1,3
                  CPENU(IBLK)= CPENU(IBLK) +
     a                 HDIP(JCRD,JL) *
     a                 CRDS(JCRD)*C11ON
c*
c     a                 SCLBND(IL,JL)
 153           CONTINUE
C********************************************************
C  END GAUSSIAN MONOPOLE-POINT DIPOLE SUBBLOCK
C********************************************************
C     WE ARE NOT DONE WITH THE NU VECTOR YET!
C     WE NEED TO GET OUT OF THIS JL LOOP AND THEN
C     ADD IN THE EXTERNAL APPLIED FIELD CONTRIBUTION
C     --- END NU VECTOR SECTION
            END IF
C     END IF IL == IJ -----------------
 200     CONTINUE
C     END JLOOP--------------------------------
C     -----------------------------------------
C     OK.  WE NEED TO FINISH UP THE NU VECTOR BY ADDING IN
C     THE CONTRIBUTIONS DUE TO APPLIED ELECTRIC FIELDS.
C     FIFI SHOULD ALREADY BY IN ATOMIC UNITS.
         DO 155 ICRD = 1,3
C     MONOPOLE GAUSSIAN WITH THE FIELD IS THE FIELD*COORDINATE
C     THE MINUS SIGN IN FRONT OF THE FIELD HAS TO DO WITH
C     THE ISSUE OF WORKING IN TERMS OF CHARGE VERSUS DENSITY
            CPENU(IBLK)=CPENU(IBLK)-CPEPS *
     1           FIFI(ICRD)*COORD(ICRD,IL)/A0
C     DIPOLE GAUSSIAN WITH THE FIELD IS JUST THE FIELD
C     THE MINUS SIGN IN FRONT OF THE FIELD HAS TO DO WITH
C     THE ISSUE OF WORKING IN TERMS OF CHARGE VERSUS DENSITY
            CPENU(IBLK+ICRD)=CPENU(IBLK+ICRD)-CPEPS*FIFI(ICRD)
 155     CONTINUE
C     AND ADD IN THE MU0 PARAMETER TO THE MONOPOLE FUNCTIONS
         CPENU(IBLK)=CPENU(IBLK)+CPEM0(NATIL)
C     ADD THE ETA0 PARAMETERS TO THE DIAGONAL ELEMENTS
C     OF THE HARDNESS MATRIX
C     HERE IS THE MONOPOLE PARAMETER
         CPEETA(IBLK,IBLK) = CPEETA(IBLK,IBLK)+CPEEM0(NATIL)
C     HERE ARE THE DIPOLE PARAMETERS
         DO 160 ICRD = 1,3
            CPEETA(IBLK+ICRD,IBLK+ICRD) =
     1          CPEETA(IBLK+ICRD,IBLK+ICRD)+CPEED0(NATIL)
 160     CONTINUE
 100  CONTINUE
C     END ILOOP--------------------------------------------
C     -----------------------------------------------------
c$$$C     NOW WE ARE GOING TO SCALE THE COUPLING SUBBLOCKS
c$$$C     BY THE RESPONSE-RESPONSE PARAMETER
c$$$
c$$$      DO 250 I=1,NUMAT,1
c$$$         IBLK = (I-1)*4
c$$$         RRI  = CPERR(NAT(I))
c$$$         DO 260 J=1,NUMAT,1
c$$$            JBLK = (J-1)*4
c$$$            IF ( J.EQ.I ) THEN
c$$$               RRJ = ONE
c$$$               RRI = ONE
c$$$            ELSE
c$$$               RRJ = CPERR(NAT(J))
c$$$            END IF
c$$$            DO 255 K=1,4,1
c$$$               DO 256 L=1,4,1
c$$$                  CPEETA(IBLK+K,JBLK+L)=RRI*RRJ*CPEETA(IBLK+K,JBLK+L)
c$$$ 256           CONTINUE
c$$$ 255        CONTINUE
c$$$ 260     CONTINUE
c$$$ 250  CONTINUE
      IF(NPTCPE.GT.0) THEN
C     PRINT THE ETA MATRIX
         WRITE(6,*)
         WRITE(6,*)"----------------------------"
         WRITE(6,*)"       CPE ETA MATRIX"
         WRITE(6,*)"----------------------------"
         WRITE(6,'(A/)')"**** SELF ENERGY BLOCKS **** "
         DO 270 I=1,NUMAT
            WRITE(6,'(/A,2I3/)')"SELF ENERGY BLOCK: ATOM I/Z",I,NAT(I)
            IBLK = (I-1)*4
            WRITE(6,'(15X,4(2(A,I3),A,1X))')" [",I,"|",NAT(I),"]"
     A           ," [",I,"|",NAT(I),"]"," [",I,"|",NAT(I),"]"
     B           ," [",I,"|",NAT(I),"]"
            WRITE(6,'(12X,4A11)')"S   ","P(X)","P(Y)","P(Z)"
            WRITE(6,'(2(A,I3),A,4(F11.5))')"[",I,"|",NAT(I),"] S   "
     A           ,(CPEETA(IBLK+1,IBLK+K),K=1,4)
            WRITE(6,'(2(A,I3),A,4(F11.5))')"[",I,"|",NAT(I),"] P(X)"
     A           ,(CPEETA(IBLK+2,IBLK+K),K=1,4)
            WRITE(6,'(2(A,I3),A,4(F11.5))')"[",I,"|",NAT(I),"] P(Y)"
     A           ,(CPEETA(IBLK+3,IBLK+K),K=1,4)
            WRITE(6,'(2(A,I3),A,4(F11.5))')"[",I,"|",NAT(I),"] P(Z)"
     A           ,(CPEETA(IBLK+4,IBLK+K),K=1,4)
 270     CONTINUE
         WRITE(6,'(//A/)')"**** COUPLING BLOCKS **** "
         DO 271 I=1,NUMAT-1
            IBLK = (I-1)*4
            DO 272 J=I+1,NUMAT
               JBLK = (J-1)*4
               WRITE(6,'(/A,4I3/)')"COUPLING BLOCK: ATOM I/Z,J/Z",
     A              I,NAT(I),J,NAT(J)
               WRITE(6,'(15X,4(2(A,I3),A,1X))')" [",J,"|",NAT(J),"]"
     A              ," [",J,"|",NAT(J),"]"," [",J,"|",NAT(J),"]"
     B              ," [",J,"|",NAT(J),"]"
               WRITE(6,'(12X,4A11)')"S   ","P(X)","P(Y)","P(Z)"
               WRITE(6,'(2(A,I3),A,4(F11.5))')"[",I,"|",NAT(I),"] S   "
     A              ,(CPEETA(IBLK+1,JBLK+K),K=1,4)
               WRITE(6,'(2(A,I3),A,4(F11.5))')"[",I,"|",NAT(I),"] P(X)"
     A              ,(CPEETA(IBLK+2,JBLK+K),K=1,4)
               WRITE(6,'(2(A,I3),A,4(F11.5))')"[",I,"|",NAT(I),"] P(Y)"
     A              ,(CPEETA(IBLK+3,JBLK+K),K=1,4)
               WRITE(6,'(2(A,I3),A,4(F11.5))')"[",I,"|",NAT(I),"] P(Z)"
     A              ,(CPEETA(IBLK+4,JBLK+K),K=1,4)
 272        CONTINUE
 271     CONTINUE
         WRITE(6,*)
         WRITE(6,*)"----------------------------"
         WRITE(6,*)"       CPE NU VECTOR"
         WRITE(6,*)"----------------------------"
         WRITE(6,'(2A3,4A12)')"I","Z","S","P(X)","P(Y)","P(Z)"
         DO 290 I=1,NUMAT,1
            J = (I-1)*4
            WRITE(6,'(2I3,4F12.5:)')I,NAT(I),(CPENU(J+K),K=1,4)
 290     CONTINUE
      END IF
      IF ( PTFILE ) THEN
         WRITE(58,*)NUMAT,4
         DO 292 I=1,NCPE
            WRITE(58,'(E28.20)')CPENU(I)
 292     CONTINUE
         DO 293 I=1,NCPE
            DO 294 J=1,NCPE
               WRITE(58,'(E28.20)')CPEETA(I,J)
 294        CONTINUE
 293     CONTINUE
      END IF
      END SUBROUTINE
C***********************************************************************
C***********************************************************************
      SUBROUTINE CPEZET(IA,ZA,QA)
C     --------------------------------------------------------
C     CALCULATE THE ZETA EXPONENT OF THE CPE RESPONSE FUNCTION
C     --------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (LM1=600)
      PARAMETER (LMX=9*LM1)
      PARAMETER (LMZ=86)
      COMMON
     .     /ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     .     /CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
     .     /CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     .     /CPEOPT/ NFLCPE,NPTCPE,NCPEZ
     .     /CPEPAR/ CPEZ(LMZ),CPEED0(LMZ),CPEQ0(LMZ),ATPS(LMZ),
     .     ANPS(LMZ),ATPT(LMZ), ANPT(LMZ),CPESL(LMZ),CPEEM0(LMZ),
     .     CPEM0(LMZ),CPERS(LMZ),CPERR(LMZ),CMXCAT(LMZ),CMXAN(LMZ),
     .     CPEFT(LMZ),CPERC(LMZ),CPEQM(LMZ),CPEPS,CPEZA(LMZ),CPEZB(LMZ)
C     THESE ARE JUST SOME CONSTANTS
      PT333 = ONE/THREE
      SQ2OPI = SQRT(TWO/PI)
      PREFZ = (THREE/SQ2OPI)**PT333
      CUBANG = 1/A0
      CUBANG = CUBANG**3
      IF ( NCPEZ.NE.0 ) THEN
C     WE JUST READ A CONSTANT CHARGE INDEPENDENT VALUE OF ZETA
C     THESE SHOULD BE ATOMIC UNITS (1/BOHR) SINCE I DON'T USE A CONVERSION
         ZA = CPEZ(NAT(IA))
      ELSE
C     WE CALCULATE ZETA ON THE FLY BASED ON MULLIKEN CHARGE
C     THE ATOM POLARIZABILITY CALCULATED WITH SEMIEMPERICAL
         ATS = ATPS(NAT(IA))
C     THE ANION POLARIZABILITY CALCULATED WITH SEMIEMEPERICAL
         ANS = ANPS(NAT(IA))
C     THE TARGET POLARIZABILITY OF THE ATOM
         ATT = ATPT(NAT(IA))
C     THE TARGET POLARIZABILITY OF THE ANION
         ANT = ANPT(NAT(IA))
C     CONVERT FROM ANGSTROM**3 TO CUBIC BOHR
         DPAT = (ATT-ATS)*CUBANG
         DPAN = (ANT-ANS)*CUBANG
         IF ( DPAT .LT. 1.0D-6 .OR. DPAN .LT. 1.0D-6 ) THEN
C            IF ( NPTCPE.GT.0 ) THEN
C               WRITE(6,'(A,2E12.5)')"ZETA DPAT,DPAN=",DPAT,DPAN
C            END IF
            ZA = CPEZ(NAT(IA))
            GOTO 900
         END IF
C     WE ARE GOING TO DIVIDE BY SOME NUMBER
C     WE NEED TO CHECK TO FOR A DIVIDE BY ZERO
         IF ( DPAN .LT. 1.0D-6 ) THEN
C     OUR DENOMINATOR IS REALLY SMALL, BUT IF THE CHARGE IS 0,
C     THEN THAT TERM IS JUST UNITY ANYWAY
            IF ( QA**2 .LT. 1.0D-6 ) THEN
               TERMA = 1.0D0
            ELSE
C     TO AVOID A DIVIDE BY ZERO, WE SET THE TERM TO SOME BIG NUMBER
               TERMA = 1.0D6
            END IF
         ELSE
C     OUR DENOMINATOR IS REASONABLE OF REASONABLE SIZE
C     AND WE CAN ACTUALLY DIVIDE WITH IT
            TERMA = (ONE/DPAN)**(-QA/THREE)
         END IF
C     AGAIN, WE ARE GOING TO DIVIDE BY A NUMBER AND WE NEED
C     TO CHECK IF IT IS CLOSE TO ZERO
         IF ( DPAT .LT. 1.0D-6 ) THEN
C     OUR DENOMINATOR IS REALLY SMALL, BUT IF THE CHARGE IS
C     CLOSE TO -1, THEN THAT TERM IS JUST 1 ANYWAY
            IF ( (QA+ONE)**2 .LT. 1.0D-6 ) THEN
               TERMB = 1.0D0
            ELSE
C     TO AVOID A DIVIDE BY ZERO, WE SET THE TERM TO SOME BIG NUMBER
               TERMB = 1.0D6
            END IF
         ELSE
C     OUR DENOMINATOR IS REASONABLE
            TERMB = (ONE/DPAT)**((QA+ONE)/THREE)
         END IF
C     NOW WE COMBINE THE TERMS TO GET OUR ZETA
         ZA = PREFZ * TERMA * TERMB
      END IF
 900  CONTINUE
C      IF ( NPTCPE.GT.0 ) THEN
         WRITE(6,'(A,I3,2X,E12.4)')"ATOM,ZETA = ",NAT(IA),ZA
C      END IF
      END SUBROUTINE CPEZET
C***********************************************************************
C***********************************************************************
      SUBROUTINE NEWZET(IA,ZA,QA)
C     --------------------------------------------------------
C     CALCULATE THE ZETA EXPONENT OF THE CPE RESPONSE FUNCTION
C     --------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (LM1=600)
      PARAMETER (LMX=9*LM1)
      PARAMETER (LMZ=86)
      COMMON
     .     /ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     .     /CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
     .     /CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     .     /CPEOPT/ NFLCPE,NPTCPE,NCPEZ
     .     /CPEPAR/ CPEZ(LMZ),CPEED0(LMZ),CPEQ0(LMZ),ATPS(LMZ),
     .     ANPS(LMZ),ATPT(LMZ), ANPT(LMZ),CPESL(LMZ),CPEEM0(LMZ),
     .     CPEM0(LMZ),CPERS(LMZ),CPERR(LMZ),CMXCAT(LMZ),CMXAN(LMZ),
     .     CPEFT(LMZ),CPERC(LMZ),CPEQM(LMZ),CPEPS,CPEZA(LMZ),CPEZB(LMZ)
C     THESE ARE JUST SOME CONSTANTS
      PT333 = ONE/THREE
C      SQ2OPI = SQRT(TWO/PI)
C      PREFZ = (THREE/SQ2OPI)**PT333
C      CUBANG = ONE/A0
C      CUBANG = CUBANG**3
      CUBANG = (0.1889725982D1)**3
      IF ( NCPEZ.NE.0 ) THEN
C     WE JUST READ A CONSTANT CHARGE INDEPENDENT VALUE OF ZETA
C     THESE SHOULD BE ATOMIC UNITS (1/BOHR) SINCE I DON'T USE A CONVERSION
         ZA = CPEZ(NAT(IA))
      ELSE
         ZA = CPEZA(NAT(IA))*CUBANG
C     WE CALCULATE ZETA ON THE FLY BASED ON MULLIKEN CHARGE
C         WRITE(6,'(A,2F12.5)')"#  8",QA,EXP(CPEZB(NAT(IA))*QA)
         ZA = (THREE*SQRT(PI/TWO)/ZA)**PT333 * EXP(CPEZB(NAT(IA))*QA)
      END IF
C     THIS ENSURES THAT THE POLLARIZBILITY OF THE SITE IS INDEPENDENT OF
C     THE EFFECTIVE DIELECTRIC PREFACTOR, CPEPS.
      ZA = ZA*(CPEPS**PT333)
      END SUBROUTINE
c$$$C***********************************************************************
c$$$C***********************************************************************
c$$$      SUBROUTINE CPEDMP(DMP,RIJ,SL,Q0,RC)
c$$$C     --------------------------------------------------------
c$$$C     CALCULATE THE DAMPING OF THE RESPONSE BASIS WITH THE
c$$$C     STATIC CHARGE.
c$$$C     --------------------------------------------------------
c$$$      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c$$$      COMMON
c$$$
c$$$     .     /CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
c$$$     .     /CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
c$$$
c$$$C     THIS IS THR NATURAL LOG OF 0.05
c$$$C     0.05 IS ASSOCIATED WITH THE DAMPING OF 5 PERCENT
c$$$      LOG5 = -2.9957322736D0
c$$$
c$$$C     THE DAMPING TAKES THE FORM
c$$$C     1-EXP(A*RIJ**2)
c$$$C     WHERE A IS RELATED TO THE SCREENING LENGTH
c$$$C     THE SCREENING LENGTH IS CHOSEN SUCH THAT
c$$$C     AT THAT DISTANCE, THE DAMPING IS 95% TURNED *OFF*
c$$$C     WHICH MEANS, THE DAMPING HAS A VALUE OF 0.95
c$$$C     WE NEED TO CALCULATE THE PARAMETER A, BUT WE NEED
c$$$C     TO MAKE SURE WE DO NOT DIVIDE BY ZERO
c$$$
c$$$C     IF THE SCREENING LENGTH IS REALLY SMALL, THEN
c$$$C     THE THERE IS NO DAMPING AND THE RESULT IS UNITY
c$$$      SLB   = SL/A0
c$$$      RIJC   = RIJ-RC/A0
c$$$C      RIJC  = RIJ-RCB
c$$$c      RSLC2 = (SLB-RCB)**2
c$$$C      IF ( RIJC .LT. 0.00001 ) THEN
c$$$C         DMP = ZERO
c$$$C      ELSE
c$$$C         IF ( RSLC2 .LT. 0.0001 ) THEN
c$$$C            DMP = ONE
c$$$C         ELSE
c$$$C            DMP =
c$$$
c$$$      IF( SLB .LT. 1.0D-2 .OR. RC .LT. 0.00001D0 ) THEN
c$$$         DMP = ONE
c$$$      ELSE
c$$$!         DMP = ONE - EXP(RIJ**2 * LOG5/SLB**2)
c$$$         TMP = EXP(-SLB*RIJC)
c$$$         DMP = ONE - TMP/(ONE+TMP)
c$$$      END IF
c$$$
c$$$C     NOW WE SCALE THE DAMPING BY THE SCALING FACTOR, Q0
c$$$
c$$$      DMP = DMP*Q0
c$$$
c$$$      END SUBROUTINE CPEDMP
      SUBROUTINE CPEDMP(DMP,RIJ,NATIL,NATJL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (LM1=600)
      PARAMETER (LMX=9*LM1)
      PARAMETER (LMZ=86)
      COMMON
     .     /ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     .     /CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
     .     /CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
c     .     /CPEOPT/ NFLCPE,NPTCPE,NCPEZ
     .     /CPEPAR/ CPEZ(LMZ),CPEED0(LMZ),CPEQ0(LMZ),ATPS(LMZ),
     .     ANPS(LMZ),ATPT(LMZ), ANPT(LMZ),CPESL(LMZ),CPEEM0(LMZ),
     .     CPEM0(LMZ),CPERS(LMZ),CPERR(LMZ),CMXCAT(LMZ),CMXAN(LMZ),
     .     CPEFT(LMZ),CPERC(LMZ),CPEQM(LMZ),CPEPS,CPEZA(LMZ),CPEZB(LMZ)
      SIX = FOUR+TWO
      TEN = SIX+FOUR
      FIFTN = TEN+SIX-ONE
      RLOW  = (CPERC(NATIL)-CPESL(NATIL)/2)/A0
      RLEN  = CPESL(NATIL)/A0
      RHIGH = RLOW+RLEN
      IF ( RIJ .LE. RLOW ) THEN
         DMP = ZERO
         RETURN
      ELSE IF ( RIJ .GE. RHIGH ) THEN
         DMP = ONE
      ELSE
         X = (RIJ-RLOW)/RLEN
         X2=X*X
         X3=X2*X
         DMP = X3*(TEN-FIFTN*X+SIX*X2)
      END IF
C      WRITE(6,'(A,5F12.5)')"DMP1",DMP,RLOW,RHIGH,RIJ,X
      RLOW  = (CPERC(NATJL)-CPESL(NATJL)/2)/A0
      RLEN = CPESL(NATJL)/A0
      RHIGH = RLOW+RLEN
      IF ( RIJ .LE. RLOW ) THEN
         DMP = ZERO
         RETURN
      ELSE IF ( RIJ .GE. RHIGH ) THEN
C         RETURN
      ELSE
         X = (RIJ-RLOW)/RLEN
         X2=X*X
         X3=X2*X
         DMP = DMP * X3*(TEN-FIFTN*X+SIX*X2)
      END IF
C      WRITE(6,'(A,5F12.5)')"DMP2",DMP,RLOW,RHIGH,RIJ,X
C      DMP = DMP*CPEPS
      END SUBROUTINE
      SUBROUTINE CLRP(ETA,NU,N,NCPE)
      INTEGER N,NCPE
      DOUBLE PRECISION ETA(NCPE,NCPE),NU(NCPE)
      INTEGER I,J,IBLK,JBLK,K,L
      DO 10 I=1,N
         IBLK = (I-1)*4+1
         DO 20 J=1,N
            JBLK = (J-1)*4+1
            DO 30 K=1,3
               DO 40 L=1,3
                  IF (K.EQ.L) CYCLE
                  ETA(IBLK+K,JBLK+L) = 0.0D0
 40            CONTINUE
 30         CONTINUE
 20      CONTINUE
         DO 50 K=1,3
            NU(IBLK+K) = 0.0D0
            ETA(IBLK+K,JBLK) = 0.0D0
            ETA(IBLK,JBLK+K) = 0.0D0
 50      CONTINUE
 10   CONTINUE
      CALL PTM(NCPE,ETA,NCPE,4,4,"CLRP:ETA",6,.FALSE.)
      CALL PTV(NU,NCPE,4,"CLRP:NU ",6,.FALSE.)
      END SUBROUTINE
c$$$
c$$$
c$$$
c$$$      ELSE IF ( DOSS ) THEN
c$$$C*****************************************************
c$$$C     SUBSYSTEM NORMALIZATION
c$$$C*****************************************************
c$$$
c$$$C     CHOOSE THE TYPE OF SUBSYSTEMS
c$$$         IF ( NPTCPE.GT.0 ) THEN
c$$$           WRITE(6,'(A,A)')"SOLVING FOR RESPONSE UNDER ",
c$$$     A           "SUBSYSTEM CONSTRAINTS"
c$$$        END IF
c$$$        IF ( NFLCPE.LE.3 ) THEN
c$$$           CALL CPEMNM(PA,PB,LM4,LM2,EIGVEC,EIGVAL,PMOL)
c$$$c           CALL CPEMNM(PA,PB,LM4,LM2,PMOL)
c$$$         ELSE
c$$$            WRITE(6,*)"**********************"
c$$$            WRITE(6,*)"FATAL ERROR"
c$$$            WRITE(6,*)"FLCPE = ",NFLCPE
c$$$            WRITE(6,*)"NO SUBSYSTEM NORMALIZATION IS DEFINED"
c$$$            WRITE(6,*)"FOR THIS INPUT OPTION"
c$$$            WRITE(6,*)"THIS IS A USER ERROR, NOT A PROGRAMMING ERROR"
c$$$            WRITE(6,*)"ABNORMAL TERMINATION"
c$$$            STOP
c$$$         END IF
c$$$c------------------- PRINT SECTION
c$$$         IF ( PTFILE ) THEN
c$$$            DO 301 I=1,NUMAT
c$$$               DO 302 J=1,NUMAT
c$$$                  WRITE(58,'(E28.20)')PMOL(I,J)
c$$$ 302           CONTINUE
c$$$ 301        CONTINUE
c$$$         END IF
c$$$         IF ( NPTCPE.GT.0 ) THEN
c$$$            WRITE(6,'(/A)')"CHARGE TRANSFER PROBABILITY MATRIX"
c$$$            WRITE(6,'(3X,30(I4,3X))')(NAT(I),I=1,NUMAT)
c$$$            DO 299 I=1,NUMAT
c$$$               WRITE(6,'(I3,30F7.3)')NAT(I),(PMOL(I,J),J=1,NUMAT)
c$$$ 299        CONTINUE
c$$$            WRITE(6,'(/A)')"EIGENVALUE | EIGENVECTOR"
c$$$            WRITE(6,'(9X,30(I4,3X))')(NAT(I),I=1,NUMAT)
c$$$            DO 300 I=1,NUMAT
c$$$               WRITE(6,'(F6.3,A,30F7.3)')EIGVAL(I)," | ",
c$$$     A              (EIGVEC(J,I),J=1,NUMAT)
c$$$ 300        CONTINUE
c$$$         END IF
c$$$c--------------------
c$$$C     CONVERT THE ORTHOGONAL SUBSYSTEMS TO THE CONSTRAINT MATRIX
c$$$C     (PROGRAMMING NOTE: PMOL IS NO LONGER USED AND CAN BE REMOVED)
c$$$         CALL CNVRTP(PMOL,EIGVEC,EIGVAL,NUMAT)
c$$$c--------------------PRINT SECTION
c$$$         IF ( NPTCPE.GT.0 ) THEN
c$$$            WRITE(6,'(/A)')"CONSTRAINT MATRIX"
c$$$            WRITE(6,'(30(I4,3X))')(NAT(I),I=1,NUMAT)
c$$$            DO 297 I=1,NUMAT
c$$$               IF ( ABS(EIGVAL(I)).GT.0.0001D0 ) THEN
c$$$                  WRITE(6,'(30F7.3)')(EIGVEC(J,I),J=1,NUMAT)
c$$$               END IF
c$$$ 297        CONTINUE
c$$$         END IF
c$$$c--------------------
c$$$C     DETERMINE THE NUMBER OF CONSTRAINTS
c$$$C     ...WE ONLY USE THOSE ORTHOGONAL SUBSYSTEMS
c$$$C     WHICH HAVE A NONZERO EIGENVALUE
c$$$         DO 296 I=1,NUMAT
c$$$            IF ( EIGVAL(I) .GT. 0.0001D0 ) EXIT
c$$$ 296     CONTINUE
c$$$         LCON = NUMAT-I+1
c$$$         ICON = I-1
c$$$c         LCON = NUMAT
c$$$c         ICON = 0
c$$$         IF ( LCON .EQ. 0 ) THEN
c$$$            DOSS = .FALSE.
c$$$            GOTO 1298
c$$$         END IF
c$$$C     CALCULATE ETA INVERSE
c$$$         CALL SYMINV(CPEETA,NCPE,CPEINV)
c$$$C     SOLVE FOR THE RESPONSE COEFFICIENTS
c$$$         CALL SSCPE(CPERES,CPEINV,CPENU,EIGVEC,EIGVAL,ICON,LCON)
c$$$C     SKIP AHEAD TO THE ENERGY AND DIPOLE CALCULATION
c$$$C ********************************
c$$$      END IF
      SUBROUTINE HYBDIP(HDIP,PA,PB,LM4,PRT)
C     COMPUTES THE DIPOLE MOMENT OF AN ATOM AND STORES IT IN HDIP
C     *
C     NOTATION. I=INPUT, O=OUTPUT.
C     PA(LM4)   RHF OR UHF-ALPHA DENSITY MATRIX (I).
C     PB(LM4)   UHF-BETA DENSITY MATRIX (I).
C     PRT       LOGICAL PRINTING FLAG (I).
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (LM1=600)
      PARAMETER (LMX=9*LM1)
      PARAMETER (LMZ=86)
      CHARACTER*2 ELEMNT
      LOGICAL PRT
      COMMON
C     ./AMASS / AMS(LM1)
     ./ATOMC / COORD(3,LM1)
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
C     ./CHARGE/ QI(LM1)
C     ./CHARGP/ POP(3,LM1)
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./DIPOL / DM,DMX(3)
     ./DIPOL1/ HYF(LMZ),HYFPD(LMZ)
     ./ELEMTS/ ELEMNT(107)
C.##IF CHARMM
C.##ELSE
     ./INDEX / INDX(LMX)
C.##ENDIF
     ./INOPT2/ IN2(200)
     ./PARDER/ CORE(LMZ),EHEAT(LMZ),EISOL(LMZ)
      DIMENSION PA(LM4),PB(LM4)
      DIMENSION CG(3),DM2(4),PAO(9)
      DIMENSION HDIP(3,NUMAT)
C *** INPUT OPTIONS.
C      KHARGE = IN2(65)
C *** INITIALIZATION.
C      DM     = ZERO
C      DO 10 I=1,3
C      CG(I)  = ZERO
C   10 DMX(I) = ZERO
      DO 20 I=1,4
C     DM1(I) = ZERO
         DM2(I) = ZERO
 20   CONTINUE
C     *** LOOP OVER ALL ATOMS.
      DO 70 I=1,NUMAT
         NI     = NAT(I)
         IA     = NFIRST(I)
         IB     = NLAST(I)
         IORBS  = IB-IA+1
C     *** CHARGE CONTRIBUTIONS TO DIPOLE MOMENT.
C     CONVERSION FACTOR : 1 e  =  4.803242D-10 esu
C     TAKEN FROM J.PHYS.CHEM.REF.DATA 2, 663 (1973).
C     DO 50 J=1,3
C     50 DM1(J) = DM1(J)+4.80324D0*QI(I)*(COORD(J,I)-CG(J))
C     *** SP HYBRIDIZATION CONTRIBUTIONS TO DIPOLE MOMENT.
         IF(IORBS.EQ.1) GO TO 70
         DO 60 J=1,3
            L      = INDX(IA+J)+IA
            DM2(J) = -HYF(NI)*(PA(L)+PB(L))
 60      CONTINUE
C     *** PD HYBRIDIZATION CONTRIBUTIONS TO DIPOLE MOMENT.
         HDIP(1,I) = DM2(1) / (A0*4.80324D0)
         HDIP(2,I) = DM2(2) / (A0*4.80324D0)
         HDIP(3,I) = DM2(3) / (A0*4.80324D0)
         IF(IORBS.EQ.4) GO TO 70
         XT     = ONE/SQRT(THREE)
C     REQUIRED SUMS OF DENSITY MATRIX ELEMENTS.
C     DX     = P(XZ,Z) + P(X2Y2,X) + P(XY,Y) - ONE/SQRT(THREE)*P(Z2,X)
C     DY     = P(YZ,Z) - P(X2Y2,Y) + P(XY,X) - ONE/SQRT(THREE)*P(Z2,Y)
C     DZ     = P(XZ,X) +             P(YZ,Y) + TWO/SQRT(THREE)*P(Z2,Z)
         LL     = INDX(IA+5) + IA+3
         DX     = PA(LL)+PB(LL)
         LL     = INDX(IA+4) + IA+1
         DX     = DX+PA(LL)+PB(LL)
         LL     = INDX(IA+8) + IA+2
         DX     = DX+PA(LL)+PB(LL)
         LL     = INDX(IA+6) + IA+1
         DX     = DX-XT*(PA(LL)+PB(LL))
         LL     = INDX(IA+7) + IA+3
         DY     = PA(LL)+PB(LL)
         LL     = INDX(IA+4) + IA+2
         DY     = DY-(PA(LL)+PB(LL))
         LL     = INDX(IA+8) + IA+1
         DY     = DY+PA(LL)+PB(LL)
         LL     = INDX(IA+6) + IA+2
         DY     = DY-XT*(PA(LL)+PB(LL))
         LL     = INDX(IA+5) + IA+1
         DZ     = PA(LL)+PB(LL)
         LL     = INDX(IA+7) + IA+2
         DZ     = DZ+PA(LL)+PB(LL)
         LL     = INDX(IA+6) + IA+3
         DZ     = DZ+TWO*XT*(PA(LL)+PB(LL))
         DM2(1) = DM2(1) - DX*HYFPD(NI)
         DM2(2) = DM2(2) - DY*HYFPD(NI)
         DM2(3) = DM2(3) - DZ*HYFPD(NI)
         HDIP(1,I) = HDIP(1,I) + DM2(1) / (A0*4.80324D0)
         HDIP(2,I) = HDIP(2,I) + DM2(2) / (A0*4.80324D0)
         HDIP(3,I) = HDIP(3,I) + DM2(3) / (A0*4.80324D0)
 70   CONTINUE
      END SUBROUTINE
      SUBROUTINE BOSCL(PA,PB,LM4,LM2,NUMAT,SCLBND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (LM1=600)
      PARAMETER (LMX=9*LM1)
      PARAMETER (LMZ=86)
      INTEGER LM4,LM2
      DIMENSION PA(LM4),PB(LM4)
C      COMMON
C     .     /ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
C     .     /CPEOPT/ NFLCPE,NPTCPE,NCPEZ
C      DOUBLE PRECISION BNDMAT,MAXPIJ,R,RLO,RHI,SW
      DIMENSION SCLBND(NUMAT,NUMAT)
      DIMENSION BNDMAT(NUMAT,NUMAT)
      LOGICAL DPRINT
      INTEGER I,J,K,NUMAT
      DPRINT = .FALSE.
      CALL BNDPOP(PA,PB,LM4,LM2,BNDMAT,DPRINT)
C     CONVERT BOND ORDER MATRIX TO A PROBABILITY MATRIX
      RLO = 0.00D0
      RHI = 0.30D0
      DO 10 I=1,NUMAT-1
         DO 11 J=I+1,NUMAT
            R = BNDMAT(I,J)
            CALL GENSW( R, RLO ,RHI ,SW )
            BNDMAT(I,J) = SW
            BNDMAT(J,I) = SW
 11      CONTINUE
         BNDMAT(I,I) = 1.0D0
 10   CONTINUE
      BNDMAT(NUMAT,NUMAT) = 1.0D0
c$$$C     FIND THE MAXIMUM PROBABILITY OF BEING IN THE SAME RESIDUE
      DO 20 I=1,NUMAT-1
         DO 21 J=I+1,NUMAT
C------------------------------------------------------------------
C     1-2 interactions
            R = BNDMAT(I,J)
            IF ( R .GE. 0.99999D0 ) GOTO 999
C     1-3 interactions
            DO 22 K=1,NUMAT
               IF ( K .NE. I .AND. K .NE. J ) THEN
                  Z = BNDMAT(I,K)*BNDMAT(K,J)
                  IF ( Z .GT. R ) THEN
                     R = Z
                     IF ( R .GT. 0.99999D0 ) GOTO 999
                  END IF
               END IF
C     1-4 interactions
               DO 23 L=1,NUMAT
                  IF ( K .NE. I .AND. K .NE. J .AND.
     1                 L.NE.I.AND.L.NE.J.AND.L.NE.K ) THEN
                     Z = BNDMAT(I,K)*BNDMAT(K,L)*BNDMAT(L,J)
                     IF ( Z .GT. R ) THEN
                        R = Z
                        IF ( R .GT. 0.99999D0 ) GOTO 999
                     END IF
                  END IF
 23            CONTINUE
 22         CONTINUE
C-----------------------------------------------------------
 999        CONTINUE
            SCLBND(I,J) = R
            SCLBND(J,I) = R
 21      CONTINUE
         SCLBND(I,I) = 1.0D0
 20   CONTINUE
      SCLBND(NUMAT,NUMAT) = 1.0D0
      DO 200 I=1,NUMAT
         DO 201 J=1,NUMAT
            SCLBND(I,J) = 1.0D0 - SCLBND(I,J)
 201     CONTINUE
 200  CONTINUE
c$$$
c$$$      IF ( NPTCPE .NE. 0 ) THEN
         DPRINT = .TRUE.
c$$$      ELSE
c$$$         DPRINT = .FALSE.
c$$$      END IF
c      IF ( DPRINT ) THEN
c         WRITE(6,*)"------- BOND SCALE MATRIX ----------"
c         DO 50 I=1,NUMAT
c            WRITE(6,'(100F9.2)')(SCLBND(I,J),J=1,NUMAT)
c 50      CONTINUE
c      END IF
C     TURN OFF BOND DAMPING
C         DO 300 I=1,NUMAT
C            DO 301 J=1,NUMAT
C               SCLBND(I,J) = 1.0D0
C 301        CONTINUE
C 300     CONTINUE
      END SUBROUTINE
      SUBROUTINE GENSW( R,LO,HI,SW )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION R,LO,HI,WID,X,X2,X3,SW
      IF ( R .GE. HI ) THEN
         SW = 1.0D0
      ELSE IF ( R .LE. LO ) THEN
         SW = 0.0D0
      ELSE
         WID = HI-LO
         IF ( WID .LE. 0.00001D0 ) THEN
            WRITE(6,*)"WID TOO LOW IN GENSW",LO,HI,WID
            STOP
         END IF
C         WRITE(6,*)R,LO,HI,WID
         X = (R-LO)/WID
         X2=X*X
         X3=X2*X
         SW = X3*(10.0D0-15.0D0*X+6.0D0*X2)
      END IF
      END SUBROUTINE
