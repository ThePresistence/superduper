      SUBROUTINE NACANA (C,E,F,H,W,ONEDEN,ARRAY,CG,CC,DELTAE,
     1                   LM2,LM3,LM4,LM5,LM6)
C *** (PARTIALLY) ANALYTICAL EVALUATION OF THE SECOND TERM (MO COEFFICIENT
C     AND AO DERIVATIVES) FOR NON-ADIABATIC COUPLING ELEMENTS.
C     ARGUMENTS:
C       C:      MOLECULAR ORBITAL COEFFICIENTS (I).
C       E:      MOLECULAR ORBITAL ENERGIES, EV (I).
C       F:      FOCK MATRIX, EV (S).
C       H:      CORE HAMILTONIAN, EV (S).
C       ONEDEN: ONE-PARTICLE TRANSITION DENSITY MATRIX (I).
C       ARRAY:  UNNAMED COMMON (I,O,S).
C       CG:     FIRST TERM OF THE COUPLING WITHOUT ENERGY DENOMINATOR,
C               KCAL/(MOL*ANG) (I).
C       CC:     FULL NON-ADIABATIC COUPLING VECTOR, 1/ANG (O).
C       DELTAE: ENERGY DIFFERENCE OF STATES INVOLVED, KCAL/MOL (I).
C
      USE LIMIT, ONLY: LM1, LMX, LM1M, LMACT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./ATOMC / COORD(3,LM1)
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./CIMOS / IMOCI(LMX)
     ./GUGA1 / NCIO,NCIGAM
     ./HALFE / IODD,JODD
     ./INOPT2/ IN2(300)
     ./NBFILE/ NBF(20)
     ./ORBITS/ NUMB,NORBS,NMOS,NALPHA,NBETA
     ./PSDOPT/ DSTORE,DSTEP,RPSOPT(3:10),IPSOPT(28)
     ./QMMM1 / COORDM(3,LM1M),CHARGM(LM1M)
     ./QMMM6 / ISELCT(LM1+LM1M)
      DIMENSION C(LM2,LM3)
      DIMENSION E(LM3)
      DIMENSION F(LM4)
      DIMENSION H(LM4)
      DIMENSION W(LM6,LM6)
      DIMENSION ARRAY(LM5)
      DIMENSION ONEDEN(LMACT,NCIO)
      DIMENSION CG(3,LM1+LM1M)
      DIMENSION CC(3,LM1+LM1M)
      DIMENSION OCC(NORBS)
      ALLOCATABLE CORSAV(:,:), WMO(:,:), IMO(:), JMO(:)
      ALLOCATABLE GVEC(:), ATRANS(:,:), RHS(:), ZVEC(:), IPIV(:)
      ALLOCATABLE FSAV(:), FDERIV(:,:), FDERMO(:,:)
      ALLOCATABLE SMIPT5(:,:), SPLPT5(:,:)
      ALLOCATABLE SMI1(:,:), SMI2(:,:)
      ALLOCATABLE DSMI(:,:), SIGMA(:,:), TEMP(:,:)
      ALLOCATABLE TERM1(:,:), SUM12A(:,:), TERM2B(:,:)
      LOGICAL DEBUG
C *** LMPAIR IS THE TOTAL NUMBER OF MO PAIRS.
C     NPAIR1 IS THE NUMBER OF MO PAIRS WITH DIFFERENT OCCUPATION.
C     NPAIR2 IS THE TOTAL NUMBER OF NON-DEGENERATE MO PAIRS.
      INTEGER LMPAIR,NPAIR1,NPAIR2,I,J,INFO
C *** INITIALIZATION.
C     FILE NUMBERS.
      NB6     = NBF(6)
C     INPUT OPTIONS.
      MPRINT = IN2( 41)
      NUMATM = IN2(120)
      ICI1   = IN2(131)
      ICI2   = IN2(132)
      NACTIV = ICI1 + ICI2
      IF(DSTEP.EQ.0.D0)  DSTEP = 2.D-4
      DEBUG  = .FALSE.
*     DEBUG  = .TRUE.
      NUMALL = NUMAT + NUMATM
      LMPAIR = LM4   - NORBS
      FACTOR = 0.5D0 / DSTEP
      ALLOCATE(CORSAV(3,NUMAT))
      ALLOCATE(WMO(LM4,LM4))
      ALLOCATE(IMO(LMPAIR))
      ALLOCATE(JMO(LMPAIR))
C     GENERATE LIST OF MO PAIRS AND DETERMINE NPAIR1 ETC.
      CALL MOPAIR
      ALLOCATE(GVEC(LM4))
      ALLOCATE(ATRANS(NPAIR1,NPAIR1))
      ALLOCATE(RHS(NPAIR1))
      ALLOCATE(ZVEC(LM4))
      ALLOCATE(FSAV(LM4))
      ALLOCATE(FDERIV(NORBS,NORBS))
      ALLOCATE(FDERMO(NORBS,NORBS))
      ALLOCATE(SMIPT5(NORBS,NORBS))
      ALLOCATE(SPLPT5(NORBS,NORBS))
      ALLOCATE(SMI1(NORBS,NORBS))
      ALLOCATE(SMI2(NORBS,NORBS))
      ALLOCATE(DSMI(NORBS,NORBS))
      ALLOCATE(SIGMA(NORBS,NORBS))
      ALLOCATE(TEMP(NORBS,NORBS))
      ALLOCATE(TERM1(3,NUMALL))
      ALLOCATE(SUM12A(3,NUMALL))
      ALLOCATE(TERM2B(3,NUMALL))
      CORSAV(1:3,1:NUMAT) = COORD(1:3,1:NUMAT)
      OCC(1     :NUMB)    = 2.D0
      OCC(NUMB+1:NORBS)   = 0.D0
      IF(IODD.GT.0) OCC(IODD) = 1.D0
      IF(JODD.GT.0) OCC(JODD) = 1.D0
C     DEBUG PRINT.
      IF(DEBUG) THEN
         WRITE(NB6,'(///1X,
     1         "CARTESIAN COORDINATES (ORIGINAL VALUES)."/)')
         CALL PRTCOR(CORSAV(:,1:NUMAT))
         WRITE(NB6,'(///1X,"ONE-PARTICLE TRANSITION DENSITY MATRIX.")')
         CALL PRTMAT(ONEDEN(1:NACTIV,1:NACTIV))
      ENDIF
C *** TRANSFORM TWO ELECTRON INTEGRALS INTO MO BASIS.
      IF(DEBUG) THEN
         WRITE(NB6,'(//1X,"AO BASIS ELECTRON REPULSION INTEGRALS.")')
         CALL PRTMAT(W)
      END IF
      CALL ERI2MO(C,WMO)
      IF(DEBUG) THEN
         WRITE(NB6,'(//1X,"MO BASIS ELECTRON REPULSION INTEGRALS.")')
         CALL PRTMAT(WMO)
      END IF
C *** POWERS OF OVERLAP MATRIX.
C     FORM S**(-1/2) MATRIX.
      CALL SMICAL(SMIPT5,CORSAV,NORBS,DEBUG)
C     FORM S**(+1/2) MATRIX.
      CALL SPLCAL(SPLPT5,CORSAV,NORBS,DEBUG)
C
C *** DERIVATIVE-INDEPENDENT PART OF THE Z VECTOR PROCEDURE.
C
C     MAKE SURE Z VECTOR EQUATIONS CORRESPONDING TO
C     DEGENERATE MO PAIRS ARE SATISFIED.
      CALL ASSERT
C     FORM COEFFICIENT MATRIX OF Z VECTOR LSE.
      CALL MAKATR
C     FORM RIGHT-HAND SIDE OF Z VECTOR LSE.
      CALL MAKRHS
C     SOLVE Z VECTOR LSE.
      ALLOCATE(IPIV(NPAIR1))
      CALL DGESV(NPAIR1,1,ATRANS(1,1),NPAIR1,IPIV(1),RHS(1),NPAIR1,INFO)
      IF(INFO.NE.0) THEN
         WRITE(NB6,'(1X,"DGESV FAILED WITH INFO =",I6)') INFO
         STOP 'NACANA'
      END IF
      DEALLOCATE(IPIV)
      DO IPAIR=1,NPAIR1
         I        = IMO(IPAIR)
         J        = JMO(IPAIR)
         IJ       = KITE(I,J)
         ZVEC(IJ) = RHS(IPAIR)
      END DO
C
C     SHIFT EVERY NUCLEAR COORDINATE IN BOTH DIRECTIONS AND EACH TIME
C     DETERMINE THE TRANSFORMATION MATRIX CONVERTING THE ORIGINAL MOS
C     TO THE MOS OF THE DISTORTED STRUCTURE. THESE TRANSFORMATION
C     MATRICES SHOULD BE VERY CLOSE TO THE UNIT MATRIX. FORM THE
C     DERIVATIVE OF THE TRANSFORMATION MATRIX DIVIDING THE DIFFERENCE
C     BY TWO TIMES THE SHIFT.
C
C *** CONTRIBUTIONS FROM THE VARIATION OF THE POSITIONS OF THE QM ATOMS.
C
      DO IA=1,NUMAT
         DO IC=1,3
            IF(DEBUG) THEN
               WRITE(NB6,'(///1X,"ATOM",I4,", DIRECTION",I2,".")') IA,IC
            ENDIF
            Q0           = COORD(IC,IA)
            COORD(IC,IA) = Q0 + DSTEP
            IF(DEBUG) THEN
               WRITE(NB6,'(///1X,
     1               "CARTESIAN COORDINATES (POSITIVE SHIFT)."/)')
               CALL PRTCOR(COORD(:,1:NUMAT))
            ENDIF
            CALL COMPHF(ARRAY,LM5,IA,ESCAL)
C           F ONLY CONTAINS THE CONTRIBUTIONS FROM THE
C           TWO-ELECTRON INTEGRALS TO THE FOCK MATRIX.
            FSAV(1:LM4)  = H(1:LM4) + F(1:LM4)
            CALL SMICAL(SMI1,COORD,NORBS,DEBUG)
            COORD(IC,IA) = Q0 - DSTEP
            IF(DEBUG) THEN
               WRITE(NB6,'(///1X,
     1               "CARTESIAN COORDINATES (NEGATIVE SHIFT)."/)')
               CALL PRTCOR(COORD(:,1:NUMAT))
            ENDIF
            CALL COMPHF(ARRAY,LM5,IA,ESCAL)
            CALL SMICAL(SMI2,COORD,NORBS,DEBUG)
            COORD(IC,IA) = Q0
            IJ = 0
            DO I=1,NORBS
               DO J=1,I
                  IJ          = IJ + 1
                  FDERIV(I,J) = FACTOR * (FSAV(IJ) - H(IJ) - F(IJ))
                  FDERIV(J,I) = FDERIV(I,J)
               END DO
            END DO
            IF(DEBUG) THEN
               WRITE(NB6,'(///1X,"AO BASIS FOCK MATRIX DERIVATIVE.")')
               CALL PRTMAT(FDERIV)
            ENDIF
            CALL SIMTRA(FDERIV,C,FDERMO)
            IF(DEBUG) THEN
               WRITE(NB6,'(///1X,"MO BASIS FOCK MATRIX DERIVATIVE.")')
               CALL PRTMAT(FDERMO)
            ENDIF
            DSMI(:,:) = FACTOR * (SMI1(:,:) - SMI2(:,:))
            IF(DEBUG) THEN
               WRITE(NB6,'(///1X,"S**(-1/2) MATRIX DERIVATIVE.")')
               CALL PRTMAT(DSMI)
            ENDIF
C           DIVIDE BY ENERGY DIFFERENCE.
            TERM1(IC,IA) = CG(IC,IA) / DELTAE
C           CONTRIBUTION FROM THE VARIATION OF THE MOS.
            TMP = 0.D0
            IJ  = 0
            DO I=1,NORBS
               DO J=1,I
                  IJ  = IJ  + 1
                  TMP = TMP + ZVEC(IJ) * FDERMO(I,J)
               END DO
            END DO
            SUM12A(IC,IA) = TERM1(IC,IA) + TMP
C           FORM SIGMA MATRIX IN NON-ORTHOGONAL AO BASIS.
            CALL SIGNUM(CORSAV,COORD,SIGMA,IA,IC)
            IF(DEBUG) THEN
               WRITE(NB6,'(///1X,"AO BASIS SIGMA MATRIX.")')
               CALL PRTMAT(SIGMA)
            ENDIF
C           TRANSFORM SIGMA MATRIX INTO ORTHOGONAL AO BASIS.
C           TEMP = SMIPT5 * SIGMA.
            CALL DGEMM('N','N',NORBS,NORBS,NORBS,1.D0,SMIPT5(1,1),NORBS,
     1                 SIGMA(1,1),NORBS,0.D0,TEMP(1,1),NORBS)
C           SIGMA = TEMP * SMIPT5.
            CALL DGEMM('N','N',NORBS,NORBS,NORBS,1.D0,TEMP(1,1),NORBS,
     1                 SMIPT5(1,1),NORBS,0.D0,SIGMA(1,1),NORBS)
            IF(DEBUG) THEN
               WRITE(NB6,'(///1X,"ORTHOGONAL AO BASIS SIGMA MATRIX.")')
               CALL PRTMAT(SIGMA)
            ENDIF
C           FORM SIGMA = SIGMA + SPLPT5 * DSMI.
            CALL DGEMM('N','N',NORBS,NORBS,NORBS,1.D0,SPLPT5(1,1),NORBS,
     1                 DSMI(1,1),NORBS,1.D0,SIGMA(1,1),NORBS)
            IF(DEBUG) THEN
               WRITE(NB6,'(///1X,
     1               "ANTISYMMETRIC AO BASIS SIGMA MATRIX.")')
               CALL PRTMAT(SIGMA)
            ENDIF
C           TRANSFORM SIGMA MATRIX INTO MO BASIS.
C           TEMP = C(TRANS) * SIGMA.
            CALL DGEMM('T','N',NORBS,NORBS,NORBS,1.D0,C(1,1),NORBS,
     1                 SIGMA(1,1),NORBS,0.D0,TEMP(1,1),NORBS)
C           SIGMA = TEMP * C.
            CALL DGEMM('N','N',NORBS,NORBS,NORBS,1.D0,TEMP(1,1),NORBS,
     1                 C(1,1),NORBS,0.D0,SIGMA(1,1),NORBS)
            IF(DEBUG) THEN
               WRITE(NB6,'(///1X,"MO BASIS SIGMA MATRIX.")')
               CALL PRTMAT(SIGMA)
            ENDIF
C           SIGMA CONTRIBUTIONS TO THE NON-ADIABATIC COUPLING.
            TMP = 0.D0
*           NO DIAGONAL CONTRIBUTION FOR ANTISYMMETRIC MATRIX.
*           DO I=1,NUMB
*              TMP = TMP + 2.D0 * SIGMA(I,I)
*           END DO
*           IF(IODD(1).GT.0)  TMP = TMP - SIGMA(IODD(1), IODD(1))
*           IF(IODD(2).GT.0)  TMP = TMP - SIGMA(IODD(2), IODD(2))
            DO I=1,NACTIV
               DO J=1,NACTIV
                  TMP = TMP + ONEDEN(I,J) * SIGMA(IMOCI(I), IMOCI(J))
               END DO
            END DO
            TERM2B(IC,IA) = TMP
            CC    (IC,IA) = SUM12A(IC,IA) + TMP
         END DO
      END DO
C
C *** CONTRIBUTIONS FROM THE VARIATION OF THE POSITIONS OF THE
C     MM POINT CHARGES.
C
C     CONTRIBUTIONS FROM VARIATION OF THE MO COEFFICIENTS ONLY.
C     NO VARIATION OF THE ATOMIC ORBITALS, I.E. CONSTANT OVERLAP
C     MATRIX AND ZERO AO DERIVATIVES.
      DO IAMM=1,NUMATM
         IA = IAMM + NUMAT
         IF(ISELCT(IA).LE.0) THEN
            DO IC=1,3
               IF(DEBUG) THEN
                  WRITE(NB6,'(///1X,"POINT CHARGE",I4,
     1                              ", DIRECTION",I2,".")') IAMM,IC
               ENDIF
               Q0              = COORDM(IC,IAMM)
               COORDM(IC,IAMM) = Q0 + DSTEP
               IF(DEBUG) THEN
                  WRITE(NB6,'(///1X,
     1                  "CARTESIAN COORDINATES (POSITIVE SHIFT)."/)')
                  CALL PRTCOR(COORDM(:,1:NUMATM))
               ENDIF
               CALL COMPHF(ARRAY,LM5,IA,ESCAL)
C              F ONLY CONTAINS THE CONTRIBUTIONS FROM THE
C              TWO-ELECTRON INTEGRALS TO THE FOCK MATRIX.
               FSAV(1:LM4)     = H(1:LM4) + F(1:LM4)
               COORDM(IC,IAMM) = Q0 - DSTEP
               IF(DEBUG) THEN
                  WRITE(NB6,'(///1X,
     1                  "CARTESIAN COORDINATES (NEGATIVE SHIFT)."/)')
                  CALL PRTCOR(COORD(:,1:NUMAT))
               ENDIF
               CALL COMPHF(ARRAY,LM5,IA,ESCAL)
               COORDM(IC,IAMM) = Q0
               IJ = 0
               DO I=1,NORBS
                  DO J=1,I
                     IJ          = IJ + 1
                     FDERIV(I,J) = FACTOR * (FSAV(IJ) - H(IJ) - F(IJ))
                     FDERIV(J,I) = FDERIV(I,J)
                  END DO
               END DO
               IF(DEBUG) THEN
                  WRITE(NB6,
     1                  '(///1X,"AO BASIS FOCK MATRIX DERIVATIVE.")')
                  CALL PRTMAT(FDERIV)
               ENDIF
               CALL SIMTRA(FDERIV,C,FDERMO)
               IF(DEBUG) THEN
                  WRITE(NB6,
     1                  '(///1X,"MO BASIS FOCK MATRIX DERIVATIVE.")')
                  CALL PRTMAT(FDERMO)
               ENDIF
C              DIVIDE BY ENERGY DIFFERENCE.
               TERM1(IC,IA) = CG(IC,IA) / DELTAE
C              NO CONTRIBUTION FROM THE SCF REFERENCE DUE TO ANTISYMMETRY
C              OF THE U MATRIX DERIVATIVE (ZERO DIAGONAL).
               TMP = 0.D0
               IJ  = 0
               DO I=1,NORBS
                  DO J=1,I
                     IJ  = IJ  + 1
                     TMP = TMP + ZVEC(IJ) * FDERMO(I,J)
                  END DO
               END DO
               SUM12A(IC,IA) = TERM1 (IC,IA) + TMP
               TERM2B(IC,IA) = 0.D0
               CC    (IC,IA) = SUM12A(IC,IA)
            END DO
         ELSE
            TERM1 (:,IA) = 0.D0
            SUM12A(:,IA) = 0.D0
            TERM2B(:,IA) = 0.D0
            CC    (:,IA) = 0.D0
         ENDIF
      END DO
      IF(MPRINT.GT.0) THEN
         WRITE(NB6,'(///1X,"NON-ADIABATIC COUPLING MATRIX ELEMENTS.")')
         WRITE(NB6,'(/1X,"FIRST TERM."/)')
         CALL PRTCOR(TERM1 (:,1:NUMALL))
         WRITE(NB6,'(/1X,"FIRST TERM PLUS SECOND TERM, FIRST PART."/)')
         CALL PRTCOR(SUM12A(:,1:NUMALL))
         WRITE(NB6,'(/1X,"SECOND TERM, SECOND PART."/)')
         CALL PRTCOR(TERM2B(:,1:NUMALL))
         WRITE(NB6,'(/1X,"COMPLETE EXPRESSION."/)')
         CALL PRTCOR(CC(:,1:NUMALL))
         WRITE(NB6,'(//)')
      END IF
C     DEALLOCATE TEMPORARY STORAGE.
      DEALLOCATE (TERM2B)
      DEALLOCATE (SUM12A)
      DEALLOCATE (TERM1)
      DEALLOCATE (TEMP)
      DEALLOCATE (SIGMA)
      DEALLOCATE (DSMI)
      DEALLOCATE (SMI2)
      DEALLOCATE (SMI1)
      DEALLOCATE (SPLPT5)
      DEALLOCATE (SMIPT5)
      DEALLOCATE (FDERMO)
      DEALLOCATE (FDERIV)
      DEALLOCATE (FSAV)
      DEALLOCATE (ZVEC)
      DEALLOCATE (RHS)
      DEALLOCATE (ATRANS)
      DEALLOCATE (GVEC)
      DEALLOCATE (JMO)
      DEALLOCATE (IMO)
      DEALLOCATE (WMO)
      DEALLOCATE (CORSAV)
      RETURN


      CONTAINS


         SUBROUTINE MOPAIR
         IMPLICIT NONE
         DOUBLE PRECISION, PARAMETER :: TOLEMO = 1.D-4
         INTEGER NDOCC, NODD, NPAIR, MXPAIR
         IF(JODD.GT.0) THEN
            NODD = 2
         ELSE IF(IODD.GT.0) THEN
            NODD = 1
         ELSE
            NODD = 0
         END IF
         NDOCC = NUMB - NODD
         NPAIR = 0
C
C        MO PAIRS INVOLVING ONE DOUBLY OCCUPIED
C        AND ONE SEMI-OCCUPIED MO.
         IF(NODD.GT.0) THEN
            DO I=IODD,NUMB
               DO J=1,NDOCC
                  NPAIR      = NPAIR + 1
                  IMO(NPAIR) = I
                  JMO(NPAIR) = J
               END DO
            END DO
         END IF
C
C        MO PAIRS INVOLVING ONE VIRTUAL MO.
         DO I=NUMB+1,NORBS
            DO J=1,NUMB
               NPAIR      = NPAIR + 1
               IMO(NPAIR) = I
               JMO(NPAIR) = J
            END DO
         END DO
         NPAIR1 = NPAIR
         MXPAIR = LMPAIR
C
C        MO PAIRS INVOLVING TWO DOUBLY OCCUPIED MOS.
         DO I=2,NDOCC
            DO J=1,I-1
C              I.GT.J THUS E(I).GT.E(J).
               IF (E(I)-E(J).GT.TOLEMO) THEN
                  NPAIR       = NPAIR + 1
                  IMO(NPAIR)  = I
                  JMO(NPAIR)  = J
               ELSE
                  IMO(MXPAIR) = I
                  JMO(MXPAIR) = J
                  MXPAIR      = MXPAIR - 1
               END IF
            END DO
         END DO
C
C        MO PAIR INVOLVING TWO SINGLY OCCUPIED MOS.
         IF(NODD.EQ.2) THEN
            IF (E(JODD)-E(IODD).GT.TOLEMO) THEN
               NPAIR       = NPAIR + 1
               IMO(NPAIR)  = JODD
               JMO(NPAIR)  = IODD
            ELSE
               IMO(MXPAIR) = JODD
               JMO(MXPAIR) = IODD
               MXPAIR      = MXPAIR - 1
            END IF
         END IF
C
C        MO PAIRS INVOLVING TWO VIRTUAL MOS.
         DO I=NUMB+2,NORBS
            DO J=NUMB+1,I-1
C              I.GT.J THUS E(I).GT.E(J).
               IF (E(I)-E(J).GT.TOLEMO) THEN
                  NPAIR       = NPAIR + 1
                  IMO(NPAIR)  = I
                  JMO(NPAIR)  = J
               ELSE
                  IMO(MXPAIR) = I
                  JMO(MXPAIR) = J
                  MXPAIR      = MXPAIR - 1
               END IF
            END DO
         END DO
         NPAIR2 = NPAIR
         IF(NPAIR.NE.MXPAIR) THEN
            WRITE(NB6,'(1X,"INTERNAL ERROR: INCONSISTENT VALUES.")')
            WRITE(NB6,'(1X,"NPAIR1 =",I5,5X,"NPAIR2 =",I5,5X,
     1                     "MXPAIR =",I5,5X,"LMPAIR =",I5)')
     2        NPAIR1, NPAIR2, MXPAIR, LMPAIR
            STOP 'NACANA -- MOPAIR'
         END IF
         IF(DEBUG) THEN
            WRITE(NB6,'(//1X,"LIST OF MO PAIRS ",
     1                       "WITH DIFFERENT OCCUPATION."/)')
            DO IPAIR=1,NPAIR1
               WRITE(NB6,'(1X,I10,2I6)') IPAIR, IMO(IPAIR), JMO(IPAIR)
            END DO
            WRITE(NB6,'(//1X,"LIST OF NON-DEGENERATE MO PAIRS ",
     1                       "WITH EQUAL OCCUPATION."/)')
            DO IPAIR=NPAIR1+1,NPAIR2
               WRITE(NB6,'(1X,I10,2I6)') IPAIR, IMO(IPAIR), JMO(IPAIR)
            END DO
            WRITE(NB6,'(//1X,"LIST OF DEGENERATE MO PAIRS."/)')
            DO IPAIR=NPAIR2+1,LMPAIR
               WRITE(NB6,'(1X,I10,2I6)') IPAIR, IMO(IPAIR), JMO(IPAIR)
            END DO
         END IF
         END SUBROUTINE MOPAIR



         SUBROUTINE ERI2MO(CBLK,WMOBLK)
         IMPLICIT NONE
         DOUBLE PRECISION CBLK(:,:), WMOBLK(:,:)
         INTEGER          NMOPRO
         DOUBLE PRECISION PRODCC(:,:)
         ALLOCATABLE      PRODCC
         INTEGER          NMOBLK
         NMOBLK = UBOUND(CBLK,2)
         NMOPRO = NMOBLK * (NMOBLK+1) / 2
         ALLOCATE(PRODCC(LM6,NMOPRO))
         CALL MOPROD(CBLK,PRODCC)
         CALL SIMTRA(W,PRODCC,WMOBLK)
         DEALLOCATE(PRODCC)
         END SUBROUTINE ERI2MO



         SUBROUTINE MOPROD(CBLK,PRODCC)
         IMPLICIT NONE
         DOUBLE PRECISION CBLK(:,:), PRODCC(:,:)
         INTEGER          I, J, IJ, IA, MU, NU, MUNU
         INTEGER          NMOBLK, NMOPRO
         NMOBLK = UBOUND(CBLK  ,2)
         NMOPRO = UBOUND(PRODCC,2)
C        DOUBLE LOOP OVER ALL MO PAIRS.
         IJ     = 0
         DO I=1,NMOBLK
            DO J=1,I
               IJ   = IJ + 1
C              TRIPLE LOOP OVER ALL AO PAIRS PER ATOM.
               MUNU = 0
               DO IA=1,NUMAT
                  DO MU=NFIRST(IA), NLAST(IA)
                     DO NU=NFIRST(IA), MU-1
                        MUNU            = MUNU + 1
                        PRODCC(MUNU,IJ) = C(MU,I) * C(NU,J) +
     1                                    C(NU,I) * C(MU,J)
                     END DO
                     MUNU            = MUNU + 1
                     PRODCC(MUNU,IJ) = C(MU,I) * C(NU,J)
                  END DO
               END DO
            END DO
         END DO
         END SUBROUTINE MOPROD



         SUBROUTINE SIMTRA(AMAT, BMAT, CMAT)
         IMPLICIT NONE
         DOUBLE PRECISION AMAT(:,:), BMAT(:,:), CMAT(:,:)
C        CALCULATES C = BT * A * B.
         DOUBLE PRECISION TEMP(:,:)
         ALLOCATABLE      TEMP
         INTEGER          LDA, LDB, LDC, NA, NB, NC
C        DIMENSION OF MATRICES:
C          AMAT: NA X NA
C          BMAT: NA X NC
C          CMAT: NC X NC
C          TEMP: NC X NA
         LDA = UBOUND(AMAT,1)
         LDB = UBOUND(BMAT,1)
         LDC = UBOUND(CMAT,1)
         NA  = UBOUND(AMAT,2)
         NB  = UBOUND(BMAT,2)
         NC  = UBOUND(CMAT,2)
         IF(LDA.LT.NA.OR.LDC.LT.NC.OR.LDB.LT.NA.OR.NB.NE.NC) THEN
            WRITE(NB6,'(1X,"INCOMPATIBLE ARRAY DIMENSIONS.")')
            WRITE(NB6,'(1X,"  AMAT =",I6," x",I6)') LDA, NA
            WRITE(NB6,'(1X,"  BMAT =",I6," x",I6)') LDB, NB
            WRITE(NB6,'(1X,"  CMAT =",I6," x",I6)') LDC, NC
            STOP 'UANA -- SIMTRA'
         END IF
         ALLOCATE(TEMP(NC,NA))
C        TEMP = BT   * A
         CALL DGEMM('T', 'N', NC, NA, NA, 1.D0, BMAT(1,1), LDB,
     1              AMAT(1,1), LDA, 0.D0, TEMP(1,1), NC)
C        C    = TEMP * B
         CALL DGEMM('N', 'N', NC, NC, NA, 1.D0, TEMP(1,1), NC,
     1              BMAT(1,1), LDB, 0.D0, CMAT(1,1), LDC)
         DEALLOCATE(TEMP)
         END SUBROUTINE SIMTRA



         SUBROUTINE ASSERT
C ***    MAKE SURE THAT THE EQUATIONS CORRESPONDING TO
C        DEGENERATE MO PAIRS ARE FULFILLED.
         IMPLICIT NONE
         DOUBLE PRECISION, PARAMETER :: TOLRHS = 1D-8
         DOUBLE PRECISION            :: XRHS
         INTEGER                     :: IACT, JACT
         LOGICAL                     :: ERROR
C ***    INITIALIZE GVEC, NEEDED AGAIN IN MAKRHS.
C        PLEASE NOTE THAT THE ORDER OF THE MOS IN IMOCI(:)
C        IS NOT NECESSARILY ASCENDING, SO ATTENTION MUST
C        BE PAYED TO GET THE CORRECT SIGN.
         GVEC(:) = 0.D0
         DO IACT=2,NACTIV
            DO JACT=1,IACT-1
               I  = IMOCI(IACT)
               J  = IMOCI(JACT)
               IJ = KITE(I,J)
               IF(I.GT.J) THEN
                  GVEC(IJ) = ONEDEN(IACT,JACT) - ONEDEN(JACT,IACT)
               ELSE
                  GVEC(IJ) = ONEDEN(JACT,IACT) - ONEDEN(IACT,JACT)
               END IF
            END DO
         END DO
C ***    PERFORM CHECK.
         ERROR = .FALSE.
         DO IPAIR=NPAIR2+1,LMPAIR
            I    = IMO(IPAIR)
            J    = JMO(IPAIR)
            IJ   = KITE(I,J)
            XRHS = GVEC(IJ)
            IF(ABS(XRHS).GT.TOLRHS) THEN
               WRITE(NB6,'(1X,"Z VECTOR EQUATION NOT SATISFIED FOR "
     1                     "DEGENERATE MO PAIR",I4,",",I4,".")') J, I
               WRITE(NB6,'(1X,"RHS(I5) =",F14.6)')  IPAIR, XRHS
               ERROR = .TRUE.
            END IF
         END DO
         IF(ERROR) STOP 'NACANA -- ASSERT'
         END SUBROUTINE ASSERT



         SUBROUTINE MAKATR
C        FORM TOP LEFT BLOCK OF LSE COEFFICIENT MATRIX.
         IMPLICIT NONE
         DOUBLE PRECISION TEMP
         INTEGER          IPAIR,JPAIR,K,M,IJ,IK,IM,JK,JM,KM
         DO IPAIR=1,NPAIR1
            I  = IMO(IPAIR)
            J  = JMO(IPAIR)
            IJ = KITE(I,J)
            DO JPAIR=1,IPAIR
               K    = IMO(JPAIR)
               M    = JMO(JPAIR)
               IK   = KITE(I,K)
               IM   = KITE(I,M)
               JK   = KITE(J,K)
               JM   = KITE(J,M)
               KM   = KITE(K,M)
C              THE FOLLOWING EXPRESSION IS SYMMETRIC WITH RESPECT
C              TO INTERCHANGING BOTH MO PAIRS.
               TEMP = 0.5D0 * (WMO(IK,JM) + WMO(IM,JK)) -
     1                2.0D0 *  WMO(IJ,KM)
C              ONLY THE DIFFERENCE OF THE OCCUPATION NUMBERS
C              MAKES THE COEFFICIENT MATRIX NON-SYMMETRIC.
               ATRANS(IPAIR,JPAIR) = (OCC(J) - OCC(I)) * TEMP
               ATRANS(JPAIR,IPAIR) = (OCC(M) - OCC(K)) * TEMP
            END DO
            ATRANS(IPAIR,IPAIR) = ATRANS(IPAIR,IPAIR) - E(I) + E(J)
         END DO
         IF(DEBUG) THEN
            WRITE(NB6,'(//1X,"TRANSPOSED COEFFICIENT MATRIX.")')
            CALL PRTMAT(ATRANS)
         END IF
         END SUBROUTINE MAKATR



         SUBROUTINE MAKRHS
C        FORM RIGHT-HAND SIDE OF Z VECTOR LSE.
         IMPLICIT NONE
         DOUBLE PRECISION ZELEM, TEMP, ATR
         INTEGER          IPAIR,JPAIR,K,M,IJ,IK,IM,JK,JM,KM
C        INITIALIZE RIGHT-HAND SIDE VECTOR.
         DO IPAIR=1,NPAIR1
            I          = IMO(IPAIR)
            J          = JMO(IPAIR)
            IJ         = KITE(I,J)
            RHS(IPAIR) = GVEC(IJ)
         END DO
C        ELIMINATE CONTRIBUTIONS FROM MO PAIRS WITH
C        EQUAL OCCUPATION OF BOTH ORBITALS.
         ZVEC(:) = 0.D0
         DO JPAIR=NPAIR1+1,NPAIR2
            K     = IMO(JPAIR)
            M     = JMO(JPAIR)
            KM    = KITE(K,M)
            ZELEM = GVEC(KM) / (E(M) - E(K))
            DO IPAIR=1,NPAIR1
               I          = IMO(IPAIR)
               J          = JMO(IPAIR)
               IJ         = KITE(I,J)
               IK         = KITE(I,K)
               IM         = KITE(I,M)
               JK         = KITE(J,K)
               JM         = KITE(J,M)
               TEMP       = 0.5D0 * (WMO(IK,JM) + WMO(IM,JK)) -
     1                      2.0D0 *  WMO(IJ,KM)
               ATR        = (OCC(J) - OCC(I)) * TEMP
               RHS(IPAIR) = RHS(IPAIR) - ATR * ZELEM
            END DO
            ZVEC(KM) = ZELEM
         END DO
         END SUBROUTINE MAKRHS



         INTEGER FUNCTION KITE(I,J)
         IMPLICIT NONE
         INTEGER I,J,K,L
         K    = MAX(I,J)
         L    = MIN(I,J)
         KITE = K*(K-1)/2 + L
         END FUNCTION KITE



         SUBROUTINE SMICAL(SMIPT5,COORD,NORBS,DEBUG)
         USE LIMIT, ONLY: LM1
         IMPLICIT NONE
         DOUBLE PRECISION, DIMENSION(:,:) :: SMIPT5, COORD
         INTEGER                          :: NORBS
         LOGICAL                          :: DEBUG
         DOUBLE PRECISION                 :: OVER(NORBS,NORBS)
         DOUBLE PRECISION                 :: SEIGEN(NORBS), TMP
         DOUBLE PRECISION, ALLOCATABLE    :: WORK(:)
         INTEGER                          :: LWORK, INFO, I, J, K
C        FORM OVERLAP MATRIX.
         CALL OVRMAT(OVER,COORD)
         IF(DEBUG) THEN
            WRITE(NB6,'(///1X,"OVERLAP MATRIX.")')
            CALL PRTMAT(OVER)
         ENDIF
C        DIAGONALIZE OVERLAP MATRIX.
         LWORK = 3 * NORBS - 1
         ALLOCATE(WORK(LWORK))
         CALL DSYEV('V','U',NORBS,OVER(1,1),NORBS,
     1              SEIGEN(1),WORK(1),LWORK,INFO)
         IF(INFO.NE.0) THEN
            WRITE(NB6,'(1X,"DIAGONALIZATION FAILED WITH INFO=",I4)')
     1        INFO
            STOP 'SMICAL'
         END IF
         DEALLOCATE(WORK)
C        FORM SEIGEN**(-1/2).
         SEIGEN(:) = 1.0D0 / SQRT(SEIGEN(:))
C        FORM S**(-1/2).
         DO I=1,NORBS
            DO J=1,NORBS
               TMP = 0.D0
               DO K=1,NORBS
                  TMP = TMP + OVER(I,K) * OVER(J,K) * SEIGEN(K)
               END DO
               SMIPT5(I,J) = TMP
            END DO
         END DO
         IF(DEBUG) THEN
            WRITE(NB6,'(///1X,"S**(-1/2) MATRIX.")')
            CALL PRTMAT(SMIPT5)
         ENDIF
         END SUBROUTINE SMICAL



         SUBROUTINE SPLCAL(SPLPT5,COORD,NORBS,DEBUG)
         USE LIMIT, ONLY: LM1
         IMPLICIT NONE
         DOUBLE PRECISION, DIMENSION(:,:) :: SPLPT5, COORD
         INTEGER                          :: NORBS
         LOGICAL                          :: DEBUG
         DOUBLE PRECISION                 :: OVER(NORBS,NORBS)
         DOUBLE PRECISION                 :: SEIGEN(NORBS), TMP
         DOUBLE PRECISION, ALLOCATABLE    :: WORK(:)
         INTEGER                          :: LWORK, INFO, I, J, K
C        FORM OVERLAP MATRIX.
         CALL OVRMAT(OVER,COORD)
         IF(DEBUG) THEN
            WRITE(NB6,'(///1X,"OVERLAP MATRIX.")')
            CALL PRTMAT(OVER)
         ENDIF
C        DIAGONALIZE OVERLAP MATRIX.
         LWORK = 3 * NORBS - 1
         ALLOCATE(WORK(LWORK))
         CALL DSYEV('V','U',NORBS,OVER(1,1),NORBS,
     1              SEIGEN(1),WORK(1),LWORK,INFO)
         IF(INFO.NE.0) THEN
            WRITE(NB6,'(1X,"DIAGONALIZATION FAILED WITH INFO=",I4)')
     1        INFO
            STOP 'SPLCAL'
         END IF
         DEALLOCATE(WORK)
C        FORM SEIGEN**(+1/2).
         SEIGEN(:) = SQRT(SEIGEN(:))
C        FORM S**(+1/2).
         DO I=1,NORBS
            DO J=1,NORBS
               TMP = 0.D0
               DO K=1,NORBS
                  TMP = TMP + OVER(I,K) * OVER(J,K) * SEIGEN(K)
               END DO
               SPLPT5(I,J) = TMP
            END DO
         END DO
         IF(DEBUG) THEN
            WRITE(NB6,'(///1X,"S**(+1/2) MATRIX.")')
            CALL PRTMAT(SPLPT5)
         ENDIF
         END SUBROUTINE SPLCAL



         SUBROUTINE SIGNUM(CORSAV,COORD,SIGAO,IB,IXYZ)
         USE LIMIT, ONLY: LM1
         IMPLICIT NONE
         COMMON
     .   /ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
         DOUBLE PRECISION, DIMENSION(:,:) :: CORSAV, COORD, SIGAO
         INTEGER                          :: IB, IXYZ
         DOUBLE PRECISION, DIMENSION(4,4) :: OVER1, OVER2
         DOUBLE PRECISION                 :: Q0, FAC
         INTEGER                          :: NUMAT, NAT, NFIRST, NLAST
         INTEGER                          :: NORBSA, NORBSB, IA
         SIGAO(:,:) = 0.0D0
         NORBSB     = NLAST(IB) - NFIRST(IB) + 1
         FAC        = 0.5D0 / DSTEP
         Q0         = COORD(IXYZ,IB)
         DO IA=1,NUMAT
            NORBSA = NLAST(IA) - NFIRST(IA) + 1
            COORD(IXYZ,IB) = Q0 + DSTEP
            CALL OVRBLK(OVER1, CORSAV(:,IA), COORD(:,IB), IA, IB)
            COORD(IXYZ,IB) = Q0 - DSTEP
            CALL OVRBLK(OVER2, CORSAV(:,IA), COORD(:,IB), IA, IB)
            SIGAO(NFIRST(IA):NLAST(IA), NFIRST(IB):NLAST(IB)) =
     1        (OVER1(1:NORBSA, 1:NORBSB) -
     2         OVER2(1:NORBSA, 1:NORBSB)) * FAC
         END DO
         COORD(IXYZ,IB) = Q0
         END SUBROUTINE SIGNUM



         SUBROUTINE OVRMAT(OVER,COORD)
         USE LIMIT, ONLY: LM1
         IMPLICIT NONE
         COMMON
     .   /ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
         DOUBLE PRECISION, DIMENSION(:,:) :: OVER, COORD
         INTEGER                          :: NUMAT, NAT, NFIRST, NLAST
         INTEGER                          :: IA, IB
         DO IA=1,NUMAT
            DO IB=1,NUMAT
               CALL OVRBLK(OVER(NFIRST(IA):NLAST(IA),
     1                          NFIRST(IB):NLAST(IB)),
     2                     COORD(:,IA),COORD(:,IB),IA,IB)
            END DO
         END DO
         END SUBROUTINE OVRMAT



         SUBROUTINE OVRBLK(OVER,CORA,CORB,IA,IB)
         USE CONST, ONLY: OLDCF, A0
         USE LIMIT, ONLY: LM1, LMGS, LMGP
         IMPLICIT NONE
         COMMON
     .   /ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     .   /GAUSS1/ KSTART(LMGS),KNG(LMGS),KTYPE(LMGS),NSHELL,NPRIM
     .   /GAUSSP/ EXX(LMGP),DS(LMGP),DP(LMGP),DD(LMGP)
C        COMMON BLOCK GAUSS1: INFORMATION ABOUT GAUSSIAN SHELLS.
C        KSTART(I)  INDEX OF FIRST GAUSSIAN PRIMITIVE OF SHELL I.
C        KNG(I)     NUMBER OF GAUSSIAN PRIMITIVES IN SHELL I.
C        KTYPE(I)   TYPE OF SHELL I.
C        NSHELL     TOTAL NUMBER OF GAUSSIAN SHELLS.
C        NPRIM      TOTAL NUMBER OF PRIMITIVES.
C        COMMON BLOCK GAUSSP: INFORMATION ABOUT GAUSSIAN PRIMITIVES.
C        SHARED EXPONENTS AND S, P, D CONTRACTION COEFFICIENTS.
         DOUBLE PRECISION, DIMENSION(1:,1:) :: OVER
         DOUBLE PRECISION, DIMENSION(3)     :: CORA, CORB
         INTEGER                            :: IA, IB
         INTEGER                            :: NUMAT, NAT, NFIRST, NLAST
         INTEGER                            :: KSTART, KNG, KTYPE
         INTEGER                            :: NSHELL, NPRIM
         DOUBLE PRECISION                   :: EXX, DS, DP, DD
         DOUBLE PRECISION, DIMENSION(4,4)   :: PRIM
         DOUBLE PRECISION, DIMENSION(3)     :: DIST
         INTEGER                            :: NORBSA, NORBSB
         INTEGER                            :: I, J, K, L, KA, KB
         DOUBLE PRECISION                   :: AL, AL05, AL15
         DOUBLE PRECISION                   :: BE, BE05, BE15
         DOUBLE PRECISION                   :: REZI, REZI05, REZI15
         DOUBLE PRECISION                   :: REZI20, REZI25
         DOUBLE PRECISION                   :: DISTSQ, FEXP, DPROD
C
C        DESCRIPTION OF THE PARAMETERS:
C          OVER    BLOCK OF OVERLAP INTEGRALS BETWEEN
C                  BASIS FUNCTIONS ON ATOM A AND ATOM B.
C          CORA    CARTESIAN COORDINATES OF ATOM A.
C          CORB    CARTESIAN COORDINATES OF ATOM B.
C          IA      INDEX OF ATOM A IN THE LIST OF ATOMS.
C          IB      INDEX OF ATOM B IN THE LIST OF ATOMS.
C
C        LOCAL VARIABLES:
C          PRIM    OVERLAP OF GAUSSIAN PRIMITIVES.
C          DIST    INTERATOMIC DISTANCE VECTOR CORA MINUS CORB.
C          NORBSA  NUMBER OF ATOMIC ORBITALS ON ATOM A.
C          NORBSB  NUMBER OF ATOMIC ORBITALS ON ATOM B.
C
C        THIS SUBROUTINE ASSUMES THE FOLLOWING ORGANIZATION OF THE
C        GAUSSIAN BASIS SET. THERE IS ONE S SHELL FOR PER HYDROGEN
C        ATOM AND ONE SP SHELL PER FIRST-ROW ATOM. THE NUMBERING
C        OF THE SHELLS CORRESPONDS TO THE NUMBERING OF THE ATOMS.
C        THE NUMBER OF GAUSSIAN PRIMITIVES OF SHELL I (ATOM I) IS
C        KNG(I), THE INDEX OF THE FIRST PRIMITIVE IS KSTART(I).
C        EXPONENTS AND CONTRACTION COEFFICIENTS OF EACH PRIMITIVE
C        ARE STORED IN EXX, DS, AND DP, RESPECTIVELY.
C
C        INITIALIZE.
         CALL OLDCF
         OVER(:,:) = 0.0D0
         NORBSA    = NLAST(IA) - NFIRST(IA) + 1
         NORBSB    = NLAST(IB) - NFIRST(IB) + 1
C        DOUBLE LOOP OVER ALL PRIMITIVES.
         KA = KSTART(IA)
         DO K=1, KNG(IA)
            AL   = EXX(KA)
            AL05 = SQRT(AL)
            AL15 = AL * AL05
            KB = KSTART(IB)
            DO L=1, KNG(IB)
               BE        = EXX(KB)
               BE05      = SQRT(BE)
               BE15      = BE * BE05
               REZI      = 2.0D0 / (AL + BE)
               REZI05    = SQRT(REZI)
               REZI15    = REZI   * REZI05
               REZI20    = REZI   * REZI
               REZI25    = REZI20 * REZI05
               DIST(:)   = (CORA(:) - CORB(:)) / A0
               DISTSQ    = SUM(DIST(:) * DIST(:))
               FEXP      = EXP(-0.5D0*REZI*AL*BE*DISTSQ)
               PRIM(1,1) = (AL*BE)**0.75D0 * REZI15 * FEXP
               OVER(1,1) = OVER(1,1) + DS(KA) * DS(KB) * PRIM(1,1)
               IF(NORBSA.GT.1) THEN
                  DPROD  = DP(KA) * DS(KB)
                  DO I=2,4
                     PRIM(I,1) = -AL05*BE*REZI*PRIM(1,1)*DIST(I-1)
                     OVER(I,1) = OVER(I,1) + DPROD * PRIM(I,1)
                  END DO
               END IF
               IF(NORBSB.GT.1) THEN
                  DPROD  = DS(KA) * DP(KB)
                  DO I=2,4
                     PRIM(1,I) = AL*BE05*REZI*PRIM(1,1)*DIST(I-1)
                     OVER(1,I) = OVER(1,I) + DPROD * PRIM(1,I)
                  END DO
               END IF
               IF(NORBSA.GT.1.AND.NORBSB.GT.1) THEN
                  DPROD  = DP(KA) * DP(KB)
                  DO I=2,4
                     DO J=2,4
                        PRIM(I,J) = -AL15*BE15*REZI20*PRIM(1,1)*
     1                              DIST(I-1)*DIST(J-1)
                     END DO
                     PRIM(I,I)    = PRIM(I,I) +
     1                              (AL*BE)**1.25D0 * REZI25 * FEXP
                     DO J=2,4
                        OVER(I,J) = OVER(I,J) + DPROD * PRIM(I,J)
                     END DO
                  END DO
               END IF
               KB = KB + 1
            END DO
            KA = KA + 1
         END DO
         END SUBROUTINE OVRBLK




         SUBROUTINE PRTCOR(COR)
         IMPLICIT NONE
         DOUBLE PRECISION, DIMENSION(:,:) :: COR
         INTEGER                          :: I
         DO I=1, UBOUND(COR,2)
            WRITE(NB6,'(1X,I4,3F12.6)') I, COR(:,I)
         END DO
         END SUBROUTINE PRTCOR



         SUBROUTINE PRTMAT(A)
         IMPLICIT NONE
         DOUBLE PRECISION, DIMENSION(:,:) :: A
         INTEGER                          :: NCOLS, I, JSTART, JEND
         NCOLS = UBOUND(A,2)
         DO JSTART=1, NCOLS, 10
            JEND = MIN(JSTART+9, NCOLS)
            WRITE(NB6,'(/1X,10I12)') (I,I=JSTART,JEND)
            DO I=1, UBOUND(A,1)
               WRITE(NB6,'(1X,I4,10F12.8)') I, A(I,JSTART:JEND)
            END DO
         END DO
         END SUBROUTINE PRTMAT


      END SUBROUTINE NACANA
