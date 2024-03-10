      SUBROUTINE NACNUM (C,ONEDEN,ARRAY,CG,CC,DELTAE,
     1                   LM2,LM3,LM5,ICALL,SCFCAL)
C *** FULLY NUMERICAL EVALUATION OF THE SECOND TERM (MO COEFFICIENT
C     AND AO DERIVATIVES) FOR NON-ADIABATIC COUPLING ELEMENTS.
      USE LIMIT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./ATOMC / COORD(3,LM1)
     ./ATOMS / NUMAT,NAT(LMV)
     ./CIMOS / IMOCI(LMX)
     ./GUGA1 / NCIO,NCIGAM
     ./HALFE / IODD(2)
     ./INOPT2/ IN2(300)
     ./NBFILE/ NBF(20)
     ./ORBITS/ NUMB,NORBS,NMOS,NALPHA,NBETA
     ./PSDOPT/ DSTORE,DSTEP,RPSOPT(3:10),IPSOPT(28)
     ./QMMM1 / COORDM(3,LM1M),CHARGM(LM1M)
     ./QMMM6 / ISELCT(LM1+LM1M)
     ./SKIPA / ISKPA
      DIMENSION C(LM2,LM3)
      DIMENSION ONEDEN(LMACT,NCIO)
      DIMENSION CG(3,LM1+LM1M)
      DIMENSION CC(3,LM1+LM1M)
      DIMENSION ARRAY(LM5)
      LOGICAL   DEBUG
      ALLOCATABLE CORSAV(:,:), CSAV(:,:), SMIPT5(:,:), SPLPT5(:,:)
      ALLOCATABLE U1(:,:), U2(:,:), UDERIV(:,:), SMI1(:,:), SMI2(:,:)
      ALLOCATABLE DSMI(:,:) , SIGMA(:,:), TEMP(:,:)
      ALLOCATABLE TERM1(:,:), SUM12A(:,:), TERM2B(:,:)
C     FILE NUMBERS.
      NB6     = NBF(6)
C     INPUT OPTIONS.
      NUMATM = IN2(120)
      ICI1   = IN2(131)
      ICI2   = IN2(132)
      NACTIV = ICI1 + ICI2
      IF(DSTEP.EQ.0.D0)  DSTEP = 2.D-4
C     SAVE INPUT OPTIONS AND CHANGE THEM TEMPORARILY
C     TO AVOID CI CALCULATION AND UHF WHEN CALLING SCFCAL.
      IUHF    = IN2(70)
      NPRINT  = IN2(72)
      KCI     = IN2(77)
      ISKPA0  = ISKPA
      IN2(70) = -1
      IN2(77) =  0
      ISKPA   =  0
      ICALL   =  0
      DEBUG   = .FALSE.
*     DEBUG   = .TRUE.
      NUMALL  = NUMAT + NUMATM
      ALLOCATE(CORSAV(3,NUMAT))
      ALLOCATE(CSAV(NORBS,NORBS))
      ALLOCATE(SMIPT5(NORBS,NORBS))
      ALLOCATE(SPLPT5(NORBS,NORBS))
      ALLOCATE(U1(NORBS,NORBS))
      ALLOCATE(U2(NORBS,NORBS))
      ALLOCATE(UDERIV(NORBS,NORBS))
      ALLOCATE(SMI1(NORBS,NORBS))
      ALLOCATE(SMI2(NORBS,NORBS))
      ALLOCATE(DSMI(NORBS,NORBS))
      ALLOCATE(SIGMA(NORBS,NORBS))
      ALLOCATE(TEMP(NORBS,NORBS))
      ALLOCATE(TERM1(3,NUMALL))
      ALLOCATE(SUM12A(3,NUMALL))
      ALLOCATE(TERM2B(3,NUMALL))
      CORSAV(1:3,1:NUMAT)   = COORD(1:3,1:NUMAT)
      CSAV(1:NORBS,1:NORBS) = C(1:NORBS,1:NORBS)
C     DEBUG PRINT.
      IF(DEBUG) THEN
         WRITE(NB6,'(///1X,
     1         "CARTESIAN COORDINATES (ORIGINAL VALUES)."/)')
         CALL PRTCOR(CORSAV(:,1:NUMAT),NB6)
         WRITE(NB6,'(///1X,
     1         "CARTESIAN GRADIENT (ORIGINAL COORDINATES)."/)')
         CALL PRTCOR(CG(:,1:NUMAT),NB6)
         WRITE(NB6,'(///1X,"MO COEFFICIENTS (ORIGINAL COORDINATES).")')
         CALL PRTMAT(CSAV,NB6)
         WRITE(NB6,'(///1X,"ONE-PARTICLE TRANSITION DENSITY MATRIX.")')
         CALL PRTMAT(ONEDEN(1:NACTIV,1:NACTIV),NB6)
      ELSE
C        DECREASE NPRINT TO TURN OFF SCF OUTPUT.
         IN2(72) = -5
      ENDIF
C     FORM S**(-1/2) MATRIX.
      CALL SMICAL(SMIPT5,CORSAV,NORBS,DEBUG)
C     FORM S**(+1/2) MATRIX.
      CALL SPLCAL(SPLPT5,CORSAV,NORBS,DEBUG)
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
               CALL PRTCOR(COORD(:,1:NUMAT),NB6)
            ENDIF
            CALL SCFCAL(ARRAY,LM5,ICALL)
            IF(ICALL.LT.0) RETURN
            IF(DEBUG) THEN
               WRITE(NB6,'(///1X,"MO COEFFICIENTS (POSITIVE SHIFT).")')
               CALL PRTMAT(C,NB6)
            ENDIF
            CALL UMAT(CSAV,C,U1,LM2,NORBS,NB6)
            IF(DEBUG) THEN
               WRITE(NB6,'(///1X,"U MATRIX (POSITIVE SHIFT).")')
               CALL PRTMAT(U1,NB6)
            ENDIF
            CALL SMICAL(SMI1,COORD,NORBS,DEBUG)
            COORD(IC,IA) = Q0 - DSTEP
            IF(DEBUG) THEN
               WRITE(NB6,'(///1X,
     1               "CARTESIAN COORDINATES (NEGATIVE SHIFT)."/)')
               CALL PRTCOR(COORD(:,1:NUMAT),NB6)
            ENDIF
            CALL SCFCAL(ARRAY,LM5,ICALL)
            IF(ICALL.LT.0) RETURN
            IF(DEBUG) THEN
               WRITE(NB6,'(///1X,"MO COEFFICIENTS (NEGATIVE SHIFT).")')
               CALL PRTMAT(C,NB6)
            ENDIF
            CALL UMAT(CSAV,C,U2,LM2,NORBS,NB6)
            IF(DEBUG) THEN
               WRITE(NB6,'(///1X,"U MATRIX (NEGATIVE SHIFT).")')
               CALL PRTMAT(U2,NB6)
            ENDIF
            CALL SMICAL(SMI2,COORD,NORBS,DEBUG)
            COORD(IC,IA) = Q0
            FACTOR       = 0.5D0 / DSTEP
            UDERIV(:,:)  = (U1(:,:) - U2(:,:)) * FACTOR
            IF(DEBUG) THEN
               WRITE(NB6,'(///1X,"U MATRIX DERIVATIVE.")')
               CALL PRTMAT(UDERIV,NB6)
            ENDIF
            DSMI(:,:)  = (SMI1(:,:) - SMI2(:,:)) * FACTOR
            IF(DEBUG) THEN
               WRITE(NB6,'(///1X,"S**(-1/2) MATRIX DERIVATIVE.")')
               CALL PRTMAT(DSMI,NB6)
            ENDIF
C           DIVIDE BY ENERGY DIFFERENCE.
            TERM1(IC,IA) = CG(IC,IA) / DELTAE
C           NO CONTRIBUTION FROM THE SCF REFERENCE DUE TO ANTISYMMETRY
C           OF THE U MATRIX DERIVATIVE (ZERO DIAGONAL).
            TMP = 0.D0
            DO I=1,NACTIV
               DO J=1,NACTIV
                  TMP = TMP + ONEDEN(I,J) * UDERIV(IMOCI(I), IMOCI(J))
               END DO
            END DO
            SUM12A(IC,IA) = TERM1(IC,IA) + TMP
C           FORM SIGMA MATRIX IN NON-ORTHOGONAL AO BASIS.
            CALL SIGNUM(CORSAV,COORD,SIGMA,IA,IC)
            IF(DEBUG) THEN
               WRITE(NB6,'(///1X,"AO BASIS SIGMA MATRIX.")')
               CALL PRTMAT(SIGMA,NB6)
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
               CALL PRTMAT(SIGMA,NB6)
            ENDIF
C           FORM SIGMA = SIGMA + SPLPT5 * DSMI.
            CALL DGEMM('N','N',NORBS,NORBS,NORBS,1.D0,SPLPT5(1,1),NORBS,
     1                 DSMI(1,1),NORBS,1.D0,SIGMA(1,1),NORBS)
            IF(DEBUG) THEN
               WRITE(NB6,'(///1X,
     1               "ANTISYMMETRIC AO BASIS SIGMA MATRIX.")')
               CALL PRTMAT(SIGMA,NB6)
            ENDIF
C           TRANSFORM SIGMA MATRIX INTO MO BASIS.
C           TEMP = CSAV(TRANS) * SIGMA.
            CALL DGEMM('T','N',NORBS,NORBS,NORBS,1.D0,CSAV(1,1),NORBS,
     1                 SIGMA(1,1),NORBS,0.D0,TEMP(1,1),NORBS)
C           SIGMA = TEMP * CSAV.
            CALL DGEMM('N','N',NORBS,NORBS,NORBS,1.D0,TEMP(1,1),NORBS,
     1                 CSAV(1,1),NORBS,0.D0,SIGMA(1,1),NORBS)
            IF(DEBUG) THEN
               WRITE(NB6,'(///1X,"MO BASIS SIGMA MATRIX.")')
               CALL PRTMAT(SIGMA,NB6)
            ENDIF
C           SIGMA CONTRIBUTIONS TO THE NON-ADIABATIC COUPLING.
            TMP = 0.D0
C           NO DIAGONAL CONTRIBUTION FOR ANTISYMMETRIC MATRIX.
C           DO I=1,NUMB
C              TMP = TMP + 2.D0 * SIGMA(I,I)
C           END DO
C           IF(IODD(1).GT.0)  TMP = TMP - SIGMA(IODD(1), IODD(1))
C           IF(IODD(2).GT.0)  TMP = TMP - SIGMA(IODD(2), IODD(2))
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
                  CALL PRTCOR(COORDM(:,1:NUMATM),NB6)
               ENDIF
               CALL SCFCAL(ARRAY,LM5,ICALL)
               IF(ICALL.LT.0) RETURN
               IF(DEBUG) THEN
                  WRITE(NB6,'(///1X,
     1                  "MO COEFFICIENTS (POSITIVE SHIFT).")')
                  CALL PRTMAT(C,NB6)
               ENDIF
               CALL UMAT(CSAV,C,U1,LM2,NORBS,NB6)
               IF(DEBUG) THEN
                  WRITE(NB6,'(///1X,"U MATRIX (POSITIVE SHIFT).")')
                  CALL PRTMAT(U1,NB6)
               ENDIF
               COORDM(IC,IAMM) = Q0 - DSTEP
               IF(DEBUG) THEN
                  WRITE(NB6,'(///1X,
     1                  "CARTESIAN COORDINATES (NEGATIVE SHIFT)."/)')
                  CALL PRTCOR(COORD(:,1:NUMAT),NB6)
               ENDIF
               CALL SCFCAL(ARRAY,LM5,ICALL)
               IF(ICALL.LT.0) RETURN
               IF(DEBUG) THEN
                  WRITE(NB6,'(///1X,
     1                  "MO COEFFICIENTS (NEGATIVE SHIFT).")')
                  CALL PRTMAT(C,NB6)
               ENDIF
               CALL UMAT(CSAV,C,U2,LM2,NORBS,NB6)
               IF(DEBUG) THEN
                  WRITE(NB6,'(///1X,"U MATRIX (NEGATIVE SHIFT).")')
                  CALL PRTMAT(U2,NB6)
               ENDIF
               COORDM(IC,IAMM) = Q0
               FACTOR          = 0.5D0 / DSTEP
               UDERIV(:,:)     = (U1(:,:) - U2(:,:)) * FACTOR
               IF(DEBUG) THEN
                  WRITE(NB6,'(///1X,"U MATRIX DERIVATIVE.")')
                  CALL PRTMAT(UDERIV,NB6)
               ENDIF
C              DIVIDE BY ENERGY DIFFERENCE.
               TERM1(IC,IA) = CG(IC,IA) / DELTAE
C              NO CONTRIBUTION FROM THE SCF REFERENCE DUE TO ANTISYMMETRY
C              OF THE U MATRIX DERIVATIVE (ZERO DIAGONAL).
               TMP = 0.D0
               DO I=1,NACTIV
                  DO J=1,NACTIV
                     TMP = TMP + ONEDEN(I,J) *
     1                           UDERIV(IMOCI(I), IMOCI(J))
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
C     RESTORE ALL QUANTITIES.
      CALL SCFCAL(ARRAY,LM5,ICALL)
      IF(ICALL.LT.0) RETURN
C     PRINT RESULTS.
      WRITE(NB6,'(///1X,"NON-ADIABATIC COUPLING MATRIX ELEMENTS.")')
      WRITE(NB6,'(/1X,"FIRST TERM."/)')
      CALL PRTCOR(TERM1 (:,1:NUMALL),NB6)
      WRITE(NB6,'(/1X,"FIRST TERM PLUS SECOND TERM, FIRST PART."/)')
      CALL PRTCOR(SUM12A(:,1:NUMALL),NB6)
      WRITE(NB6,'(/1X,"SECOND TERM, SECOND PART."/)')
      CALL PRTCOR(TERM2B(:,1:NUMALL),NB6)
      WRITE(NB6,'(/1X,"COMPLETE EXPRESSION."/)')
      CALL PRTCOR(CC(:,1:NUMALL),NB6)
      WRITE(NB6,'(//)')
C     DEALLOCATE TEMPORARY STORAGE.
      DEALLOCATE (TERM2B)
      DEALLOCATE (SUM12A)
      DEALLOCATE (TERM1)
      DEALLOCATE (TEMP)
      DEALLOCATE (SIGMA)
      DEALLOCATE (DSMI)
      DEALLOCATE (SMI2)
      DEALLOCATE (SMI1)
      DEALLOCATE (UDERIV)
      DEALLOCATE (U2)
      DEALLOCATE (U1)
      DEALLOCATE (SPLPT5)
      DEALLOCATE (SMIPT5)
      DEALLOCATE (CSAV)
      DEALLOCATE (CORSAV)
C     RESTORE INPUT OPTIONS.
      IN2(70) = IUHF
      IN2(72) = NPRINT
      IN2(77) = KCI
      ISKPA   = ISKPA0
      RETURN


      CONTAINS


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
            CALL PRTMAT(OVER,NB6)
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
            CALL PRTMAT(SMIPT5,NB6)
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
            CALL PRTMAT(OVER,NB6)
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
            CALL PRTMAT(SPLPT5,NB6)
         ENDIF
         END SUBROUTINE SPLCAL



         SUBROUTINE UMAT(CSAV,C,U,LDC,N,NB6)
         IMPLICIT NONE
         DOUBLE PRECISION, DIMENSION(:,:) :: CSAV, C, U
         INTEGER                          :: LDC, N, NB6
         DOUBLE PRECISION, PARAMETER      :: THRESH = 0.9D0
         DOUBLE PRECISION                 :: H
         INTEGER                          :: I, J, K
         CALL DGEMM('T', 'N', N, N, N, 1.D0,
     1              CSAV(1,1), N, C(1,1), LDC, 0.D0, U(1,1), N)
         DO I=1,N
C           FIND MO CORRESPONDING TO ORIGINAL MO I.
            K = 1
            DO J=2,N
               IF(ABS(U(I,J)).GT.ABS(U(I,K)))  K = J;
            END DO

            IF(ABS(U(I,K)).LT.THRESH) THEN
               WRITE(NB6,'(1X,"MAPPING OF MOS FAILED.")')
               WRITE(NB6,'(1X,"LARGEST OVERLAP:",F10.6)')  ABS(U(I,K))
               STOP 'NACNUM -- UMAT'
            END IF

            IF(U(I,K).LT.0.D0)  U(:,K) = -U(:,K)

            IF(K.NE.I) THEN
C              SWITCH MOS.
               DO J=1,N
                  H      = U(J,I)
                  U(J,I) = U(J,K)
                  U(J,K) = H
               END DO
            END IF
         END DO
         END SUBROUTINE UMAT



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




         SUBROUTINE PRTCOR(COORD,NB6)
         IMPLICIT NONE
         DOUBLE PRECISION, DIMENSION(:,:) :: COORD
         INTEGER                          :: NB6, I
         DO I=1, UBOUND(COORD,2)
            WRITE(NB6,'(1X,I4,3F12.6)') I, COORD(:,I)
         END DO
         END SUBROUTINE PRTCOR



         SUBROUTINE PRTMAT(C,NB6)
         IMPLICIT NONE
         DOUBLE PRECISION, DIMENSION(:,:) :: C
         INTEGER                          :: NB6
         INTEGER                          :: NCOLS, I, JSTART, JEND
         NCOLS = UBOUND(C,2)
         DO JSTART=1, NCOLS, 10
            JEND = MIN(JSTART+9, NCOLS)
            WRITE(NB6,'(/1X,10I12)') (I,I=JSTART,JEND)
            DO I=1, UBOUND(C,1)
               WRITE(NB6,'(1X,I4,10F12.8)') I, C(I,JSTART:JEND)
            END DO
         END DO
         END SUBROUTINE PRTMAT


      END SUBROUTINE NACNUM
