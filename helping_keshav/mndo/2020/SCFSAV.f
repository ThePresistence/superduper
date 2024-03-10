      SUBROUTINE SCFSAV (NSAV13,NSAV15,ICALL1,ICALL3,ICALL4)
C     *
C     SAVE CURRENT RESULTS IN SUBROUTINE SCF.
C     USEFUL AS INTERFACE TO OTHER PROGRAMS AND FOR DEBUGGING.
C     *
C     NSAV13=2: SAVE DATA FOR POSTPROCESSING USING MOLDEN.
C     NSAV15>0: SAVE DATA FOR PROCESSING IN QM/MM CODES.
C     *
C     NOTATION. I=INPUT, O=OUTPUT.
C     NSAV13    CONTROL FLAG FOR SAVING MOLDEN DATA ON FILE NB13 (I).
C               = 0  DO NOT SAVE SUCH DATA.
C               = 2  SAVE CARTESIAN COORDINATES AND GRADIENT NORMS.
C     NSAV15    CONTROL OVER WHICH DATA ARE SAVED ON FILE NB15 (I).
C               = 0  DO NOT SAVE ANYTHING.
C               = 1  SAVE CARTESIAN COORDINATES.
C               = 2  ALSO SAVE ENERGY AND GRADIENT NORMS.
C               = 3  ALSO SAVE CARTESIAN AND INTERNAL GRADIENT.
C               = 4  ALSO SAVE ELECTROSTATIC POTENTIAL AND FIELD.
C               = 5  EQUIVALENT TO 4, BUT WRITES UNFORMATTED FILE.
C               = 9  SAVE EVERYTHING AND DO NOT REWIND FILE NB15.
C     ICALL3    AVAILABILITY OF CARTESIAN GRADIENT (I).
C               = 0  NOT AVAILABLE.
C               = 1  AVAILABLE
C     ICALL4    AVAILABILITY OF INTERNAL GRADIENT (I).
C               = 0  NOT AVAILABLE.
C               = 1  AVAILABLE
C     *
      USE LIMIT, ONLY: LM1, LMV, LM1M
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*2  ELEMNT
      CHARACTER*80 KTITLE,KOMENT
      COMMON
     ./ATOMC / COORD(3,LM1)
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./CGRAD / CG(3,LM1+LM1M)
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./DFP   / X(LMV),NVAR
     ./ELEMTS/ ELEMNT(107)
     ./ERG   / ENERGY,G(LMV),GNORM,CNORM
     ./FLAG1 / KTITLE,KOMENT
     ./INOPT2/ IN2(300)
     ./NBFILE/ NBF(20)
     ./PARM3 / LOC(LMV),NVAR3
     ./QMMM1 / COORDM(3,LM1M),CHARGM(LM1M)
     ./QMMM2 / LINK(LM1)
     ./QMMM4 / ELPOT(LM1M),ESF(3,LM1M)
C
C *** SAVE DATA FOR POSTPROCESSING USING MOLDEN.
      IF(NSAV13.NE.2) GO TO 100
C     SKIP MOLDEN OUTPUT DURING FORCE CONSTANT CALCULATION.
      IF(ICALL1.EQ.3) GO TO 100
      NB13 = NBF(13)
C     CARTESIAN GRADIENT NORM.
      IF(ICALL3.EQ.1) THEN
         GGMAX  = ZERO
         SUM    = ZERO
         DO 20 I=1,NUMAT
         DO 10 J=1,3
         GGMAX  = MAX(GGMAX,ABS(CG(J,I)))
         SUM    = SUM + CG(J,I)*CG(J,I)
   10    CONTINUE
   20    CONTINUE
         GGRMS = SQRT(SUM) / (THREE*NUMAT)
         WRITE(NB13,400) '[GEOCONV]  CG', ENERGY, GGRMS, GGMAX
         WRITE(NB13,410) '[GEOMETRIES]  CG   XYZ', NUMAT
         WRITE(NB13,400) KTITLE
         DO 30 I=1,NUMAT
         WRITE(NB13,420) ELEMNT(NAT(I)),(COORD(J,I),J=1,3)
   30    CONTINUE
      ENDIF
C     INTERNAL GRADIENT NORM.
      IF(ICALL4.EQ.1) THEN
         GGMAX  = ZERO
         SUM    = ZERO
         DO 40 I=1,NVAR
         GGMAX  = MAX(GGMAX,ABS(G(I)))
         SUM    = SUM + G(I)*G(I)
   40    CONTINUE
         GGRMS = SQRT(SUM) / NVAR
         WRITE(NB13,400) '[GEOCONV]  G ', ENERGY, GGRMS, GGMAX
         WRITE(NB13,410) '[GEOMETRIES]  G    XYZ', NUMAT
         WRITE(NB13,400) KTITLE
         DO 50 I=1,NUMAT
         WRITE(NB13,420) ELEMNT(NAT(I)),(COORD(J,I),J=1,3)
   50    CONTINUE
      ENDIF
C
C *** SAVE DATA FOR PROCESSING IN OTHER CODES.
  100 CONTINUE
      IF(NSAV15.LE.0) RETURN
C     INPUT OPTIONS.
      MMINP  = IN2(12)
      NUMATM = IN2(120)
      MMCOUP = IN2(121)
      MMPOT  = IN2(122)
      NB15   = NBF(15)

      IF (NSAV15.EQ.5) THEN
         CALL write_bin
         RETURN
      END IF

      IF(NSAV15.LT.9) REWIND NB15
C *** SAVE CURRENT CARTESIAN COORDINATES.
C     REAL QM ATOMS.
      WRITE(NB15,500) NUMAT
      DO 110 I=1,NUMAT
      WRITE(NB15,510) I,NAT(I),(COORD(J,I),J=1,3)
  110 CONTINUE
      WRITE(NB15,520)
C     EXTERNAL POINTS OR MM ATOMS.
      IF(MMINP.GT.0) THEN
         WRITE(NB15,530) NUMATM
         DO 120 I=1,NUMATM
         K = NUMAT+I
         WRITE(NB15,510) K,I,(COORDM(J,I),J=1,3),CHARGM(I)
  120    CONTINUE
         WRITE(NB15,520)
      ENDIF
      IF(NSAV15.LE.1) RETURN
C *** SAVE ENERGY AND GRADIENT NORMS (ALWAYS DEFINED).
      WRITE(NB15,540) ENERGY,CNORM,GNORM
      WRITE(NB15,520)
      IF(NSAV15.LE.2) RETURN
C *** SAVE CARTESIAN GRADIENT (IF DEFINED).
      IF(ICALL3.EQ.1) THEN
C        REAL QM ATOMS.
         IF(MMINP.LE.0) THEN
            WRITE(NB15,550)
            DO 130 I=1,NUMAT
            WRITE(NB15,510) I,NAT(I),(CG(J,I),J=1,3)
  130       CONTINUE
            WRITE(NB15,520)
C        REAL QM ATOMS AND EXTERNAL POINTS OR MM ATOMS.
         ELSE
            WRITE(NB15,550)
            DO 140 I=1,NUMAT
            WRITE(NB15,560) I,NAT(I),(CG(J,I),J=1,3),LINK(I)
  140       CONTINUE
            WRITE(NB15,520)
            WRITE(NB15,570)
            DO 150 I=1,NUMATM
            K = NUMAT+I
            WRITE(NB15,510) K,I,(CG(J,K),J=1,3)
  150       CONTINUE
            WRITE(NB15,520)
         ENDIF
      ENDIF
C *** SAVE INTERNAL COORDINATES AND GRADIENTS (IF DEFINED).
      IF(ICALL4.EQ.1) THEN
         WRITE(NB15,580) NVAR
         DO 160 I=1,NVAR
         K   = (LOC(I)+2)/3
         L   = LOC(I)-3*(K-1)
         WRITE(NB15,510) K,L,X(I),G(I)
  160    CONTINUE
         WRITE(NB15,520)
      ENDIF
      IF(NSAV15.LE.3) RETURN
C *** SAVE ELECTROSTATIC POTENTIAL.
      IF(MMINP.EQ.1) THEN
         WRITE(NB15,590) MMPOT
         DO 170 I=1,NUMATM
         K = NUMAT+I
         WRITE(NB15,510) K,I,ELPOT(I)
  170    CONTINUE
         WRITE(NB15,520)
      ENDIF
C *** SAVE ELECTRIC FIELD.
      IF(MMINP.GT.0 .AND. MMCOUP.EQ.3) THEN
         WRITE(NB15,600) MMPOT
         DO 180 I=1,NUMATM
         K = NUMAT+I
         WRITE(NB15,510) K,I,(ESF(J,I),J=1,3)
  180    CONTINUE
         WRITE(NB15,520)
      ENDIF
      RETURN

  400 FORMAT(A,3F20.10)
  410 FORMAT(A,I5)
  420 FORMAT(A,5X,3F20.10)
  500 FORMAT(' CARTESIAN COORDINATES: NUMAT =',I5)
  510 FORMAT(2I5,3F20.10,F10.5)
  520 FORMAT(' ')
  530 FORMAT(' CARTESIAN COORDINATES AND CHARGES OF MM ATOMS:',
     1       '  NUMATM =',I5)
  540 FORMAT(' ENERGY, CARTESIAN AND INTERNAL GRADIENT NORM',
     1       / 10X,3F20.10)
  550 FORMAT(' CARTESIAN GRADIENT')
  560 FORMAT(2I5,3F20.10,I5)
  570 FORMAT(' CARTESIAN GRADIENT OF MM ATOMS')
  580 FORMAT(' INTERNAL COORDINATES AND GRADIENTS: NVAR =',I5)
  590 FORMAT(' ELECTROSTATIC POTENTIAL AT THE MM ATOMS: MMPOT =',I5)
  600 FORMAT(' ELECTRIC FIELD AT THE MM ATOMS: MMPOT =',I5)

      CONTAINS

      SUBROUTINE write_bin
       !---------------------------------------------------------------------
       ! Write unformatted (binary) file for data transfer to ChemShell
       ! Data written corresponds to formatted file obtained with NSAV=4
       ! Format of the file:
       ! I:    Length of content vector (= number of sections in file)
       ! I(n): Vector containing number of items in each section
       ! R:    n sections of length specified in content vector
       !       All items are reals
       !---------------------------------------------------------------------
        IMPLICIT NONE
        INTEGER(4),    PARAMETER      :: Nsec = 7 ! # of sections in file
        CHARACTER(32), PARAMETER      :: fname = 'fort.15'
        INTEGER(4)                    :: items(Nsec) = 0
        REAL(8),       ALLOCATABLE    :: MMdata(:,:)
        !--------------------------------------------------------------------
        OPEN(UNIT=NB15, FILE=TRIM(fname), STATUS='REPLACE',
     $        FORM='UNFORMATTED', ACTION='WRITE')

        ! Construct and write content vector
        items(1) = 3*NUMAT                  ! Coords of QM atoms
        IF (MMINP.GT.0) items(2) = 4*NUMATM ! Coords and charge of MM atoms
        items(3) = 3                        ! Energy, cart. and internal grad norm
        IF (ICALL3.EQ.1) THEN               ! Cart. gradient
           IF (MMINP.LE.0) THEN
              items(4) = 3*NUMAT
           ELSE
              items(4) = 3*(NUMAT + NUMATM)
           END IF
        END IF
        IF (ICALL4.EQ.1) items(5) = 2*NVAR  ! Internal coords and gradient
        IF (MMINP.EQ.1)  items(6) = NUMATM  ! El.stat. potential at MM atoms
        IF (MMINP.GT.0 .AND. MMCOUP.EQ.3) items(7) = 3*NUMATM  ! Electric field at MM atoms

        WRITE(NB15) Nsec
        WRITE(NB15) items

        ! SAVE CURRENT CARTESIAN COORDINATES OF REAL QM ATOMS
        ! N.B.: We have to specify upper bounds explicitly because arrays are defined larger
        !       than effectively used.
        WRITE(NB15) COORD(:,1:NUMAT)

        ! SAVE COORDS AND CHARGES OF EXTERNAL POINTS OR MM ATOMS
        IF (MMINP.GT.0) THEN
           ALLOCATE(MMdata(4,NUMATM))
           MMdata(1:3,:) = COORDM(:,1:NUMATM)
           MMdata(4,:)   = CHARGM(1:NUMATM)
           WRITE(NB15) MMdata
           DEALLOCATE(MMdata)
        END IF

        ! SAVE ENERGY AND GRADIENT NORMS (ALWAYS DEFINED)
        WRITE(NB15) ENERGY, CNORM, GNORM

        ! SAVE CARTESIAN GRADIENT (IF DEFINED)
        ! CG is the full grad matrix, first QM, then MM
        IF (ICALL3.EQ.1) WRITE(NB15) CG(:,1:NUMAT+NUMATM)

        ! SAVE INTERNAL COORDINATES AND GRADIENTS (IF DEFINED).
        IF (ICALL4.EQ.1) WRITE(NB15) X(1:NVAR), G(1:NVAR)

        ! SAVE ELECTROSTATIC POTENTIAL
        IF (MMINP.EQ.1) WRITE(NB15) ELPOT(1:NUMATM)

        ! SAVE ELECTRIC FIELD.
        IF(MMINP.GT.0 .AND. MMCOUP.EQ.3) WRITE(NB15) ESF(:,NUMATM)

        CLOSE(nb15)
        RETURN
      END SUBROUTINE write_bin

      END SUBROUTINE SCFSAV
