      SUBROUTINE MULSAV (NSAV13,NSAV15,ICALL3,ICALL4)
C     *
C     SAVE CURRENT RESULTS IN SUBROUTINE MULSAV.
C     ANALOGOUS TO SCFSAV FOR SINGLE SURFACE CALCULATIONS.
C     *
C     NSAV13=2: SAVE DATA FOR POSTPROCESSING USING MOLDEN.
C     NSAV15>0: SAVE DATA FOR PROCESSING IN QM/MM CODES. 
C     *
C     NOTATION. I=INPUT, O=OUTPUT.
C     NSAV13    CONTROL FLAG FOR SAVING MOLDEN DATA ON FILE NB13 (I).
C               = 0  DO NOT SAVE SUCH DATA.
C               = 2  SAVE CARTESIAN COORDINATES AND GRADIENT NORMS.
C                    THE NORMS SAVED ARE THOSE OF THE STATE OF
C                    INTEREST DEFINED BY LROOT, PROVIDED THAT LROOT IS
C                    IN THE LIST OF NCIGRD GRADIENTS.
C     NSAV15    CONTROL OVER WHICH DATA ARE SAVED ON FILE NB15 (I).
C               = 0  DO NOT SAVE ANYTHING.
C               = 1  SAVE CARTESIAN COORDINATES.
C               = 2  ALSO SAVE MULTIPLE STATE ENERGIES AND GRADIENT
C                    NORMS.
C               = 3  ALSO SAVE CARTESIAN AND INTERNAL GRADIENTS FOR
C                    EACH STATE AND INTERSTATE COUPLING GRADIENTS IF
C                    AVAILABLE.
C               = 4  EQUIVALENT TO 3, AS THE ELECTROSTATIC POTENTIAL
C                    AND FIELD CANNOT CURRENTLY BE CALCULATED USING
C                    THE GUGA-CI METHOD.
C               = 5  EQUIVALENT TO 4, BUT WRITES UNFORMATTED FILE.
C               = 6  EQUIVALENT TO 4, BUT ALSO SAVES CI COMPONENTS
C                    FOR SIMPLIFIED SURFACE HOPPING ETC.
C               = 7  EQUIVALENT TO 6, BUT WRITES UNFORMATTED FILE.
C               = 9  SAVE EVERYTHING AND DO NOT REWIND FILE NB15.
C     ICALL3    AVAILABILITY OF CARTESIAN GRADIENT (I).
C               = 0  NOT AVAILABLE.
C               = 1  MULTIPLE STATE GRADIENTS AVAILABLE
C               = 2  MULTIPLE STATE GRADIENTS AND INTERSTATE COUPLING
C                    GRADIENTS AVAILABLE
C     ICALL4    AVAILABILITY OF INTERNAL GRADIENT (I).
C               = 0  NOT AVAILABLE.
C               = 1  MULTIPLE STATE GRADIENTS AVAILABLE
C               = 2  MULTIPLE STATE GRADIENTS AND INTERSTATE COUPLING
C                    GRADIENTS AVAILABLE
C     *
      USE LIMIT, ONLY: LM1, LMV, LM1M, LMGRD, LMNAC, LMCONF
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*2  ELEMNT
      CHARACTER*80 KTITLE,KOMENT
      COMMON
     ./ATOMC / COORD(3,LM1)
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./CGRAD / CG(3,LM1+LM1M)
     ./CGRADM/ CGX(3,LM1+LM1M,LMGRD),CNORMX(LMGRD)
     ./CNACVM/ CCX(3,LM1+LM1M,LMNAC),CCNRMX(LMNAC)
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./DFP   / X(LMV),NVAR
     ./DYNVA/  CICOMP(LMCONF,6),NCICONF
     ./ELEMTS/ ELEMNT(107)
     ./ERG   / ENERGY,G(LMV),GNORM,CNORM
     ./ERGM  / ENERGX(LMGRD),GRADX(LMV,LMGRD),GNORMX(LMGRD)
     ./FLAG1 / KTITLE,KOMENT
     ./GNACVM/ GCX(LMV,LMNAC),GCNRMX(LMNAC)
     ./GRDORG/ IGRST(LMGRD),ISTATE,JSTATE
     ./INOPT2/ IN2(300)
     ./NBFILE/ NBF(20)
     ./PARM3 / LOC(LMV),NVAR3
     ./QMMM1 / COORDM(3,LM1M),CHARGM(LM1M)
     ./QMMM2 / LINK(LM1)
     ./QMMM4 / ELPOT(LM1M),ESF(3,LM1M)
C
C *** SAVE DATA FOR POSTPROCESSING USING MOLDEN.
C     NO CHANGE FROM SCFSAV VERSION
      IF(NSAV13.NE.2) GO TO 100
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
      NCIGRD = IN2(159)
      NB15   = NBF(15)

      IF (NSAV15.EQ.5 .OR. NSAV15.EQ.7) THEN
         CALL write_multistate_bin
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
C *** SAVE ELECTRONIC STATES, ENERGIES, AND GRADIENT NORMS.
      WRITE(NB15,540) NCIGRD
      DO 122 I=1,NCIGRD
      WRITE(NB15,542) IGRST(I),ENERGX(I),CNORMX(I),GNORMX(I)
 122  CONTINUE
      WRITE(NB15,520)
      IF(NSAV15.LE.2) RETURN
C *** SAVE CARTESIAN GRADIENTS (IF DEFINED).
      IF(ICALL3.GE.1) THEN
C        LOOP OVER STATES
         DO 152 IS=1,NCIGRD
C        REAL QM ATOMS.
         IF(MMINP.LE.0) THEN
            WRITE(NB15,550) IGRST(IS)
            DO 130 I=1,NUMAT
            WRITE(NB15,510) I,NAT(I),(CGX(J,I,IS),J=1,3)
  130       CONTINUE
            WRITE(NB15,520)
C        REAL QM ATOMS AND EXTERNAL POINTS OR MM ATOMS.
         ELSE
            WRITE(NB15,550) IGRST(IS)
            DO 140 I=1,NUMAT
            WRITE(NB15,560) I,NAT(I),(CGX(J,I,IS),J=1,3),LINK(I)
  140       CONTINUE
            WRITE(NB15,520)
            WRITE(NB15,570) IGRST(IS)
            DO 150 I=1,NUMATM
            K = NUMAT+I
            WRITE(NB15,510) K,I,(CGX(J,K,IS),J=1,3)
  150       CONTINUE
            WRITE(NB15,520)
         ENDIF
  152    CONTINUE
      ENDIF
C *** SAVE INTERNAL COORDINATES AND GRADIENTS (IF DEFINED).
      IF(ICALL4.GE.1) THEN
C        FIRST, INTERNAL COORDINATES
         WRITE(NB15,580) NVAR
         DO 160 I=1,NVAR
         K   = (LOC(I)+2)/3
         L   = LOC(I)-3*(K-1)
         WRITE(NB15,510) K,L,X(I)
  160    CONTINUE
         WRITE(NB15,520)
C        LOOP OVER STATES
         DO 164 IS=1,NCIGRD
         WRITE(NB15,582) IGRST(IS)
         DO 162 I=1,NVAR
         K   = (LOC(I)+2)/3
         L   = LOC(I)-3*(K-1)
         WRITE(NB15,510) K,L,GRADX(I,IS)
  162    CONTINUE
         WRITE(NB15,520)
  164    CONTINUE  
      ENDIF
C *** SAVE INTERSTATE COUPLING NORMS (IF DEFINED).
      IF(ICALL3.EQ.2 .OR. ICALL4.EQ.2) THEN
         WRITE(NB15,700)
         DO 210 IS=2,NCIGRD
            DO 200 JS=1,IS-1
               IJ = (((IS-1)*(IS-2))/2) + JS
               WRITE(NB15,710) IGRST(IS), IGRST(JS),
     1                         CCNRMX(IJ), GCNRMX(IJ)
 200        CONTINUE
 210     CONTINUE
         WRITE(NB15,520)
      ENDIF
C *** SAVE CARTESIAN INTERSTATE COUPLING GRADIENTS (IF DEFINED).
      IF(ICALL3.EQ.2) THEN
C        LOOP OVER STATES
         DO 260 IS=2,NCIGRD
         DO 250 JS=1,IS-1
         IJ = (((IS-1)*(IS-2))/2) + JS
C        REAL QM ATOMS.
         IF(MMINP.LE.0) THEN
            WRITE(NB15,720) IGRST(IS), IGRST(JS)
            DO 220 I=1,NUMAT
            WRITE(NB15,510) I,NAT(I),(CCX(J,I,IJ),J=1,3)
  220       CONTINUE
            WRITE(NB15,520)
C        REAL QM ATOMS AND EXTERNAL POINTS OR MM ATOMS.
         ELSE
            WRITE(NB15,720) IGRST(IS), IGRST(JS)
            DO 230 I=1,NUMAT
            WRITE(NB15,560) I,NAT(I),(CCX(J,I,IJ),J=1,3),LINK(I)
  230       CONTINUE
            WRITE(NB15,520)
            WRITE(NB15,730) IGRST(IS), IGRST(JS)
            DO 240 I=1,NUMATM
            K = NUMAT+I
            WRITE(NB15,510) K,I,(CCX(J,K,IJ),J=1,3)
  240       CONTINUE
            WRITE(NB15,520)
         ENDIF
  250    CONTINUE
  260    CONTINUE
      ENDIF
C *** SAVE INTERNAL INTERSTATE COUPLING GRADIENTS (IF DEFINED).
      IF(ICALL4.EQ.2) THEN
C        LOOP OVER STATES
         DO 290 IS=2,NCIGRD
         DO 280 JS=1,IS-1
         IJ = (((IS-1)*(IS-2))/2) + JS
         WRITE(NB15,740) IGRST(IS), IGRST(JS)
         DO 270 I=1,NVAR
         K   = (LOC(I)+2)/3
         L   = LOC(I)-3*(K-1)
         WRITE(NB15,510) K,L,GCX(I,IJ)
  270    CONTINUE
         WRITE(NB15,520)
  280    CONTINUE 
  290    CONTINUE
      ENDIF
C *** SAVE CI COMPONENTS
      IF(NSAV15.EQ.6) THEN
         WRITE(NB15,800) NCICONF
C        LOOP OVER STATES
         DO 300 IS=1,NCIGRD
            WRITE(NB15,810) IGRST(IS)
            DO 310 I=1,NCICONF
               WRITE(NB15,820) I,CICOMP(I,IGRST(IS))
  310       CONTINUE    
            WRITE(NB15,520)
  300    CONTINUE 
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
  540 FORMAT(' STATES, ENERGIES, CARTESIAN AND',
     1       ' INTERNAL GRADIENT NORMS: NCIGRD =',I5)
  542 FORMAT(I5,5X,3F20.10)
  550 FORMAT(' CARTESIAN GRADIENT FOR STATE',I5)
  560 FORMAT(2I5,3F20.10,I5)
  570 FORMAT(' CARTESIAN GRADIENT OF MM ATOMS FOR STATE',I5)
  580 FORMAT(' INTERNAL COORDINATES: NVAR =',I5)
  582 FORMAT(' INTERNAL GRADIENT FOR STATE',I5) 
  700 FORMAT(' CARTESIAN AND INTERNAL INTERSTATE COUPLING NORMS')
  710 FORMAT(2I5,20X,2F20.10) 
  720 FORMAT(' CARTESIAN INTERSTATE COUPLING GRADIENT FOR STATES',2I5)
  730 FORMAT(' CARTESIAN INTERSTATE COUPLING GRADIENT OF MM ATOMS',
     1       ' FOR STATES ',2I5)   
  740 FORMAT(' INTERNAL INTERSTATE COUPLING GRADIENT FOR STATES',2I5) 
  800 FORMAT(' CI COMPONENTS: NCICONF =',I5) 
  810 FORMAT(' CI COMPONENTS FOR STATE',I5)
  820 FORMAT(I5,5X,F20.10)

      CONTAINS

      SUBROUTINE write_multistate_bin
       !---------------------------------------------------------------------
       ! Write unformatted (binary) file for data transfer to ChemShell
       ! Adapted for mutiple states from the original write_bin in SCFSAV
       ! Data written corresponds to formatted file obtained with NSAV=4
       ! Format of the file:
       ! I:    Length of content vector (= number of sections in file)
       ! I(n): Vector containing number of items in each section
       ! R:    n sections of length specified in content vector
       !       All items are reals
       !       (integer info such as atomic numbers and electronic states
       !        are already known by ChemShell)
       !---------------------------------------------------------------------
        IMPLICIT NONE        
        INTEGER(4),    PARAMETER      :: Nsec = 9 ! # of sections in file
        CHARACTER(32), PARAMETER      :: fname = 'fort.15'
        INTEGER(4)                    :: items(Nsec) = 0
        REAL(8),       ALLOCATABLE    :: MMdata(:,:)
        REAL(8),       ALLOCATABLE    :: cicmpx(:,:)
        INTEGER(4)                    :: numIcs = 0 ! no. of interstate coupling pairs
        INTEGER(4)                    :: ivar
        !--------------------------------------------------------------------
        OPEN(UNIT=NB15, FILE=TRIM(fname), STATUS='REPLACE',
     $        FORM='UNFORMATTED', ACTION='WRITE')
      
        ! Construct and write content vector
        items(1) = 3*NUMAT                  ! Coords of QM atoms
        IF (MMINP.GT.0) items(2) = 4*NUMATM ! Coords and charge of MM atoms
        items(3) = 3*NCIGRD                 ! Energies, cart. and internal grad norms
        IF (ICALL3.GE.1) THEN               ! Cart. gradients
           IF (MMINP.LE.0) THEN
              items(4) = 3*NUMAT*NCIGRD
           ELSE
              items(4) = 3*(NUMAT + NUMATM)*NCIGRD
           END IF
        END IF
        IF (ICALL4.GE.1) items(5) = NVAR*(1+NCIGRD)  ! Internal coords and gradients

        ! Note items 6, 7 and 8 here correspond to interstate coupling gradients
        ! Cf the original write_bin where 6 and 7 are electrostatic pot. and electric field
        numIcs = (NCIGRD*(NCIGRD-1))/2
        IF (ICALL3.EQ.2 .OR. ICALL4.EQ.2) items(6) = 2*numIcs ! gradient norms
        IF (ICALL3.EQ.2) THEN
           IF (MMINP.LE.0) THEN
              items(7) = 3*NUMAT*numIcs
           ELSE
              items(7) = 3*(NUMAT + NUMATM)*numIcs
           END IF
        END IF
        IF (ICALL4.EQ.2) items(8) = NVAR*numIcs
        IF (NSAV15.EQ.7) items(9) = NCICONF*NCIGRD

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

        ! SAVE ENERGIES AND GRADIENT NORMS (ALWAYS DEFINED)
        WRITE(NB15) ENERGX(1:NCIGRD), CNORMX(1:NCIGRD), GNORMX(1:NCIGRD)

        ! SAVE CARTESIAN GRADIENTS (IF DEFINED)
        ! CG is the full multiple state grad matrix, first QM, then MM
        IF (ICALL3.GE.1) WRITE(NB15) CGX(:,1:NUMAT+NUMATM,1:NCIGRD)

        ! SAVE INTERNAL COORDINATES AND GRADIENTS (IF DEFINED).
        IF (ICALL4.GE.1) THEN 
           WRITE(NB15) X(1:NVAR), GRADX(1:NVAR,1:NCIGRD)
        ENDIF

        ! SAVE INTERSTATE COUPLING NORMS (IF DEFINED).
        IF (ICALL3.EQ.2) WRITE(NB15) CCNRMX(1:numIcs), GCNRMX(1:numIcs)
 
        ! SAVE CARTESIAN INTERSTATE COUPLING GRADIENTS (IF DEFINED).
        IF (ICALL3.EQ.2) WRITE(NB15) CCX(:,1:NUMAT+NUMATM,1:numIcs)

        ! SAVE INTERNAL INTERSTATE COUPLING GRADIENTS (IF DEFINED).
        IF (ICALL4.EQ.2) WRITE(NB15) GCX(1:NVAR,1:numIcs)

        ! SAVE CI COMPONENTS IF REQUESTED
        ! CICOMP contains all NROOT components, so we have to do a bit
        ! of shuffling to get the needed NCIGRD components into a 
        ! single record.
        IF (NSAV15.EQ.7) THEN
           ALLOCATE(cicmpx(NCICONF,NCIGRD))
           DO ivar = 1, NCIGRD
              cicmpx(1:NCICONF,ivar) = cicomp(1:NCICONF,IGRST(ivar))
           END DO
           WRITE(NB15) cicmpx
           DEALLOCATE(cicmpx)
        ENDIF

        CLOSE(nb15)
        RETURN
      END SUBROUTINE write_multistate_bin
      
      END SUBROUTINE MULSAV
