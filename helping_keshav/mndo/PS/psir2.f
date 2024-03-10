C     ******************************************************************
C
C     Integral transformation: Rotate orbitals on the second atom.
C
C     ******************************************************************
      SUBROUTINE PSIR2(IMOD,IORBA,IORBB,AINT,AINT1,AINT2,OUT,OUT1,OUT2)
C
C   Second atom integrals rotation, with explicit use of zero
C   columns
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IMOD   - Computation type
C               0 = compute rotated integrals
C               1 = .... and first derivatives
C               2 = ......... and second derivatives
C      IORBA  - Number of orbitals on first atom
C               (0-based)
C      IORBB  - Number of orbitals on second atom
C               (0-based)
C      AINT   - Partially transformed integrals
C      AINT1  - Partially transformed first derivatives
C               (Not accessed if IMOD.LE.0)
C               First derivatives are in order X,Y,Z
C      AINT2  - Partially transformed second derivatives
C               (Not accessed if IMOD.LE.1)
C               Second derivatives are in order
C               XX, YY, ZZ, XY, XZ, YZ
C      OUT    - Transformed integrals
C      OUT1   - Transformed derivatives
C               (Not accessed if IMOD.LE.0)
C      OUT2   - Transformed second derivatives
C               (Not accessed if IMOD.LE.1)
C
C   Accessed common blocks:
C
C      PSROT0
C      PSROT1 - For IMOD.GE.1
C      PSROT2 - For IMOD.GE.2
C  
C   Modified common blocks:
C
C      None.
C
C   Local storage:
C
C      2120 DOUBLE PRECISION cells.
C
C      This could have been completely eliminated by
C      passing extra argument from PSINTS, which have
C      available temporary storage at the time of call.
C
C   Module logic:
C
C      PSIR2 calls sibling routine PSIR2X to perform
C      elementary transformation and combines its output
C      using ordinary differentiation rules (see PSGVML
C      for comments). We cannot perform rotation in-place
C      (although we could do it in principle for elementary
C      rotation), becase we need partially transformed
C      integrals and derivatives to deal with second
C      derivatives. We could have made TEMP array
C      unnecessary by processing integrals only when all
C      derivatives are already processed, though.
C
C      The following routines are called to perform the actual
C      computations:
C
C          PSIR2C - transforms SP and SD distributions
C          PSIR2X - transforms PP, PD and DD distributions
C                   using two-stage formulas
C          PSIR2Y - transforms PP, PD and DD distributions
C                   using direct formulas
C
C      Different routines are used to process integrals and
C      derivatives. For integrals, the two-stage transformation
C      of PP, PD and DD has the best operation count, while
C      for first and second derivatives the direct transformation
C      is better. Direct transformation routine (PSIR2Y) could
C      be used to transform integrals, too, albeit less
C      efficiently. Two-step transformation routine (PSIR2X)
C      should never be used to transform derivatives, because
C      it does not implement chain rules for two-stage
C      rotations.
C
C      Rotation matrix used (implicitly) for rotation of
C      SS and Core distributions is quite simple (single 1,
C      with obviously zero derivatives), so that SS distibutions
C      are processed (i.e., just copied) separately.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON 
     ./PSROT0/ RM (3,3), RMD (5,5), RMPP(6,6), RMPD(15,15), RMDD(15,15)
     ./PSROT1/ RM1(3,3,3), RMD1(5,5,3), RMPP1(6,6,3), RMPD1(15,15,3), 
     .         RMDD1(15,15,3)
     ./PSROT2/ RM2(3,3,6), RMD2(5,5,6), RMPP2(6,6,6), RMPD2(15,15,6), 
     .         RMDD2(15,15,6)
      SAVE /PSROT0/, /PSROT1/, /PSROT2/

      DIMENSION AINT(46,46), AINT1(46,46,3), AINT2(46,46,6)
      DIMENSION OUT (46,46), OUT1 (46,46,3), OUT2 (46,46,6)

      DIMENSION TEMP(46,46)
C
      IF( IORBA.LE.0 ) THEN
          IPAIRA =  2
      ELSE IF( IORBA.LE.3 ) THEN
          IPAIRA = 11
      ELSE 
          IPAIRA = 46
      ENDIF
C
C   copy Core and SS integrals
C
      DO 1020 IB = 1, 2
          DO 1010 IA = 1, IPAIRA
              OUT(IA,IB) = AINT(IA,IB)
 1010     CONTINUE
 1020 CONTINUE
      IF( IMOD .GE. 1 ) THEN
          DO 1050 JB = 1, 3
              DO 1040 IB = 1, 2
                  DO 1030 IA = 1, IPAIRA
                      OUT1(IA,IB,JB) = AINT1(IA,IB,JB)
 1030             CONTINUE
 1040         CONTINUE
 1050     CONTINUE
          IF( IMOD .GE. 2 ) THEN
              DO 1080 JB = 1, 6
                  DO 1070 IB = 1, 2
                      DO 1060 IA = 1, IPAIRA
                          OUT2(IA,IB,JB) = AINT2(IA,IB,JB)
 1060                 CONTINUE
 1070             CONTINUE
 1080         CONTINUE
          ENDIF
      ENDIF
C
      IF( IORBB.LE.0 ) RETURN
C
C   Integrals themselves
C
      CALL PSIR2C(IPAIRA,IORBB,AINT,RM,RMD,OUT)
      CALL PSIR2X(IPAIRA,IORBB,AINT,RM,RMD,OUT)
      IF( IMOD.LE.0 ) RETURN
C
      IF( IORBB.LE.0 ) THEN
          IPAIRB =  2
      ELSE IF( IORBB.LE.3 ) THEN
          IPAIRB = 11
      ELSE 
          IPAIRB = 46
      ENDIF
C
C   First derivatives
C
      DO 10 I=1,3
          CALL PSIR2C(IPAIRA,IORBB,AINT,RM1(1,1,I),RMD1(1,1,I),
     .                OUT1(1,1,I))
          CALL PSIR2Y(IPAIRA,IORBB,AINT,RMPP1(1,1,I),RMPD1(1,1,I),
     .                RMDD1(1,1,I),OUT1(1,1,I))
          CALL PSIR2C(IPAIRA,IORBB,AINT1(1,1,I),RM,RMD,TEMP)
          CALL PSIR2Y(IPAIRA,IORBB,AINT1(1,1,I),RMPP,RMPD,RMDD,TEMP)
          DO 296 JB=3,IPAIRB
              DO 295 JA=1,IPAIRA
                  OUT1(JA,JB,I) = OUT1(JA,JB,I) + TEMP(JA,JB)
  295         CONTINUE
  296     CONTINUE
   10 CONTINUE
C
      IF( IMOD.LE.1 ) RETURN
C
C   Second derivatives
C
C   Fxy*G and F*Gxy contributions
C
      DO 20 I=1,6
          CALL PSIR2C(IPAIRA,IORBB,AINT,RM2(1,1,I),RMD2(1,1,I),
     .                OUT2(1,1,I))
          CALL PSIR2Y(IPAIRA,IORBB,AINT,RMPP2(1,1,I),RMPD2(1,1,I),
     .                RMDD2(1,1,I),OUT2(1,1,I))
          CALL PSIR2C(IPAIRA,IORBB,AINT2(1,1,I),RM,RMD,TEMP)
          CALL PSIR2Y(IPAIRA,IORBB,AINT2(1,1,I),RMPP,RMPD,RMDD,TEMP)
          DO 396 JB=3,IPAIRB
              DO 395 JA=1,IPAIRA
                  OUT2(JA,JB,I) = OUT2(JA,JB,I) + TEMP(JA,JB)
  395         CONTINUE
  396     CONTINUE
   20 CONTINUE
C
C   Fx*Gx contributions - diagonal second derivatives
C
      DO 30 I=1,3
          CALL PSIR2C(IPAIRA,IORBB,AINT1(1,1,I),RM1(1,1,I),
     .                RMD1(1,1,I), TEMP)
          CALL PSIR2Y(IPAIRA,IORBB,AINT1(1,1,I),RMPP1(1,1,I),
     .                RMPD1(1,1,I), RMDD1(1,1,I), TEMP)
          DO 496 JB=3,IPAIRB
              DO 495 JA=1,IPAIRA
                  OUT2(JA,JB,I) = OUT2(JA,JB,I) + 2.0D0 * TEMP(JA,JB)
  495         CONTINUE
  496     CONTINUE
   30 CONTINUE
C
C   Fx*Gy and Fy*Gx contributions - off-diagonal second derivatives
C
C   I1,I2 pair successfully assumes values (1,2), (1,3), (2,3)
C
      DO 40 I=1,3
          I1 = MAX(I-1,1)
          I2 = MIN(I+1,3)
C
C   Fx*Gy contribution
C
          CALL PSIR2C(IPAIRA,IORBB,AINT1(1,1,I1),RM1(1,1,I2),
     .                RMD1(1,1,I2), TEMP)
          CALL PSIR2Y(IPAIRA,IORBB,AINT1(1,1,I1),RMPP1(1,1,I2),
     .                RMPD1(1,1,I2), RMDD1(1,1,I2), TEMP)
          DO 596 JB=3,IPAIRB
              DO 595 JA=1,IPAIRA
                  OUT2(JA,JB,I+3) = OUT2(JA,JB,I+3) + TEMP(JA,JB)
  595         CONTINUE
  596     CONTINUE
C
C   Fy*Gx contribution
C
          CALL PSIR2C(IPAIRA,IORBB,AINT1(1,1,I2),RM1(1,1,I1),
     .                RMD1(1,1,I1), TEMP)
          CALL PSIR2Y(IPAIRA,IORBB,AINT1(1,1,I2),RMPP1(1,1,I1),
     .                RMPD1(1,1,I1), RMDD1(1,1,I1), TEMP)
          DO 696 JB=3,IPAIRB
              DO 695 JA=1,IPAIRA
                  OUT2(JA,JB,I+3) = OUT2(JA,JB,I+3) + TEMP(JA,JB)
  695         CONTINUE
  696     CONTINUE
   40 CONTINUE
      RETURN
      END
C     ******************************************************************
      SUBROUTINE PSIR2C(IPAIRA,IORBB,AINT,RM,RMD,OUT)
C
C   Elementary second atom integrals rotation, distributions
C   SP and SD
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IPAIRA - Number of charge distributions on atom A
C               (1-based)
C      IORBB  - Number of orbitals on atom B
C               (0-based)
C      AINT   - Partially transformed integrals
C      RM     - P-orbitals "rotation" matrix
C      RMD    - D-orbitals "rotation" matrix
C      OUT    - Transformed integrals
C
C   Accessed common blocks:
C
C      None.
C  
C   Modified common blocks:
C
C      None.
C
C   Local storage:
C
C      None.
C
C   Module logic:
C
C      Although DAXPY semantics would fit perfectly for
C      this problem, I use loops instead, because
C      calling DAXPY will probably hurt more than it
C      could help for vectors that short.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION AINT(46,46), OUT(46,46)
      DIMENSION RM(3,3), RMD(5,5)
C
C   SP distributions
C
      DO 40 JB = 3, 5
          IB = JB - 2
          DO 30 IA = 1, IPAIRA
              OUT(IA,JB) = RM(1,IB)*AINT(IA,3) + 
     .                     RM(2,IB)*AINT(IA,4) + 
     .                     RM(3,IB)*AINT(IA,5)
   30     CONTINUE
   40 CONTINUE
C
      IF( IORBB.LE.3 ) RETURN
C
C   SD distributions
C
      DO 180 JB = 12, 16
          IB = JB - 11
          DO 170 IA = 1, IPAIRA
              OUT(IA,JB) = RMD(1,IB)*AINT(IA,12) + 
     .                     RMD(2,IB)*AINT(IA,13) + 
     .                     RMD(3,IB)*AINT(IA,14) +
     .                     RMD(4,IB)*AINT(IA,15) +
     .                     RMD(5,IB)*AINT(IA,16)
  170     CONTINUE
  180 CONTINUE
C
      RETURN
      END
C     ******************************************************************
      SUBROUTINE PSIR2X(IPAIRA,IORBB,AINT,RM,RMD,OUT)
C
C   Elementary second atom integrals rotation, with explicit 
C   use of zero columns, of PP, PD and DD distributions,
C   two-stage procedure.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IPAIRA - Number of charge distributions on atom A
C               (1-based)
C      IORBB  - Number of orbitals on atom B
C               (0-based)
C      AINT   - Partially transformed integrals
C      RM     - P-orbitals "rotation" matrix
C      RMD    - D-orbitals "rotation" matrix
C      OUT    - Transformed integrals
C
C   Accessed common blocks:
C
C      None.
C  
C   Modified common blocks:
C
C      None.
C
C   Local storage:
C
C      1150 DOUBLE PRECISION cells
C
C   Module logic:
C
C      Although DAXPY semantics would fit perfectly for
C      this problem, I use loops instead, because
C      calling DAXPY will probably hurt more than it
C      could help for vectors that short.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION AINT(46,46), OUT(46,46)
      DIMENSION RM(3,3), RMD(5,5)
      DIMENSION TEMP(46,25)
C
C   PP distributions, first stage
C
      DO 60 JB = 1, 3
          DO 50 IA = 1, IPAIRA
              TEMP(IA,JB) = RM(1,JB) * AINT(IA,6) +
     .                      RM(2,JB) * AINT(IA,7) +
     .                      RM(3,JB) * AINT(IA,9)
   50     CONTINUE
   60 CONTINUE
C
      DO 80 JB = 4, 6
          IB = JB - 3
          DO 70 IA = 1, IPAIRA
              TEMP(IA,JB) = RM(1,IB) * AINT(IA,7) +
     .                      RM(2,IB) * AINT(IA,8) +
     .                      RM(3,IB) * AINT(IA,10)
   70     CONTINUE
   80 CONTINUE
C
      DO 100 JB = 7, 9
          IB = JB - 6
          DO 90 IA = 1, IPAIRA
              TEMP(IA,JB) = RM(1,IB) * AINT(IA,9) +
     .                      RM(2,IB) * AINT(IA,10) +
     .                      RM(3,IB) * AINT(IA,11)
   90     CONTINUE
  100 CONTINUE
C
C   PP distributions, second stage
C
      DO 110 IA = 1, IPAIRA
          OUT(IA, 6) = RM(1,1) * TEMP(IA,1) +
     .                 RM(2,1) * TEMP(IA,4) +
     .                 RM(3,1) * TEMP(IA,7)  
  110 CONTINUE
C
      DO 120 IA = 1, IPAIRA
          OUT(IA, 7) = RM(1,2) * TEMP(IA,1) +
     .                 RM(2,2) * TEMP(IA,4) +
     .                 RM(3,2) * TEMP(IA,7)  
  120 CONTINUE
C
      DO 130 IA = 1, IPAIRA
          OUT(IA, 9) = RM(1,3) * TEMP(IA,1) +
     .                 RM(2,3) * TEMP(IA,4) +
     .                 RM(3,3) * TEMP(IA,7)  
  130 CONTINUE
C
      DO 140 IA = 1, IPAIRA
          OUT(IA, 8) = RM(1,2) * TEMP(IA,2) +
     .                 RM(2,2) * TEMP(IA,5) +
     .                 RM(3,2) * TEMP(IA,8)  
  140 CONTINUE
C
      DO 150 IA = 1, IPAIRA
          OUT(IA,10) = RM(1,3) * TEMP(IA,2) +
     .                 RM(2,3) * TEMP(IA,5) +
     .                 RM(3,3) * TEMP(IA,8)  
  150 CONTINUE
C
      DO 160 IA = 1, IPAIRA
          OUT(IA,11) = RM(1,3) * TEMP(IA,3) +
     .                 RM(2,3) * TEMP(IA,6) +
     .                 RM(3,3) * TEMP(IA,9)  
  160 CONTINUE
C
      IF( IORBB.LE.3 ) RETURN
C
C   PD distributions, first stage
C
      DO 200 JB = 1, 3
          DO 190 IA = 1, IPAIRA
              TEMP(IA,JB) = RM(1,JB)*AINT(IA,17) +
     .                      RM(2,JB)*AINT(IA,18)
  190     CONTINUE
  200 CONTINUE
C
      DO 220 JB = 4, 6
          IB = JB - 3
          DO 210 IA = 1, IPAIRA
              TEMP(IA,JB) = RM(1,IB)*AINT(IA,20) +
     .                      RM(3,IB)*AINT(IA,22)  
  210     CONTINUE
  220 CONTINUE
C
      DO 240 JB = 7, 9
          IB = JB - 6
          DO 230 IA = 1, IPAIRA
              TEMP(IA,JB) = RM(1,IB)*AINT(IA,23) +
     .                      RM(2,IB)*AINT(IA,24) +
     .                      RM(3,IB)*AINT(IA,25)  
  230     CONTINUE
  240 CONTINUE
C
      DO 260 JB = 10, 12
          IB = JB - 9
          DO 250 IA = 1, IPAIRA
              TEMP(IA,JB) = RM(2,IB)*AINT(IA,27) +
     .                      RM(3,IB)*AINT(IA,28)  
  250     CONTINUE
  260 CONTINUE
C
      DO 280 JB = 13, 15
          IB = JB - 12
          DO 270 IA = 1, IPAIRA
              TEMP(IA,JB) = RM(1,IB)*AINT(IA,29) +
     .                      RM(2,IB)*AINT(IA,30)  
  270     CONTINUE
  280 CONTINUE
C
C   PD distributions, second stage
C
      DO 310 IB = 1, 5
          DO 300 JB = 16, 18
              KB = JB+3*IB-2
              DO 290 IA = 1, IPAIRA
                  OUT(IA,KB) = RMD(1,IB)*TEMP(IA,JB-15) +
     .                         RMD(2,IB)*TEMP(IA,JB-12) +
     .                         RMD(3,IB)*TEMP(IA,JB- 9) +
     .                         RMD(4,IB)*TEMP(IA,JB- 6) +
     .                         RMD(5,IB)*TEMP(IA,JB- 3)
  290         CONTINUE
  300     CONTINUE
  310 CONTINUE
C
C   DD distributions, first stage
C
      DO 330 JB = 1, 5
          DO 320 IA = 1, IPAIRA
              TEMP(IA,JB) = RMD(1,JB)*AINT(IA,32) +
     .                      RMD(2,JB)*AINT(IA,33) +
     .                      RMD(3,JB)*AINT(IA,35) +
     .                      RMD(4,JB)*AINT(IA,38)
  320     CONTINUE
  330 CONTINUE
C
      DO 350 JB = 6, 10
          IB = JB - 5
          DO 340 IA = 1, IPAIRA
              TEMP(IA,JB) = RMD(1,IB)*AINT(IA,33) +
     .                      RMD(2,IB)*AINT(IA,34) +
     .                      RMD(3,IB)*AINT(IA,36) +
     .                      RMD(4,IB)*AINT(IA,39) +
     .                      RMD(5,IB)*AINT(IA,43)
  340     CONTINUE
  350 CONTINUE
C
      DO 370 JB = 11, 15
          IB = JB - 10
          DO 360 IA = 1, IPAIRA
              TEMP(IA,JB) = RMD(1,IB)*AINT(IA,35) +
     .                      RMD(2,IB)*AINT(IA,36) +
     .                      RMD(3,IB)*AINT(IA,37) +
     .                      RMD(4,IB)*AINT(IA,40) +
     .                      RMD(5,IB)*AINT(IA,44)
  360     CONTINUE
  370 CONTINUE
C
      DO 390 JB = 16, 20
          IB = JB - 15
          DO 380 IA = 1, IPAIRA
              TEMP(IA,JB) = RMD(1,IB)*AINT(IA,38) +
     .                      RMD(2,IB)*AINT(IA,39) +
     .                      RMD(3,IB)*AINT(IA,40) +
     .                      RMD(4,IB)*AINT(IA,41) +
     .                      RMD(5,IB)*AINT(IA,45)
  380     CONTINUE
  390 CONTINUE
C
      DO 410 JB = 21, 25
          IB = JB - 20
          DO 400 IA = 1, IPAIRA
              TEMP(IA,JB) = RMD(2,IB)*AINT(IA,43) +
     .                      RMD(3,IB)*AINT(IA,44) +
     .                      RMD(4,IB)*AINT(IA,45) +
     .                      RMD(5,IB)*AINT(IA,46)
  400     CONTINUE
  410 CONTINUE
C
C   DD distributions, second stage
C
      DO 420 IA = 1, IPAIRA
          OUT(IA,32) = RMD(1,1)*TEMP(IA, 1) +
     .                 RMD(2,1)*TEMP(IA, 6) +
     .                 RMD(3,1)*TEMP(IA,11) +
     .                 RMD(4,1)*TEMP(IA,16) +
     .                 RMD(5,1)*TEMP(IA,21)
  420 CONTINUE
C
      DO 430 IA = 1, IPAIRA
          OUT(IA,33) = RMD(1,2)*TEMP(IA, 1) +
     .                 RMD(2,2)*TEMP(IA, 6) +
     .                 RMD(3,2)*TEMP(IA,11) +
     .                 RMD(4,2)*TEMP(IA,16) +
     .                 RMD(5,2)*TEMP(IA,21)
  430 CONTINUE
C
      DO 440 IA = 1, IPAIRA
          OUT(IA,35) = RMD(1,3)*TEMP(IA, 1) +
     .                 RMD(2,3)*TEMP(IA, 6) +
     .                 RMD(3,3)*TEMP(IA,11) +
     .                 RMD(4,3)*TEMP(IA,16) +
     .                 RMD(5,3)*TEMP(IA,21)
  440 CONTINUE
C
      DO 450 IA = 1, IPAIRA
          OUT(IA,38) = RMD(1,4)*TEMP(IA, 1) +
     .                 RMD(2,4)*TEMP(IA, 6) +
     .                 RMD(3,4)*TEMP(IA,11) +
     .                 RMD(4,4)*TEMP(IA,16) +
     .                 RMD(5,4)*TEMP(IA,21)
  450 CONTINUE
C
      DO 460 IA = 1, IPAIRA
          OUT(IA,42) = RMD(1,5)*TEMP(IA, 1) +
     .                 RMD(2,5)*TEMP(IA, 6) +
     .                 RMD(3,5)*TEMP(IA,11) +
     .                 RMD(4,5)*TEMP(IA,16) +
     .                 RMD(5,5)*TEMP(IA,21)
  460 CONTINUE
C
      DO 470 IA = 1, IPAIRA
          OUT(IA,34) = RMD(1,2)*TEMP(IA, 2) +
     .                 RMD(2,2)*TEMP(IA, 7) +
     .                 RMD(3,2)*TEMP(IA,12) +
     .                 RMD(4,2)*TEMP(IA,17) +
     .                 RMD(5,2)*TEMP(IA,22)
  470 CONTINUE
C
      DO 480 IA = 1, IPAIRA
          OUT(IA,36) = RMD(1,3)*TEMP(IA, 2) +
     .                 RMD(2,3)*TEMP(IA, 7) +
     .                 RMD(3,3)*TEMP(IA,12) +
     .                 RMD(4,3)*TEMP(IA,17) +
     .                 RMD(5,3)*TEMP(IA,22)
  480 CONTINUE
C
      DO 490 IA = 1, IPAIRA
          OUT(IA,39) = RMD(1,4)*TEMP(IA, 2) +
     .                 RMD(2,4)*TEMP(IA, 7) +
     .                 RMD(3,4)*TEMP(IA,12) +
     .                 RMD(4,4)*TEMP(IA,17) +
     .                 RMD(5,4)*TEMP(IA,22)
  490 CONTINUE
C
      DO 500 IA = 1, IPAIRA
          OUT(IA,43) = RMD(1,5)*TEMP(IA, 2) +
     .                 RMD(2,5)*TEMP(IA, 7) +
     .                 RMD(3,5)*TEMP(IA,12) +
     .                 RMD(4,5)*TEMP(IA,17) +
     .                 RMD(5,5)*TEMP(IA,22)
  500 CONTINUE
C
      DO 510 IA = 1, IPAIRA
          OUT(IA,37) = RMD(1,3)*TEMP(IA, 3) +
     .                 RMD(2,3)*TEMP(IA, 8) +
     .                 RMD(3,3)*TEMP(IA,13) +
     .                 RMD(4,3)*TEMP(IA,18) +
     .                 RMD(5,3)*TEMP(IA,23)
  510 CONTINUE
C
      DO 520 IA = 1, IPAIRA
          OUT(IA,40) = RMD(1,4)*TEMP(IA, 3) +
     .                 RMD(2,4)*TEMP(IA, 8) +
     .                 RMD(3,4)*TEMP(IA,13) +
     .                 RMD(4,4)*TEMP(IA,18) +
     .                 RMD(5,4)*TEMP(IA,23)
  520 CONTINUE
C
      DO 530 IA = 1, IPAIRA
          OUT(IA,44) = RMD(1,5)*TEMP(IA, 3) +
     .                 RMD(2,5)*TEMP(IA, 8) +
     .                 RMD(3,5)*TEMP(IA,13) +
     .                 RMD(4,5)*TEMP(IA,18) +
     .                 RMD(5,5)*TEMP(IA,23)
  530 CONTINUE
C
      DO 540 IA = 1, IPAIRA
          OUT(IA,41) = RMD(1,4)*TEMP(IA, 4) +
     .                 RMD(2,4)*TEMP(IA, 9) +
     .                 RMD(3,4)*TEMP(IA,14) +
     .                 RMD(4,4)*TEMP(IA,19) +
     .                 RMD(5,4)*TEMP(IA,24)
  540 CONTINUE
C
      DO 550 IA = 1, IPAIRA
          OUT(IA,45) = RMD(1,5)*TEMP(IA, 4) +
     .                 RMD(2,5)*TEMP(IA, 9) +
     .                 RMD(3,5)*TEMP(IA,14) +
     .                 RMD(4,5)*TEMP(IA,19) +
     .                 RMD(5,5)*TEMP(IA,24)
  550 CONTINUE
C
      DO 560 IA = 1, IPAIRA
          OUT(IA,46) = RMD(1,5)*TEMP(IA, 5) +
     .                 RMD(2,5)*TEMP(IA,10) +
     .                 RMD(3,5)*TEMP(IA,15) +
     .                 RMD(4,5)*TEMP(IA,20) +
     .                 RMD(5,5)*TEMP(IA,25)
  560 CONTINUE
C
      RETURN
      END
C     ******************************************************************
      SUBROUTINE PSIR2Y(IPAIRA,IORBB,AINT,RMPP,RMPD,RMDD,OUT)
C
C   Elementary second atom integrals rotation, with explicit 
C   use of zero columns, of PP, PD and DD distributions,
C   direct procedure.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IPAIRA - Number of charge distributions on atom A
C               (1-based)
C      IORBB  - Number of orbitals on atom B
C               (0-based)
C      AINT   - Partially transformed integrals
C      RMPP   - PP-pairs "rotation" matrix
C      RMPD   - PD-pairs "rotation" matrix
C      RMDD   - DD-pairs "rotation" matrix
C      OUT    - Transformed integrals
C
C   Accessed common blocks:
C
C      None.
C  
C   Modified common blocks:
C
C      None.
C
C   Module logic:
C
C      Although DAXPY semantics would fit perfectly for
C      this problem, I use loops instead, because
C      calling DAXPY will probably hurt more than it
C      could help for vectors that short.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION AINT(46,46), OUT(46,46)
      DIMENSION RMPP(6,6), RMPD(15,15), RMDD(15,15)
C
C   PP distributions
C
      DO 40 JB = 6, 11
          IB = JB - 5
          DO 30 IA = 1, IPAIRA
              OUT(IA,JB) = RMPP(1,IB) * AINT(IA, 6) +
     .                     RMPP(2,IB) * AINT(IA, 7) +
     .                     RMPP(3,IB) * AINT(IA, 8) +
     .                     RMPP(4,IB) * AINT(IA, 9) +
     .                     RMPP(5,IB) * AINT(IA,10) +
     .                     RMPP(6,IB) * AINT(IA,11)
   30     CONTINUE
   40 CONTINUE
C
      IF( IORBB.LE.3 ) RETURN
C
C   PD distributions
C
      DO 60 JB = 17, 31
          IB = JB - 16
          DO 50 IA = 1, IPAIRA
              OUT(IA,JB) = RMPD( 1,IB) * AINT(IA, 17) +
     .                     RMPD( 2,IB) * AINT(IA, 18) +
     .                     RMPD( 4,IB) * AINT(IA, 20) +
     .                     RMPD( 6,IB) * AINT(IA, 22) +
     .                     RMPD( 7,IB) * AINT(IA, 23) +
     .                     RMPD( 8,IB) * AINT(IA, 24) +
     .                     RMPD( 9,IB) * AINT(IA, 25) +
     .                     RMPD(11,IB) * AINT(IA, 27) +
     .                     RMPD(12,IB) * AINT(IA, 28) +
     .                     RMPD(13,IB) * AINT(IA, 29) +
     .                     RMPD(14,IB) * AINT(IA, 30)
   50     CONTINUE
   60 CONTINUE
C
C   DD distributions
C
      DO 80 JB = 32, 46
          IB = JB - 31
          DO 70 IA = 1, IPAIRA
              OUT(IA,JB) = RMDD( 1,IB) * AINT(IA, 32) +
     .                     RMDD( 2,IB) * AINT(IA, 33) +
     .                     RMDD( 3,IB) * AINT(IA, 34) +
     .                     RMDD( 4,IB) * AINT(IA, 35) +
     .                     RMDD( 5,IB) * AINT(IA, 36) +
     .                     RMDD( 6,IB) * AINT(IA, 37) +
     .                     RMDD( 7,IB) * AINT(IA, 38) +
     .                     RMDD( 8,IB) * AINT(IA, 39) +
     .                     RMDD( 9,IB) * AINT(IA, 40) +
     .                     RMDD(10,IB) * AINT(IA, 41) +
     .                     RMDD(12,IB) * AINT(IA, 43) +
     .                     RMDD(13,IB) * AINT(IA, 44) +
     .                     RMDD(14,IB) * AINT(IA, 45) +
     .                     RMDD(15,IB) * AINT(IA, 46)
   70     CONTINUE
   80 CONTINUE
C
      RETURN
      END
