C     ******************************************************************
C
C     Integral transformation: Rotate orbitals on first atom.
C
C     ******************************************************************
      SUBROUTINE PSIR1(IMOD,IORBA,IORBB,AINT,AINT1,AINT2,OUT,OUT1,OUT2)
C
C   First atom integrals rotation with explicit 
C   use of zero elements.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IMOD   - Computation type
C               0 = compute generalized product
C               1 = .... and first derivatives
C               2 = ......... and second derivatives
C      IORBA  - Number of orbitals on first atom (0-based)
C               -1 is acceptable for pure cores.
C      IORBB  - Number of orbitals on second atom (0-based)
C               -1 is acceptable for pure cores.
C      AINT   - Packed array of integrals in
C               local coordinate system
C      AINT1  - Packed array of first derivatives
C               (Not accessed if IMOD.LE.0)
C               First derivatives are in order X,Y,Z
C      AINT2  - Packed array of second derivatives
C               (Not accessed if IMOD.LE.1)
C               Second derivatives are in order
C               XX, YY, ZZ, XY, XZ, YZ
C      OUT    - Partially transformed integrals
C      OUT1   - Partially transformed derivatives
C               (Not accessed if IMOD.LE.0)
C      OUT2   - Partially transformed second derivatives
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
C      This could be eliminated completely by passing extra
C      argument from PSINTS (which do have extra storage 
C      at the time PSIR1 works), but it will make the code
C      less modular...
C
C   Module logic:
C
C      This is special case of PSGVML, see it for comments
C      Probably it could be worthwile to split this function
C      into three separate functions instead of doing indexed
C      GOTO's on each function call. Unfortunately, F90's
C      function pointers are of no use here, because our
C      routines have different footprints.
C
C      PSIR1 calls three sibling routines: PSIR1S, PSIR1P
C      and PSIR1D, which differ by number of orbitals
C      on the first atom. PSIR1P and PSIR1D do not attempt
C      to transform integrals involving Core and SS distributions
C      on atom A, because transformation coefficient (1.0D0)
C      implicitly used to transform them do not fit into
C      general scheme and is treated separately
C
C   Bugs:
C
C      Pure cores are treated as S-type atoms.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON 
     ./PSROT0/ RM (3,3), RMD (5,5), RMPP(6,6), RMPD(15,15), RMDD(15,15)
     ./PSROT1/ RM1(3,3,3), RMD1(5,5,3), RMPP1(6,6,3), RMPD1(15,15,3), 
     .         RMDD1(15,15,3)
     ./PSROT2/ RM2(3,3,6), RMD2(5,5,6), RMPP2(6,6,6), RMPD2(15,15,6), 
     .         RMDD2(15,15,6)
      SAVE /PSROT0/, /PSROT1/, /PSROT2/

      DIMENSION AINT(461), AINT1(461,3), AINT2(461,6)
      DIMENSION OUT(46,46), OUT1(46,46,3), OUT2(46,46,6)

      DIMENSION TEMP(46,46)
C
C   Special treatment of Core and SS integrals
C
      CALL PSIR1S(IORBB,AINT,OUT)
      IF( IMOD.GE.1 ) THEN
          DO 1010 I=1,3
              CALL PSIR1S(IORBB,AINT1(1,I),OUT1(1,1,I))
 1010     CONTINUE
          IF( IMOD.GE.2 ) THEN
              DO 1020 I=1,6
                  CALL PSIR1S(IORBB,AINT2(1,I),OUT2(1,1,I))
 1020         CONTINUE
          ENDIF
      ENDIF
      IF( IORBA.LE.0 ) RETURN
C
C   General-case treatment of all orbitals with explicitly
C   computed transformation matrices
C
C
C   Integrals themselves
C
      IF( IORBA.LE.3 ) THEN
          CALL PSIR1P(IORBB,AINT,OUT,RM,RMPP)
      ELSE
          CALL PSIR1D(IORBB,AINT,OUT,RM,RMPP,RMD,RMPD,RMDD)
      ENDIF
      IF( IMOD.LE.0 ) RETURN
C
      IF( IORBA.EQ.3 ) IPAIRA = 11
      IF( IORBA.EQ.8 ) IPAIRA = 46
      IF( IORBB.LE.0 ) IPAIRB =  2
      IF( IORBB.EQ.3 ) IPAIRB = 11
      IF( IORBB.EQ.8 ) IPAIRB = 46
C
C   First derivatives
C
      DO 10 I=1,3
          IF( IORBA.LE.3 ) THEN
              CALL PSIR1P(IORBB,AINT,OUT1(1,1,I),RM1(1,1,I),
     .             RMPP1(1,1,I))
              CALL PSIR1P(IORBB,AINT1(1,I),TEMP,RM,RMPP)
          ELSE
              CALL PSIR1D(IORBB,AINT,OUT1(1,1,I),RM1(1,1,I),
     .             RMPP1(1,1,I), RMD1(1,1,I),RMPD1(1,1,I),RMDD1(1,1,I))
              CALL PSIR1D(IORBB,AINT1(1,I),TEMP,RM,RMPP,RMD,RMPD,RMDD)
          ENDIF
          DO 296 JB=1,IPAIRB
              DO 295 JA=3,IPAIRA
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
          IF( IORBA.LE.3 ) THEN
              CALL PSIR1P(IORBB,AINT,OUT2(1,1,I),RM2(1,1,I),
     .             RMPP2(1,1,I))
              CALL PSIR1P(IORBB,AINT2(1,I),TEMP,RM,RMPP)
          ELSE
              CALL PSIR1D(IORBB,AINT,OUT2(1,1,I),RM2(1,1,I),
     .             RMPP2(1,1,I), RMD2(1,1,I),RMPD2(1,1,I),RMDD2(1,1,I))
              CALL PSIR1D(IORBB,AINT2(1,I),TEMP,RM,RMPP,RMD,RMPD,RMDD)
          ENDIF
          DO 396 JB=1,IPAIRB
              DO 395 JA=3,IPAIRA
                  OUT2(JA,JB,I) = OUT2(JA,JB,I) + TEMP(JA,JB)
  395         CONTINUE
  396     CONTINUE
   20 CONTINUE
C
C   Fx*Gx contributions - diagonal second derivatives
C
      DO 30 I=1,3
          IF( IORBA.LE.3 ) THEN
              CALL PSIR1P(IORBB,AINT1(1,I),TEMP,RM1(1,1,I),
     .             RMPP1(1,1,I))
          ELSE
              CALL PSIR1D(IORBB,AINT1(1,I),TEMP,RM1(1,1,I),
     .             RMPP1(1,1,I), RMD1(1,1,I),RMPD1(1,1,I),RMDD1(1,1,I))
          ENDIF
          DO 496 JB=1,IPAIRB
              DO 495 JA=3,IPAIRA
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
          IF( IORBA.LE.3 ) THEN
              CALL PSIR1P(IORBB,AINT1(1,I1),TEMP,RM1(1,1,I2),
     .             RMPP1(1,1,I2))
          ELSE
              CALL PSIR1D(IORBB,AINT1(1,I1),TEMP,RM1(1,1,I2),
     .             RMPP1(1,1,I2),RMD1(1,1,I2),RMPD1(1,1,I2),
     .             RMDD1(1,1,I2))
          ENDIF
          DO 596 JB=1,IPAIRB
              DO 595 JA=3,IPAIRA
                  OUT2(JA,JB,I+3) = OUT2(JA,JB,I+3) + TEMP(JA,JB)
  595         CONTINUE
  596     CONTINUE
C
C   Fy*Gx contribution
C
          IF( IORBA.LE.3 ) THEN
              CALL PSIR1P(IORBB,AINT1(1,I2),TEMP,RM1(1,1,I1),
     .             RMPP1(1,1,I1))
          ELSE
              CALL PSIR1D(IORBB,AINT1(1,I2),TEMP,RM1(1,1,I1),
     .             RMPP1(1,1,I1),RMD1(1,1,I1),RMPD1(1,1,I1),
     .             RMDD1(1,1,I1))
          ENDIF
          DO 696 JB=1,IPAIRB
              DO 695 JA=3,IPAIRA
                  OUT2(JA,JB,I+3) = OUT2(JA,JB,I+3) + TEMP(JA,JB)
  695         CONTINUE
  696     CONTINUE
   40 CONTINUE
      RETURN
      END
