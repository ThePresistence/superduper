C     ******************************************************************
C
C     Reorder two-center integrals.
C
C     ******************************************************************
      SUBROUTINE PSISTF(IMOD,IORBA,IORBB,AINT,AINT1,AINT2)
C
C   Repack integrals in the format requered by external
C   world.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IMOD   - Type of computation:
C           0 = repack integrals
C           1 = repack integrals and first derivatives
C           2 = repack integrals, first and second derivatives
C      IORBA  - Number of orbitals on first atom (0-based)
C      IORBB  - Number of orbitals on second atom (0-based)
C      AINT   - Input integrals (internal order)
C               Output integrals (external order)
C      AINT1  - Input first derivatives (internal order)
C               Output first derivatives (external order)
C               (Not accessed if IMOD <= 0)
C      AINT2  - Input second derivatives (internal order)
C               Output second derivatives (external order)
C               (Not accessed if IMOD <= 1)
C
C   Accessed common blocks:
C
C      PSISTC - control string for permutor, accessed by
C               PSISTS, PSISTP and PSISTD. Defined by
C               block data PSISTB
C  
C   Modified common blocks:
C
C      None.
C
C   Local storage:
C
C      None
C
C   Module logic:
C
C      Decomposition of permutation of integrals to the set of
C      smallest cycles is used to transform integrals to a new order.
C      Therefore, only single row of scratch space is necessary.
C      Actual transformation is performed by routines PSISTS, PSISTP 
C      and PSISTD, which are computer-generated. Permutation functions
C      are driven by control string in IPERM array, which have the
C      following format:
C
C          Number of elements in cycle1 - 1 or -1
C          First element of cycle1
C          ....
C          Last element of cycle1
C          Number of elements in cycle2 - 1
C          ....
C          If (-1) is encounteres in count field,
C          it is followed by termination condition on IORBB
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AINT(46,46), AINT1(46,46,3), AINT2(46,46,6)
C
C  Shortcut: if both A and B have only S-type orbitals, they are
C  already in order.
C
      IF( IORBA.LE.0 .AND. IORBB.LE.0 ) RETURN
      IF( IORBA.LE.0 ) THEN
          CALL PSISTS( IORBB, AINT )
          IF( IMOD.LE.0 ) RETURN
          DO 110 I=1,3
              CALL PSISTS( IORBB, AINT1(1,1,I) )
  110     CONTINUE
          IF( IMOD.LE.1 ) RETURN
          DO 120 I=1,6
              CALL PSISTS( IORBB, AINT2(1,1,I) )
  120     CONTINUE
      ELSE IF( IORBA.LE.3 ) THEN
          CALL PSISTP( IORBB, AINT )
          IF( IMOD.LE.0 ) RETURN
          DO 210 I=1,3
              CALL PSISTP( IORBB, AINT1(1,1,I) ) 
  210     CONTINUE
          IF( IMOD.LE.1 ) RETURN
          DO 220 I=1,6
              CALL PSISTP( IORBB, AINT2(1,1,I) )
  220     CONTINUE
      ELSE 
          CALL PSISTD( IORBB, AINT )
          IF( IMOD.LE.0 ) RETURN
          DO 310 I=1,3
              CALL PSISTD( IORBB, AINT1(1,1,I) )
  310     CONTINUE
          IF( IMOD.LE.1 ) RETURN
          DO 320 I=1,6
              CALL PSISTD( IORBB, AINT2(1,1,I) )
  320     CONTINUE
      ENDIF
      RETURN
      END
