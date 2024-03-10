C     ******************************************************************
C
C     Driver for two-center integrals and their derivatives (final)
C     in the molecular coordinate system.
C
C     ******************************************************************
      SUBROUTINE PSINTS(IMOD,IA,IB,W,W1,W2)
C
C   Compute two-center two-electron integrals,
C   electron-core repulsion integrals and core-core
C   repulsion integral along with first and second 
C   derivatives in molecular coordinate system
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IMOD   - Computation type
C           0 = Compute integrals
C           1 = Compute rotation matrix and first derivatives
C           2 = Compute rotation matrix with 1st and 2nd 
C               derivatives
C      IA     - (If positive) Number of first atom
C               (If negative) Number of first point charge with
C                             sign reversed.
C      IB     - (If positive) Number of second atom
C               (If negative) Number if second point charge with
C                             sign reversed.
C      W      - Output integrals
C      W1     - Output integrals first derivatives
C               (not accessed if IMOD.LE.0)
C      W2     - Output integrals second derivatives
C               (not accessed if IMOD.LE.1)
C
C   Accessed common blocks:
C
C      ATOMS  - Indices of first/last orbital on atom I
C      ATOMC  - Cartesian coordinates of atom I
C      QMM1   - Cartesian coordinates of point charges
C      PSROT0,
C      PSROT1,
C      PSROT2 - used to pass information between sibling
C               routines and should be defined here to ensure
C               they are not discarded by too smart compiler
C
C      PSDIST - used to pass interatomic distance to other
C               interested modules.
C
C   Local storage:
C
C      26700 DOUBLE PRECISION cells.
C      This could easily be decreased by 5.5K cells by
C      using output arrays for WLx and WLxP temporaries,
C      but I do not care _that_ much about memory use...
C
C   Module logic:
C
C      See separate documentation for call tree, COMMON
C      blocks, detailed timings, and possible speedups.
C  
C   Bugs:
C
C     Some common blocks initialized during PSINTS
C     execution are used elsewhere to save on
C     operations count. These are:
C
C         PSROT0, PSROT1, PSROT2, PSDROT and PSDIST
C
C     They should be explicitly SAVE'd in all modules
C     using them to be conformant to F77 standard
C
C     RMPD rotation matrix (and it's derivatives) are
C     computed by PSRM in more cases then strictly
C     necessary to rotate electron integrals - it is 
C     also used by PSOVER to rotate overlap integrals
C
      USE LIMIT, ONLY: LM1, LM1M
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (BOHR=0.529167D0)

      COMMON 
     ./PSDIST/ DIST
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./ATOMC / COORD(3,LM1)
     ./QMMM1 / COORDM(3,LM1M),CHARGM(LM1M)
      SAVE /PSDIST/
C
      DIMENSION W(46,46), W1(46,46,3), W2(46,46,6)
C
      DIMENSION WL(461), WL1(461), WL2(461)
      DIMENSION WL1P(461,3), WL2P(461,6)
      DIMENSION WRT(46,46), WRT1(46,46,3), WRT2(46,46,6)
      DIMENSION R(3,2)
C
      IF( IA.GT.0 ) THEN
          IORBA  = NLAST(IA) - NFIRST(IA)
          ITYPA  = NAT(IA)
          R(1,1) = COORD(1,IA)
          R(2,1) = COORD(2,IA)
          R(3,1) = COORD(3,IA)
      ELSE
          IORBA  = -1
          ITYPA  = 0
          R(1,1) = COORDM(1,-IA)
          R(2,1) = COORDM(2,-IA)
          R(3,1) = COORDM(3,-IA)
      ENDIF
      IF( IB.GT.0 ) THEN
          IORBB  = NLAST(IB) - NFIRST(IB)
          ITYPB  = NAT(IB)
          R(1,2) = COORD(1,IB)
          R(2,2) = COORD(2,IB)
          R(3,2) = COORD(3,IB)
      ELSE
          IORBB  = -1
          ITYPB  = 0
          R(1,2) = COORDM(1,-IB)
          R(2,2) = COORDM(2,-IB)
          R(3,2) = COORDM(3,-IB)
      ENDIF
      DIST  = SQRT( (R(1,1)-R(1,2))**2 + (R(2,1)-R(2,2))**2 +
     .              (R(3,1)-R(3,2))**2 ) / BOHR
C
      CALL PSRM  (1,2,IORBA,IORBB,R,IMOD)
      CALL PSIL  (ITYPA,ITYPB,IORBA,IORBB,DIST,IMOD,WL,WL1,WL2)
      CALL PSIDR (IMOD,IORBA,IORBB,WL1,WL1P,WL2,WL2P)
      CALL PSIR1 (IMOD,IORBA,IORBB,WL,WL1P,WL2P,WRT,WRT1,WRT2)
      CALL PSIR2 (IMOD,IORBA,IORBB,WRT,WRT1,WRT2,W,W1,W2)
      CALL PSISTF(IMOD,IORBA,IORBB,W,W1,W2)
      RETURN
      END
C
