C     ******************************************************************
C
C     Driver for two-center integrals and their derivatives (local).
C
C     ******************************************************************
      SUBROUTINE PSIL(ITYPA,ITYPB,IORBA,IORBB,R,IMOD,W,W1,W2)
C
C   Compute two-center integrals and their derivatives in
C   local coordinate system. This is mostly dispatch
C   routine for PSILSS...PSILDD which do actial computations
C
C   Coded by: Serge Pachkovsky
C             Sibling routines coded by Mathematica program
C             coded by Serge Pachkovsky :-)
C
C   Parameters:
C
C      ITYPA  - Atom type of the first (left hand side) atom.
C               Zero is acceptable for point charge.
C      ITYPB  - Atom type of the second (right hand side) atom.
C               Zero is acceptable for point charge.
C      IORBA  - Number of orbitals on the first atom - 1.
C               -1 is acceptable for pure core.
C      IORBB  - Number of orbitals on the second atom - 1.
C               -1 is acceptable for pure core.
C      R      - Interatomic distance, atomic units 
C      IMOD   - Computation type
C           0 = Compute integrals
C           1 = Compute integrals and first derivatives
C           2 = Compute integrals with 1st and 2nd derivatives
C      W      - Output array, integrals
C      W1     - Output array, first derivatives of integrals
C               along A-B axis
C      W2     - Output array, second derivatives of integrals
C               along A-B axis
C
C   Only non-zero integrals are stored in W, W1 and W2.
C
C   Accessed common blocks:
C
C      MULTIP   Using both DD and PO.
C  
C   Modified common blocks:
C
C      None.
C
C   Module logic:
C
C      Actual computations are performed by one of 12 sibling
C      modules (generated automatically from the same set of
C      integrals definitions). Task of PSIL is to fetch necessary
C      parameters from common blocks and select one of PSILXX
C      modules, where XX stands for:
C           
C         SS - S orbitals on A, S orbitals on B
C         SX - S orbitals on A, S orbitals on B, same parameters
C         SP - S orbitals on A, P orbitals on B
C         SD - S orbitals on A, D orbitals on B
C         PS - P orbitals on A, S orbitals on B
C         PP - P orbitals on A, P orbitals on B
C         PP - P orbitals on A, P orbitals on B, same parameters
C         PD - P orbitals on A, D orbitals on B
C         DS - D orbitals on A, S orbitals on B
C         DP - D orbitals on A, P orbitals on B
C         DD - D orbitals on A, D orbitals on B
C         DD - D orbitals on A, D orbitals on B, same parameters
C
C      It also enforced MNDO convention RO0PP = ROCORE = ROSS.
C
C   Bugs:
C
C      IBM AIX XL FORTRAN Compiler/6000  Version 02.03.0000.0000
C      runs out of "system resources" (whatever it means) compiling
C      psildd.f with optimisation at -O2 level. Version 02.03.0000.0029 
C      works just fine (albeit sloooow).
C
      USE LIMIT, ONLY: LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (SQRT2=1.414213562373095048801688724209698078569672D0)
      PARAMETER (ZERO=0.D0)
      COMMON
     ./MULTIP/ DD(6,0:LMZ),PO(9,0:LMZ)
      DIMENSION W(461), W1(461), W2(461)
      DIMENSION ARO(9), BRO(9), AD(5), BD(5)

      ICASE  = 1

      ARO(1) = PO(9,ITYPA)
      ARO(2) = PO(1,ITYPA)

      IF( IORBA .GT. 0 ) THEN
          ICASE  = ICASE + 1
          AD(1)  = DD(2,ITYPA)
          AD(2)  = DD(3,ITYPA) * SQRT2
          ARO(3) = PO(2,ITYPA)
          ARO(4) = PO(3,ITYPA)
          ARO(5) = PO(7,ITYPA)
          IF( IORBA .GT. 3 ) THEN
              ICASE  = ICASE + 1
              AD(4)  = DD(4,ITYPA)
              AD(3)  = DD(5,ITYPA)
              AD(5)  = DD(6,ITYPA)
              ARO(8) = PO(4,ITYPA)
              ARO(7) = PO(5,ITYPA)
              ARO(9) = PO(6,ITYPA)
              ARO(6) = PO(8,ITYPA)
          ENDIF
          IF( ARO(5).EQ.ZERO ) ARO(5) = ARO(2)
      ENDIF
      IF( ARO(1).EQ.ZERO ) ARO(1) = ARO(2)
C
      IF( ITYPA .NE. ITYPB ) GOTO 190
C
C   Special case: atom A is of the same type as atom B,
C   use simplified expressions.
C
      GOTO (110,120,130), ICASE
110   CALL PSILSX(R,ARO,IMOD,W,W1,W2)
      RETURN
120   CALL PSILPX(R,ARO,AD,IMOD,W,W1,W2)
      RETURN
130   CALL PSILDX(R,ARO,AD,IMOD,W,W1,W2)
      RETURN

190   CONTINUE

      BRO(1) = PO(9,ITYPB)
      BRO(2) = PO(1,ITYPB)

      IF( IORBB .GT. 0 ) THEN
          ICASE  = ICASE + 3
          BD(1)  = DD(2,ITYPB)
          BD(2)  = DD(3,ITYPB) * SQRT2
          BRO(3) = PO(2,ITYPB)
          BRO(4) = PO(3,ITYPB)
          BRO(5) = PO(7,ITYPB)
          IF( IORBB .GT. 3 ) THEN
              ICASE  = ICASE + 3
              BD(4)  = DD(4,ITYPB)
              BD(3)  = DD(5,ITYPB)
              BD(5)  = DD(6,ITYPB)
              BRO(1) = PO(9,ITYPB)
              BRO(8) = PO(4,ITYPB)
              BRO(7) = PO(5,ITYPB)
              BRO(9) = PO(6,ITYPB)
              BRO(6) = PO(8,ITYPB)
          ENDIF
          IF( BRO(5).EQ.ZERO ) BRO(5) = BRO(2)
      ENDIF
      IF( BRO(1).EQ.ZERO ) BRO(1) = BRO(2)
C
C   This horrid GOTO should simulate table of pointers to
C   functions :-)
C
      GOTO (210, 220, 230, 240, 250, 260, 270, 280, 290), ICASE
210   CALL PSILSS(R,ARO,BRO,IMOD,W,W1,W2)
      RETURN
220   CALL PSILPS(R,ARO,BRO,AD,IMOD,W,W1,W2)
      RETURN
230   CALL PSILDS(R,ARO,BRO,AD,IMOD,W,W1,W2)
      RETURN
240   CALL PSILSP(R,ARO,BRO,BD,IMOD,W,W1,W2)
      RETURN
250   CALL PSILPP(R,ARO,BRO,AD,BD,IMOD,W,W1,W2)
      RETURN
260   CALL PSILDP(R,ARO,BRO,AD,BD,IMOD,W,W1,W2)
      RETURN
270   CALL PSILSD(R,ARO,BRO,BD,IMOD,W,W1,W2)
      RETURN
280   CALL PSILPD(R,ARO,BRO,AD,BD,IMOD,W,W1,W2)
      RETURN
290   CALL PSILDD(R,ARO,BRO,AD,BD,IMOD,W,W1,W2)
      RETURN
      END
