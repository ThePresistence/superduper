C     ******************************************************************
C
C     Compute rotation matrices and their derivatives.
C
C     ******************************************************************
      SUBROUTINE PSRM(IA,IB,IORBA,IORBB,COORD,IMODI)
C
C   Compute rotation matrices necessary for integrals
C   transformations.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IA     - Number of first atom
C      IB     - Number of second atom
C      IORBA  - Number of orbitals on the first atom - 1.
C               -1 is acceptable for pure cores, and is
C               treated the same way as 0.
C      IORBB  - Number of orbitals on the second atom - 1.
C               -1 is acceptable for pure cores, and is
C               treated the same way as 0.
C      COORD  - Coordinates of atoms
C      IMODI  - Computation type
C           0 = Compute rotation matrix RM
C           1 = Compute rotation matrix and first derivatives
C           2 = Compute rotation matrix with 1st and 2nd 
C               derivatives
C        +100 = (i.e., 100, 101, or 102) requests only rotation
C               matrices necessary for overlap-type integrals.
C               This is useful for OM1/OM2 methods, where ERI
C               integrals come from a different source. (Actually,
C               all rotation matrices are needed anyway, but it's
C               nice to feel smart).
C
C   Accessed common blocks:
C
C      None.
C  
C   Modified common blocks:
C
C      PSROT0 - Always modified by siblings. 
C      PSROT1 - Modified if IMOD .GE. 1
C      PSROT2 - Modified if IMOD .GE. 2
C
C          RM   is computed unless IMOD.EQ.0 and both atoms have
C               only S-type orbitals.
C          RMD  is computed if at least one of the atoms have 
C               D-type orbitals.
C          RMPP is computed if first atom have P-type orbitals
C               or second atom have P-type orbitals and derivatives
C               were requested.
C          RMPD is computed if first atom have D-type orbitals,
C               or second atom have D-type orbitals and derivatives
C               were requested, or first atom have P-type orbitals
C               and second atom have D-type orbitals.
C          RMDD is computed if first atom have D-type orbitals, or
C               second atom have D-type orbitals and derivatives 
C               were requested.
C               
C          First derivatives of rotation matrix are in X,Y,Z order
C          Second derivatives are in XX, YY, ZZ, XY, XZ, YZ order
C
C      PSDROT - Projection coefficients for gerivatives, see PSIDR
C               for definitions of coefficients and expressions.
C
C   Local storage:
C
C      None.
C
C   Module logic:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (BOHR=0.529167D0)
C
      COMMON
     ./PSROT0/ RM (3,3), RMD (5,5), RMPP(6,6), RMPD(15,15), RMDD(15,15)
     ./PSROT1/ RM1(3,3,3), RMD1(5,5,3), RMPP1(6,6,3), RMPD1(15,15,3), 
     .         RMDD1(15,15,3)
     ./PSDROT/ D1(3), D2A(6), D2B(6)
      SAVE /PSROT0/, /PSROT1/, /PSDROT/
C
      LOGICAL NEEDPP, NEEDPD, NEEDDD
      DIMENSION COORD(3,*)
C
      IMOD   = MOD(IMODI,100)
      IF( IORBA.LE.0 .AND. IORBB.LE.0 .AND. IMOD.EQ.0 ) RETURN
C
C     Compute interatomic distance (we need 'em in Bohrs)
C
      DX = (COORD(1,IA) - COORD(1,IB))/BOHR
      DY = (COORD(2,IA) - COORD(2,IB))/BOHR
      DZ = (COORD(3,IA) - COORD(3,IB))/BOHR
C    Compute basic rotation matrix and its derivatives
      CALL PSRMP(DX,DY,DZ,IMOD)
      IF( IORBA.GT.3 .OR. IORBB.GT.3 ) THEN
          CALL PSRMD(IMOD)
      ENDIF
C    Compute gradient projection coefficients
      IF( IMOD.GT.0 ) THEN
          D1(1) = -RM(3,1)
          D1(2) = -RM(3,2)
          D1(3) = -RM(3,3)
          IF( IMOD.GT.1 ) THEN
              D2A(1) = RM(3,1) ** 2
              D2A(2) = RM(3,2) ** 2
              D2A(3) = RM(3,3) ** 2
              D2A(4) = RM(3,1) * RM(3,2)
              D2A(5) = RM(3,1) * RM(3,3)
              D2A(6) = RM(3,2) * RM(3,3)
              D2B(1) = - RM1(3,1,1)
              D2B(2) = - RM1(3,2,2)
              D2B(3) = - RM1(3,3,3)
              D2B(4) = - RM1(3,1,2)
              D2B(5) = - RM1(3,1,3)
              D2B(6) = - RM1(3,2,3)
          ENDIF
      ENDIF
C    Compute necessary combination rotation matrices
      NEEDPP = .FALSE.
      NEEDPD = .FALSE.
      NEEDDD = .FALSE.
      IF( IORBA.GT.0 ) THEN
          NEEDPP = .TRUE.
          IF( IORBB.GT.3 ) NEEDPD = .TRUE.
      ENDIF
      IF( IORBA.GT.3 ) THEN
          NEEDPD = .TRUE.
          NEEDDD = .TRUE.
      ELSE
          IF( IMOD .GE.1 ) THEN
              IF( IORBB.GT.0 ) NEEDPP = .TRUE.
              IF( IORBB.GT.3 ) THEN
                  NEEDPD = .TRUE.
                  NEEDDD = .TRUE.
              ENDIF
          ENDIF
      ENDIF
      IF( NEEDPP ) CALL PSRMPP(IMOD)
      IF( NEEDPD ) CALL PSRMPD(IMOD)
      IF( NEEDDD ) CALL PSRMDD(IMOD)
      RETURN
      END
C
      SUBROUTINE PSRMP(DX,DY,DZ,IMOD)
C
C   Compute basic rotation matrix and it's derivatives
C
C   Coded by: Serge Pachkovsky, with his faithful
C             Mathematica :-)
C
C   Parameters:
C
C      DX, DY, DZ - coordinates of nucleus A as seen
C                   from the nucleus B in Bohrs
C      IMOD        - Computation type
C                0 = Compute rotation matrix RM
C                1 = Compute rotation matrix and first derivatives
C                2 = Compute rotation matrix with 1st and 2nd 
C                    derivatives
C
C   Accessed common blocks:
C
C      None.
C
C   Local storage:
C
C      None.
C  
C   Modified common blocks:
C
C      PSROT0, PSROT1, PSROT2 - rotation matrix RM and, depending
C          on IMOD, it's derivatives RM1 and RM2 are computed.
C          First derivatives of rotation matrix are in X,Y,Z order
C          Second derivatives are in XX, YY, ZZ, XY, XZ, YZ order
C
C   Module logic:
C
C      Standard rotation matrix used elsewhere in program have
C      unbounded derivatives then DY and DZ simultaneously approach
C      zero. Therefore, PSRM will use the "old" rotation matrix 
C      form then DY .GE. DX, and alternative one otherwise. As a
C      consequence, computed rotation matrix is discontinuous 
C      along plane X = Y. This have no bearing on the validity
C      of rotated integrals, but could cause confusion in attempts
C      to verify self-consistency of rotation matrix code by
C      evaluating it's derivatives numerically.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON 
     ./PSROT0/ RM (3,3), RMD (5,5), RMPP(6,6), RMPD(15,15), RMDD(15,15)
     ./PSROT1/ RM1(3,3,3), RMD1(5,5,3), RMPP1(6,6,3), RMPD1(15,15,3), 
     .         RMDD1(15,15,3)
     ./PSROT2/ RM2(3,3,6), RMD2(5,5,6), RMPP2(6,6,6), RMPD2(15,15,6), 
     .         RMDD2(15,15,6)
      SAVE /PSROT0/, /PSROT1/, /PSROT2/

      RN2  = 1.D0 / ( DX**2 + DY**2 + DZ**2 )
      RN   = SQRT( RN2 )
      X    = DX * RN
      Y    = DY * RN 
      Z    = DZ * RN
      X2   = X ** 2
      Y2   = Y ** 2
      Z2   = Z ** 2
      SXZ2 = X2 + Z2
      SYZ2 = Y2 + Z2
      SYZ  = SQRT( SYZ2 )
      SXZ  = SQRT( SXZ2 )
      IF( IMOD .GT. 0 ) THEN
          SXY2 = X2+Y2
          X4   = X2**2
          Y4   = Y2**2
          Z4   = Z2**2
      ENDIF

      IF( SYZ .GT. SXZ ) THEN
C
C         Case A (no zero Syz) matrix
C         Rotation matrix as a whole:
C
C
C         Syz      0        -x
C
C           x y     z
C         -(---)   ---
C           Syz    Syz      -y
C
C           x z       y
C         -(---)   -(---)
C           Syz      Syz    -z
C
          RM(1,1) = SYZ
          RM(2,1) = 0
          RM(3,1) = -X
          RM(1,2) = -(X*Y/SYZ)
          RM(2,2) = Z/SYZ
          RM(3,2) = -Y
          RM(1,3) = -(X*Z/SYZ)
          RM(2,3) = -(Y/SYZ)
          RM(3,3) = -Z
          IF( IMOD.LE.0 ) RETURN
C
          SYZ3 = SYZ * SYZ2
          RM1(1,1,1) = -(RN*SYZ*X)
          RM1(1,1,2) = RN*X2*Y/SYZ
          RM1(1,1,3) = RN*X2*Z/SYZ
          RM1(2,1,1) = 0
          RM1(2,1,2) = 0
          RM1(2,1,3) = 0
          RM1(3,1,1) = -(RN*SYZ2)
          RM1(3,1,2) = RN*X*Y
          RM1(3,1,3) = RN*X*Z
          RM1(1,2,1) = -(RN*SYZ*Y)
          RM1(1,2,2) = RN*X*(Y4 - Z2 + Y2*Z2)/SYZ3
          RM1(1,2,3) = RN*X*(2 - X2)*Y*Z/SYZ3
          RM1(2,2,1) = 0
          RM1(2,2,2) = -(RN*Y*Z/SYZ3)
          RM1(2,2,3) = RN*Y2/SYZ3
          RM1(3,2,1) = RN*X*Y
          RM1(3,2,2) = -(RN*SXZ2)
          RM1(3,2,3) = RN*Y*Z
          RM1(1,3,1) = -(RN*SYZ*Z)
          RM1(1,3,2) = RN*X*(2 - X2)*Y*Z/SYZ3
          RM1(1,3,3) = RN*X*(-Y2 + Y2*Z2 + Z4)/SYZ3
          RM1(2,3,1) = 0
          RM1(2,3,2) = -(RN*Z2/SYZ3)
          RM1(2,3,3) = RN*Y*Z/SYZ3
          RM1(3,3,1) = RN*X*Z
          RM1(3,3,2) = RN*Y*Z
          RM1(3,3,3) = -(RN*SXY2)
          IF( IMOD.LE.1 ) RETURN
C
          SYZ5 = SYZ3 * SYZ2
          RM2(1,1,1) = RN2*SYZ*(-1 + 3*X2)
          RM2(1,1,2) = RN2*X2*(1 - X2 - 4*Y2 + 3*X2*Y2)/SYZ3
          RM2(1,1,3) = RN2*X2*(1 - X2 - 4*Z2 + 3*X2*Z2)/SYZ3
          RM2(1,1,4) = RN2*X*(2 - 3*X2)*Y/SYZ
          RM2(1,1,5) = RN2*X*(2 - 3*X2)*Z/SYZ
          RM2(1,1,6) = RN2*X2*(-4 + 3*X2)*Y*Z/SYZ3
          RM2(2,1,1) = 0
          RM2(2,1,2) = 0
          RM2(2,1,3) = 0
          RM2(2,1,4) = 0
          RM2(2,1,5) = 0
          RM2(2,1,6) = 0
          RM2(3,1,1) = 3*RN2*SYZ2*X
          RM2(3,1,2) = RN2*X*(SXZ2 - 2*Y2)
          RM2(3,1,3) = RN2*X*(SXY2 - 2*Z2)
          RM2(3,1,4) = RN2*(SYZ2 - 2*X2)*Y
          RM2(3,1,5) = RN2*(SYZ2 - 2*X2)*Z
          RM2(3,1,6) = -3*RN2*X*Y*Z
          RM2(1,2,1) = 3*RN2*SYZ*X*Y
          RM2(1,2,2) = RN2*X*Y*(6 - 9*X2 + 3*X4 - 8*Y2 + 8*X2*Y2 - 
     1                 3*X4*Y2)/ SYZ5
          RM2(1,2,3) = RN2*X*Y*(2 - 3*X2 + X4 - 8*Z2 + 8*X2*Z2 - 
     1                 3*X4*Z2)/ SYZ5
          RM2(1,2,4) = RN2*(-1 + X2 + 2*Y2 - 3*X2*Y2)/SYZ
          RM2(1,2,5) = RN2*(2 - 3*X2)*Y*Z/SYZ
          RM2(1,2,6) = RN2*X*(2 - 3*X2 + X4 - 8*Y2 + 8*X2*Y2 - 
     1                 3*X4*Y2)*Z/ SYZ5
          RM2(2,2,1) = 0
          RM2(2,2,2) = RN2*(-1 + X2 + 3*Y2)*Z/SYZ5
          RM2(2,2,3) = -3*RN2*Y2*Z/SYZ5
          RM2(2,2,4) = 0
          RM2(2,2,5) = 0
          RM2(2,2,6) = RN2*Y*(-1 + X2 + 3*Z2)/SYZ5
          RM2(3,2,1) = RN2*(SYZ2 - 2*X2)*Y
          RM2(3,2,2) = 3*RN2*SXZ2*Y
          RM2(3,2,3) = RN2*Y*(SXY2 - 2*Z2)
          RM2(3,2,4) = RN2*X*(SXZ2 - 2*Y2)
          RM2(3,2,5) = -3*RN2*X*Y*Z
          RM2(3,2,6) = RN2*(SXZ2 - 2*Y2)*Z
          RM2(1,3,1) = 3*RN2*SYZ*X*Z
          RM2(1,3,2) = RN2*X*(2 - 3*X2 + X4 - 8*Y2 + 8*X2*Y2 -
     1                 3*X4*Y2)*Z/ SYZ5
          RM2(1,3,3) = RN2*X*Z*(6 - 9*X2 + 3*X4 - 8*Z2 + 8*X2*Z2 -
     1                 3*X4*Z2)/ SYZ5
          RM2(1,3,4) = RN2*(2 - 3*X2)*Y*Z/SYZ
          RM2(1,3,5) = RN2*(-1 + X2 + 2*Z2 - 3*X2*Z2)/SYZ
          RM2(1,3,6) = RN2*X*Y*(2 - 3*X2 + X4 - 8*Z2 + 8*X2*Z2 - 
     1                 3*X4*Z2)/ SYZ5
          RM2(2,3,1) = 0
          RM2(2,3,2) = 3*RN2*Y*Z2/SYZ5
          RM2(2,3,3) = RN2*Y*(Y2 - 2*Z2)/SYZ5
          RM2(2,3,4) = 0
          RM2(2,3,5) = 0
          RM2(2,3,6) = RN2*Z*(-2*Y2 + Z2)/SYZ5
          RM2(3,3,1) = RN2*(SYZ2 - 2*X2)*Z
          RM2(3,3,2) = RN2*(SXZ2 - 2*Y2)*Z
          RM2(3,3,3) = 3*RN2*SXY2*Z
          RM2(3,3,4) = -3*RN2*X*Y*Z
          RM2(3,3,5) = RN2*X*(SXY2 - 2*Z2)
          RM2(3,3,6) = RN2*Y*(SXY2 - 2*Z2)
      ELSE
C
C         Case B (no zero Sxz) matrix
C         Rotation matrix as a whole:
C
C         x y       z
C         ---      ---
C         Sxz      Sxz      -x
C
C
C         -Sxz     0        -y
C
C         y z         x
C         ---      -(---)
C         Sxz        Sxz    -z
C
C
          RM(1,1) = X*Y/SXZ
          RM(2,1) = Z/SXZ
          RM(3,1) = -X
          RM(1,2) = -SXZ
          RM(2,2) = 0
          RM(3,2) = -Y
          RM(1,3) = Y*Z/SXZ
          RM(2,3) = -(X/SXZ)
          RM(3,3) = -Z
          IF( IMOD .LE. 0 ) RETURN
C
          SXZ3 = SXZ * SXZ2
          RM1(1,1,1) = RN*Y*(-X4 + Y2*Z2 + Z4)/SXZ3
          RM1(1,1,2) = RN*SXZ*X
          RM1(1,1,3) = RN*X*Y*(-2 + Y2)*Z/SXZ3
          RM1(2,1,1) = -(RN*X*Z/SXZ3)
          RM1(2,1,2) = 0
          RM1(2,1,3) = RN*X2/SXZ3
          RM1(3,1,1) = -(RN*SYZ2)
          RM1(3,1,2) = RN*X*Y
          RM1(3,1,3) = RN*X*Z
          RM1(1,2,1) = -(RN*X*Y2/SXZ)
          RM1(1,2,2) = RN*SXZ*Y
          RM1(1,2,3) = -(RN*Y2*Z/SXZ)
          RM1(2,2,1) = 0
          RM1(2,2,2) = 0
          RM1(2,2,3) = 0
          RM1(3,2,1) = RN*X*Y
          RM1(3,2,2) = -(RN*SXZ2)
          RM1(3,2,3) = RN*Y*Z
          RM1(1,3,1) = RN*X*Y*(-2 + Y2)*Z/SXZ3
          RM1(1,3,2) = RN*SXZ*Z
          RM1(1,3,3) = RN*Y*(X4 + X2*Y2 - Z4)/SXZ3
          RM1(2,3,1) = -(RN*Z2/SXZ3)
          RM1(2,3,2) = 0
          RM1(2,3,3) = RN*X*Z/SXZ3
          RM1(3,3,1) = RN*X*Z
          RM1(3,3,2) = RN*Y*Z
          RM1(3,3,3) = -(RN*SXY2)
          IF( IMOD .LE. 1 ) RETURN
C
          SXZ5 = SXZ3 * SXZ2
          RM2(1,1,1) = RN2*X*Y*(-6 + 8*X2 + 9*Y2 - 8*X2*Y2 - 3*Y4 + 
     1                 3*X2*Y4)/SXZ5
          RM2(1,1,2) = -3*RN2*SXZ*X*Y
          RM2(1,1,3) = RN2*X*Y*(-2 + 3*Y2 - Y4 + 8*Z2 - 8*Y2*Z2 + 
     1                 3*Y4*Z2)/ SXZ5
          RM2(1,1,4) = RN2*(1 - 2*X2 - Y2 + 3*X2*Y2)/SXZ
          RM2(1,1,5) = RN2*Y*(-2 + 8*X2 + 3*Y2 - 8*X2*Y2 - Y4 + 
     1                 3*X2*Y4)*Z/ SXZ5
          RM2(1,1,6) = RN2*X*(-2 + 3*Y2)*Z/SXZ
          RM2(2,1,1) = RN2*(-1 + 3*X2 + Y2)*Z/SXZ5
          RM2(2,1,2) = 0
          RM2(2,1,3) = -3*RN2*X2*Z/SXZ5
          RM2(2,1,4) = 0
          RM2(2,1,5) = RN2*X*(-1 + Y2 + 3*Z2)/SXZ5
          RM2(2,1,6) = 0
          RM2(3,1,1) = 3*RN2*SYZ2*X
          RM2(3,1,2) = RN2*X*(SXZ2 - 2*Y2)
          RM2(3,1,3) = RN2*X*(SXY2 - 2*Z2)
          RM2(3,1,4) = RN2*(SYZ2 - 2*X2)*Y
          RM2(3,1,5) = RN2*(SYZ2 - 2*X2)*Z
          RM2(3,1,6) = -3*RN2*X*Y*Z
          RM2(1,2,1) = RN2*Y2*(-1 + 4*X2 + Y2 - 3*X2*Y2)/SXZ3
          RM2(1,2,2) = RN2*SXZ*(SXZ2 - 2*Y2)
          RM2(1,2,3) = RN2*Y2*(-1 + Y2 + 4*Z2 - 3*Y2*Z2)/SXZ3
          RM2(1,2,4) = RN2*X*Y*(-2 + 3*Y2)/SXZ
          RM2(1,2,5) = RN2*X*(4 - 3*Y2)*Y2*Z/SXZ3
          RM2(1,2,6) = RN2*Y*(-2 + 3*Y2)*Z/SXZ
          RM2(2,2,1) = 0
          RM2(2,2,2) = 0
          RM2(2,2,3) = 0
          RM2(2,2,4) = 0
          RM2(2,2,5) = 0
          RM2(2,2,6) = 0
          RM2(3,2,1) = RN2*(SYZ2 - 2*X2)*Y
          RM2(3,2,2) = 3*RN2*SXZ2*Y
          RM2(3,2,3) = RN2*Y*(SXY2 - 2*Z2)
          RM2(3,2,4) = RN2*X*(SXZ2 - 2*Y2)
          RM2(3,2,5) = -3*RN2*X*Y*Z
          RM2(3,2,6) = RN2*(SXZ2 - 2*Y2)*Z
          RM2(1,3,1) = RN2*Y*(-2 + 8*X2 + 3*Y2 - 8*X2*Y2 - Y4 + 
     1                 3*X2*Y4)*Z/ SXZ5
          RM2(1,3,2) = -3*RN2*SXZ*Y*Z
          RM2(1,3,3) = RN2*Y*Z*(-6 + 9*Y2 - 3*Y4 + 8*Z2 - 8*Y2*Z2 + 
     1                 3*Y4*Z2)/SXZ5
          RM2(1,3,4) = RN2*X*(-2 + 3*Y2)*Z/SXZ
          RM2(1,3,5) = RN2*X*Y*(-2 + 3*Y2 - Y4 + 8*Z2 - 8*Y2*Z2 + 
     1                 3*Y4*Z2)/ SXZ5
          RM2(1,3,6) = RN2*(1 - Y2 - 2*Z2 + 3*Y2*Z2)/SXZ
          RM2(2,3,1) = 3*RN2*X*Z2/SXZ5
          RM2(2,3,2) = 0
          RM2(2,3,3) = RN2*X*(X2 - 2*Z2)/SXZ5
          RM2(2,3,4) = 0
          RM2(2,3,5) = RN2*Z*(-2*X2 + Z2)/SXZ5
          RM2(2,3,6) = 0
          RM2(3,3,1) = RN2*(SYZ2 - 2*X2)*Z
          RM2(3,3,2) = RN2*(SXZ2 - 2*Y2)*Z
          RM2(3,3,3) = 3*RN2*SXY2*Z
          RM2(3,3,4) = -3*RN2*X*Y*Z
          RM2(3,3,5) = RN2*X*(SXY2 - 2*Z2)
          RM2(3,3,6) = RN2*Y*(SXY2 - 2*Z2)

      ENDIF
      RETURN
      END
C
      SUBROUTINE PSRMD(IMOD)
C
C   Compute D-orbitals rotation matrix and it's derivatives
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IMOD        - Computation type
C                0 = Compute rotation matrix RMD
C                1 = Compute rotation matrix and first derivatives
C                2 = Compute rotation matrix with 1st and 2nd 
C                    derivatives
C
C   Accessed common blocks:
C
C      PSROT0, PSROT1, PSROT2 - rotation matrix RM and, depending
C          on IMOD, it's derivatives RM1 and RM2 are used.
C  
C   Modified common blocks:
C
C      PSROT0, PSROT1, PSROT2 - rotation matrix RMD and, depending
C          on IMOD, it's derivatives RMD1 and RMD2 are computed.
C          D orbitals are in order X2Y2, XZ, Z2, YZ, XY
C          First derivatives of rotation matrix are in X,Y,Z order
C          Second derivatives are in XX, YY, ZZ, XY, XZ, YZ order
C
C   Local storage:
C
C          150 DOUBLE PRECISION cells
C
C   Module logic:
C
C          Basic subroutine PSRMDX is used to compute derivative
C          components of D-rotation matrix (or, as the degenerate
C          case, D-rotation matrix itself), which are then combined 
C          into actual derivatives by PSGVML.
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON 
     ./PSROT0/ RM (3,3), RMD (5,5), RMPP(6,6), RMPD(15,15), RMDD(15,15)
     ./PSROT1/ RM1(3,3,3), RMD1(5,5,3), RMPP1(6,6,3), RMPD1(15,15,3), 
     .         RMDD1(15,15,3)
     ./PSROT2/ RM2(3,3,6), RMD2(5,5,6), RMPP2(6,6,6), RMPD2(15,15,6), 
     .         RMDD2(15,15,6)
      SAVE /PSROT0/, /PSROT1/, /PSROT2/

      DIMENSION TEMP(5,5,6)
      EXTERNAL PSRMDX

      CALL PSGVML(IMOD,PSRMDX,9,9,25,RM,RM,RMD,RM1,RM1,RMD1,
     1            RM2,RM2,RMD2,TEMP)
      RETURN
      END
C
      SUBROUTINE PSRMDX(A,B,D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(3,3),B(3,3),D(5,5)
      PARAMETER (TWO =2.0D0)
      PARAMETER (HALF=0.5D0)
      PARAMETER (RS12=0.28867513459481288225457439025097872782380088D0)
      PARAMETER (RS3 =0.57735026918962576450914878050195745564760175D0)
C
C   First, few handy functions...
C
      CXY(I,J) = A(1,I)*B(2,J) + A(2,I)*B(1,J)
      CXZ(I,J) = A(1,I)*B(3,J) + A(3,I)*B(1,J)
      CYZ(I,J) = A(2,I)*B(3,J) + A(3,I)*B(2,J)
      CX2(I,J) = A(1,I)*B(1,J) - A(2,I)*B(2,J)
      CZ2(I,J) = RS3*(TWO*A(3,I)*B(3,J) 
     1                - A(1,I)*B(1,J) - A(2,I)*B(2,J))
C
C   And now, code itself - XY, XZ and YZ are directly computed in terms
C   of out auxiliary functions, Z2 and X2-Y2 are a bit more complex.
C

C  X2-Y2
      D(1,1) = HALF * ( CX2(1,1) - CX2(2,2) ) 
      D(2,1) = HALF * ( CXZ(1,1) - CXZ(2,2) ) 
      D(3,1) = HALF * ( CZ2(1,1) - CZ2(2,2) ) 
      D(4,1) = HALF * ( CYZ(1,1) - CYZ(2,2) ) 
      D(5,1) = HALF * ( CXY(1,1) - CXY(2,2) )
C  XZ
      D(1,2) = CX2(1,3) 
      D(2,2) = CXZ(1,3) 
      D(3,2) = CZ2(1,3) 
      D(4,2) = CYZ(1,3) 
      D(5,2) = CXY(1,3)
C  Z2
      D(1,3) = RS12 * ( TWO * CX2(3,3) - CX2(1,1) - CX2(2,2) ) 
      D(2,3) = RS12 * ( TWO * CXZ(3,3) - CXZ(1,1) - CXZ(2,2) ) 
      D(3,3) = RS12 * ( TWO * CZ2(3,3) - CZ2(1,1) - CZ2(2,2) ) 
      D(4,3) = RS12 * ( TWO * CYZ(3,3) - CYZ(1,1) - CYZ(2,2) ) 
      D(5,3) = RS12 * ( TWO * CXY(3,3) - CXY(1,1) - CXY(2,2) )
C  YZ
      D(1,4) = CX2(2,3) 
      D(2,4) = CXZ(2,3) 
      D(3,4) = CZ2(2,3) 
      D(4,4) = CYZ(2,3) 
      D(5,4) = CXY(2,3)
C  XY
      D(1,5) = CX2(1,2) 
      D(2,5) = CXZ(1,2) 
      D(3,5) = CZ2(1,2) 
      D(4,5) = CYZ(1,2) 
      D(5,5) = CXY(1,2)
      RETURN
      END
C
      SUBROUTINE PSRMPD(IMOD)
C
C   Compute PD-distribution rotation matrix and it's derivatives
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IMOD        - Computation type
C                0 = Compute rotation matrix RMPD
C                1 = Compute rotation matrix and first derivatives
C                    (RMPD and RMPD1)
C                2 = Compute rotation matrix with 1st and 2nd 
C                    derivatives (RMPD,RMPD1,RMPD2)
C
C   Accessed common blocks:
C
C      PSROT0, PSROT1, PSROT2 - rotation matrices RM, RMD and, 
C      depending on IMOD, respective derivatives RM1, RMD1, RM2
C      and RMD2 are used.
C  
C   Modified common blocks:
C
C      PSROT0, PSROT1, PSROT2 - rotation matrix RMPD and, depending
C          on IMOD, it's derivatives RMPD1 and RMPD2 are computed.
C          D orbitals are in order X2Y2, XZ, Z2, YZ, XY
C          First derivatives of rotation matrix are in X,Y,Z order
C          Second derivatives are in XX, YY, ZZ, XY, XZ, YZ order
C
C   Local storage:
C
C      1350 DOUBLE PRECISION cells
C
C   Module logic:
C
C          None :)
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON 
     ./PSROT0/ RM (3,3), RMD (5,5), RMPP(6,6), RMPD(15,15), RMDD(15,15)
     ./PSROT1/ RM1(3,3,3), RMD1(5,5,3), RMPP1(6,6,3), RMPD1(15,15,3), 
     .         RMDD1(15,15,3)
     ./PSROT2/ RM2(3,3,6), RMD2(5,5,6), RMPP2(6,6,6), RMPD2(15,15,6), 
     .         RMDD2(15,15,6)
      SAVE /PSROT0/, /PSROT1/, /PSROT2/

      DIMENSION TEMP(15,15,6)
      EXTERNAL PSRMPE

      CALL PSGVML(IMOD,PSRMPE,9,25,225,RM,RMD,RMPD,RM1,RMD1,RMPD1,
     1            RM2,RMD2,RMPD2,TEMP)
      
      RETURN
      END
C
      SUBROUTINE PSRMPE(P,D,PD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION P(3,3), D(5,5), PD(15,15)

      DO 10 ID=1,5
          DO 20 IP=1,3
              DO 30 JD=1,5
                  DO 40 JP=1,3
                      PD((JD-1)*3+JP,(ID-1)*3+IP) = P(JP,IP)*D(JD,ID)
   40             CONTINUE
   30         CONTINUE
   20     CONTINUE
   10 CONTINUE
      RETURN
      END
C 
      SUBROUTINE PSRMPP(IMOD)
C
C   Compute PP-distribution rotation matrix and it's derivatives
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IMOD        - Computation type
C                0 = Compute rotation matrix RMPP
C                1 = Compute rotation matrix and first derivatives
C                    (RMPP and RMPP1)
C                2 = Compute rotation matrix with 1st and 2nd 
C                    derivatives (RMPP,RMPP1,RMPP2)
C
C   Accessed common blocks:
C
C      PSROT0, PSROT1, PSROT2 - rotation matrix RM and, 
C      depending on IMOD, derivatives matrices RM1 and RM2
C
C   Local storage:
C
C      220 DOUBLE PRECISION cells
C  
C   Modified common blocks:
C
C      PSROT0, PSROT1, PSROT2 - rotation matrix RMPP and, depending
C          on IMOD, it's derivatives RMPP1 and RMPP2 are computed.
C          First derivatives of rotation matrix are in X,Y,Z order
C
C   Module logic:
C
C          None :)
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON 
     ./PSROT0/ RM (3,3), RMD (5,5), RMPP(6,6), RMPD(15,15), RMDD(15,15)
     ./PSROT1/ RM1(3,3,3), RMD1(5,5,3), RMPP1(6,6,3), RMPD1(15,15,3), 
     .         RMDD1(15,15,3)
     ./PSROT2/ RM2(3,3,6), RMD2(5,5,6), RMPP2(6,6,6), RMPD2(15,15,6), 
     .         RMDD2(15,15,6)
      SAVE /PSROT0/, /PSROT1/, /PSROT2/

      DIMENSION TEMP(6,6,6)
      EXTERNAL PSRMPQ

      CALL PSGVML(IMOD,PSRMPQ,9,9,36,RM,RM,RMPP,RM1,RM1,RMPP1,
     1            RM2,RM2,RMPP2,TEMP)
      
      RETURN
      END
C
      SUBROUTINE PSRMPQ(P,Q,PP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION P(3,3), Q(3,3), PP(6,6)

      IND2 = 1
      DO 10 J2 = 1,3
          DO 20 I2 = 1, J2
              IND1 = 1
              DO 30 J1 = 1, 3
                  DO 40 I1 = 1, J1-1
                      PP(IND1,IND2)=P(I1,I2)*Q(J1,J2)+P(J1,I2)*Q(I1,J2)
                      IND1=IND1+1
   40             CONTINUE
                  PP(IND1,IND2)=P(I1,I2)*Q(I1,J2)
                  IND1 = IND1 + 1
   30         CONTINUE
              IND2 = IND2 + 1
   20     CONTINUE
   10 CONTINUE
      RETURN
      END
C 
      SUBROUTINE PSRMDD(IMOD)
C
C   Compute DD-distribution rotation matrix and it's derivatives
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IMOD        - Computation type
C                0 = Compute rotation matrix RMDD
C                1 = Compute rotation matrix and first derivatives
C                    (RMDD and RMDD1)
C                2 = Compute rotation matrix with 1st and 2nd 
C                    derivatives (RMDD,RMDD1,RMDD2)
C
C   Accessed common blocks:
C
C      PSROT0, PSROT1, PSROT2 - rotation matrix RMD and, 
C      depending on IMOD, it's derivatives RMD1 and RMD2
C  
C   Modified common blocks:
C
C      PSROT0, PSROT1, PSROT2 - rotation matrix RMDD and, depending
C          on IMOD, it's derivatives RMDD1 and RMDD2 are computed.
C          D orbitals are in order X2Y2, XZ, Z2, YZ, XY
C          First derivatives of rotation matrix are in X,Y,Z order
C          Second derivatives are in XX, YY, ZZ, XY, XZ, YZ order
C
C   Local storage:
C
C      1350 DOUBLE PRECISION cells
C
C   Module logic:
C
C          See PSRMPP
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON 
     ./PSROT0/ RM (3,3), RMD (5,5), RMPP(6,6), RMPD(15,15), RMDD(15,15)
     ./PSROT1/ RM1(3,3,3), RMD1(5,5,3), RMPP1(6,6,3), RMPD1(15,15,3), 
     .         RMDD1(15,15,3)
     ./PSROT2/ RM2(3,3,6), RMD2(5,5,6), RMPP2(6,6,6), RMPD2(15,15,6), 
     .         RMDD2(15,15,6)
      SAVE /PSROT0/, /PSROT1/, /PSROT2/

      DIMENSION TEMP(15,15,6)
      EXTERNAL PSRMDE

      CALL PSGVML(IMOD,PSRMDE,25,25,225,RMD,RMD,RMDD,RMD1,RMD1,RMDD1,
     1            RMD2,RMD2,RMDD2,TEMP)
      
      RETURN
      END
C
      SUBROUTINE PSRMDE(D,E,DD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION D(5,5), E(5,5), DD(15,15)

      IND2 = 1
      DO 10 J2 = 1,5
          DO 20 I2 = 1, J2
              IND1 = 1
              DO 30 J1 = 1, 5
                  DO 40 I1 = 1, J1-1
                      DD(IND1,IND2)=D(I1,I2)*E(J1,J2)+D(J1,I2)*E(I1,J2)
                      IND1=IND1+1
   40             CONTINUE
                  DD(IND1,IND2)=D(I1,I2)*E(I1,J2)
                  IND1 = IND1 + 1
   30         CONTINUE
              IND2 = IND2 + 1
   20     CONTINUE
   10 CONTINUE
      RETURN
      END
C 
