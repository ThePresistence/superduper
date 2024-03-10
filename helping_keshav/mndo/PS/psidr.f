C     ******************************************************************
C
C     Projection of derivative integrals in local coordinates
C     on the axes of the molecular coordinate system.
C
C     ******************************************************************
      SUBROUTINE PSIDR(IMOD,IORBA,IORBB,AINT1,OUT1,AINT2,OUT2)
C
C   Compute projections of derivatives of integrals in local
C   coordinate system on axes of molecular coordinate system
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IMOD   - Computation type
C           1 = Project first derivatives
C           2 = Project both first and second derivatives
C      IORBA  - Number of AOs on first atom, zero-based
C      IORBB  - Number of AOs on second atom, zero-based
C      AINT1  - Derivatives of integrals in local system
C      OUT1   - Computed projections of derivatives to
C               molecular coordinate system (that is, 
C               "molecular derivativies of local integrals"
C               if you wanna be cryptic :-)
C               First derivatives are in X,Y,Z order
C      AINT2  - Second derivatives of integrals in local
C               system
C      OUT2   - Computed projections of second derivatives
C               Second derivatives are in XX, YY, ZZ, XY, 
C               XZ, YZ order
C
C   Accessed common blocks:
C
C      PSROT0,
C      PSROT1 - Get basic rotation matrix RM and it's 
C               first derivatives (RM1) from here. Note
C               that we don't need RM2 - first derivatives
C               of rotation matrix are enougth to project
C               second derivatives.
C  
C   Modified common blocks:
C
C      PSDROT - Scratch common, used to pass projection
C               coefficients to sibling routines.
C
C   Local storage:
C
C          None.
C
C   Module logic:
C
C      First derivatives are computed according to
C      formulas:
C
C          F3D(1) = -( F(1) * RM(3,1), F(1) * RM(3,2), F(1) * RM(3,3)
C
C      Minus sign comes from the inversion of coordinate
C      system during transition from local to molecular
C      coordinate system. We take advantage of the fact
C      RM(-1) = RM(T) to obtain transformation coefficients.
C
C      Second derivatives are computed according to the following
C      formulas:
C
C          F3D(2,XX) = F(2) * RM(3,1)**2 - F(1) * RM(1,X)(3,1)
C          F3D(2,YY) = F(2) * RM(3,2)**2 - F(1) * RM(1,Y)(3,2)
C          F3D(2,ZZ) = F(2) * RM(3,3)**2 - F(1) * RM(1,Z)(3,3)
C          F3D(2,XY) = F(2) * RM(3,1)*RM(3,2) - F(1) * RM(1,Y)(3,1)
C          F3D(2,XZ) = F(2) * RM(3,1)*RM(3,3) - F(1) * RM(1,Z)(3,1)
C          F3D(2,YZ) = F(2) * RM(3,2)*RM(3,3) - F(1) * RM(1,Z)(3,2)
C
C      This routine serves as a dispatch for the following 
C      special-case procedures:
C
C          PSIDRS - first derivatives, atom A has only S orbitals
C          PSIDRP - first derivatives, atom A has S and P orbitals
C          PSIDRD - first derivatives, atom A has S,P and D orbitals
C
C          PSIDRT - second derivatives, atom A has only S orbitals
C          PSIDRQ - second derivatives, atom A has S and P orbitals
C          PSIDRE - second derivatives, atom A has S,P and D orbitals
C
C      It should minimize number of conditional branches and
C      thus number of pipeline stalls
C
C      Note that the gradient projection coefficients are now evaluated
C      in PSRM along with the rest of the transformation matrices.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION AINT1(461), AINT2(461)
      DIMENSION OUT1(461,3), OUT2(461,6)

C
C  Do not have derivatives to project, bail out
C
      IF(IMOD.LE.0) RETURN
C
      IF( IORBA.LE.0 ) THEN
          CALL PSIDRS(IORBB,AINT1,OUT1)
      ELSE IF( IORBA.LE.3 ) THEN
          CALL PSIDRP(IORBB,AINT1,OUT1)
      ELSE
          CALL PSIDRD(IORBB,AINT1,OUT1)
      ENDIF
C
      IF(IMOD.LE.1) RETURN
C
      IF( IORBA.LE.0 ) THEN
          CALL PSIDRT(IORBB,AINT1,AINT2,OUT2)
      ELSE IF( IORBA.LE.3 ) THEN
          CALL PSIDRQ(IORBB,AINT1,AINT2,OUT2)
      ELSE
          CALL PSIDRE(IORBB,AINT1,AINT2,OUT2)
      ENDIF
      RETURN
      END
