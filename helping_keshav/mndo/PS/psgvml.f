C    ******************************************************************
C
C     Multiplication of function values with derivatives.
C
C     ******************************************************************
      SUBROUTINE PSGVML(IMOD,FUNC,NINF,NING,NOUT,
     1                  F,G,R,F1,G1,R1,F2,G2,R2,TEMP)
C
C   Generalized 2nd order vector multiplication routine.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IMOD   - Computation type
C               0 = compute generalized product
C               1 = .... and first derivatives
C               2 = ......... and second derivatives
C      FUNC   - Vector multiplication procedure
C      NINF   - Size of single component of F
C      NING   - Size of single component of G
C      NOUT   - Size of single component of R
C      F,G    - Input vectors
C      R      - Output vector
C      F1,G1  - First derivatives of input vectors
C               These arrays are not accessed if
C               IMOD .LE. 0
C      R1     - First derivative of output
C               This array is not touched if IMOD .LE.0
C      F2,G2  - Second derivatives of input vectors
C               Not accessed if IMOD.LE.1
C      R2     - Second derivative of output vector
C               Not accessed if IMOD.LE.1
C      TEMP   - Scratch space. 
C               Not touched if IMOD.LE.0
C               At least 3*NOUT if IMOD.EQ.1
C               At least 6*NOUT if IMOD.EQ.2
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
C      Generalized zero-order vector multiplication routine FUNC
C      is used to compute components of derivative expressions,
C      which are then combined according to standard differentiation
C      rules:
C
C          D (f*g,x)   = D (f,x)*g + f*D(g,x) 
C          D2(f*g,x,y) = D2(f,x,y)*g + D(f,x)*D(g,y) + D(f,y)*D(g,x) +
C                        f*D2(g,x,y)
C          D2(f*g,x,x) = D2(f,x,x)*g + 2*D(f,x)*D(g,x) + f*D2(g,x,x)
C
C      For derivatives, as much components of expressions are first
C      computed to the TEMP array and then combined with a single 
C      call to DAXPY. 
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ONE=1.0D0)
      PARAMETER (TWO=2.0D0)

      DIMENSION F(NINF),F1(NINF,3),F2(NINF,6)
      DIMENSION G(NING),G1(NING,3),G2(NING,6)
      DIMENSION R(NOUT),R1(NOUT,3),R2(NOUT,6)
      DIMENSION TEMP(NOUT,6)
      EXTERNAL  FUNC

      CALL FUNC(F,G,R)
      IF( IMOD.LE.0 ) RETURN
C
C   First derivatives
C
      DO 10 I=1,3
          CALL FUNC(F,G1(1,I),R1  (1,I))
          CALL FUNC(F1(1,I),G,TEMP(1,I))
   10 CONTINUE
      CALL DAXPY(3*NOUT,ONE,TEMP,1,R1,1)
C
      IF( IMOD.LE.1 ) RETURN
C
C   Second derivatives
C
C   Fxy*G and F*Gxy contributions
C
      DO 20 I=1,6
          CALL FUNC(F,G2(1,I),R2  (1,I))
          CALL FUNC(F2(1,I),G,TEMP(1,I))
   20 CONTINUE
      CALL DAXPY(6*NOUT,ONE,TEMP,1,R2,1)
C
C   Fx*Gx contributions - diagonal second derivatives
C
      DO 30 I=1,3
          CALL FUNC(F1(1,I),G1(1,I),TEMP(1,I))
   30 CONTINUE
      CALL DAXPY(3*NOUT,TWO,TEMP,1,R2,1)
C
C   Fx*Gy contributions - off-diagonal second derivatives
C
      CALL FUNC(F1(1,1),G1(1,2),TEMP(1,1))
      CALL FUNC(F1(1,1),G1(1,3),TEMP(1,2))
      CALL FUNC(F1(1,2),G1(1,3),TEMP(1,3))
      CALL DAXPY(3*NOUT,ONE,TEMP,1,R2(1,4),1)
C
C   Fy*Gx contributions - off-diagonal second derivatives
C
      CALL FUNC(F1(1,2),G1(1,1),TEMP(1,1))
      CALL FUNC(F1(1,3),G1(1,1),TEMP(1,2))
      CALL FUNC(F1(1,3),G1(1,2),TEMP(1,3))
      CALL DAXPY(3*NOUT,ONE,TEMP,1,R2(1,4),1)
      RETURN
      END
C
      SUBROUTINE PSXML(IMOD,A,B,OUT)
C
C   Multiply two function values along with derivatives 
C   up to the second order according to chaining rules.
C   This code is not intended to be called repeatedly,
C   please use PSGVML to multiply matrices.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IMOD   - Computation type
C               0 = compute product
C               1 = .... and first derivatives
C               2 = ......... and second derivatives
C      A,B    - Function values to multiply
C      C      - Output result
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
C      Computations are performed according to standard 
C      differentiation rules:
C
C          D (f*g,x)   = D (f,x)*g + f*D(g,x) 
C          D2(f*g,x,y) = D2(f,x,y)*g + D(f,x)*D(g,y) + D(f,y)*D(g,x) +
C                        f*D2(g,x,y)
C          D2(f*g,x,x) = D2(f,x,x)*g + 2*D(f,x)*D(g,x) + f*D2(g,x,x)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (TWO=2.0D0)

      DIMENSION A(10), B(10), OUT(10)

      OUT( 1) = A(1) * B(1)
      IF( IMOD.LE.0 ) RETURN
C
      OUT( 2) = A(1)*B(2) + A(2)*B(1)
      OUT( 3) = A(1)*B(3) + A(3)*B(1)
      OUT( 4) = A(1)*B(4) + A(4)*B(1)
      IF( IMOD.LE.1 ) RETURN
C
      OUT( 5) = A(1)*B(5) + TWO*A(2)*B(2) + A(5)*B(1)
      OUT( 6) = A(1)*B(6) + TWO*A(3)*B(3) + A(6)*B(1)
      OUT( 7) = A(1)*B(7) + TWO*A(4)*B(4) + A(7)*B(1)
      OUT( 8) = A(1)*B( 8) + A(2)*B(3) + A(3)*B(2) + A( 8)*B(1)
      OUT( 9) = A(1)*B( 9) + A(2)*B(4) + A(4)*B(2) + A( 9)*B(1)
      OUT(10) = A(1)*B(10) + A(3)*B(4) + A(4)*B(3) + A(10)*B(1)
C
      RETURN
      END
C
      SUBROUTINE PSXIV(IMOD,A,AI)
C
C   Computes 1/f(x) along with it's derivatives. This function
C   is not intended to be called repeatedly.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IMOD   - (Input) Computation type
C               0 = compute inverse
C               1 = .... and first derivatives
C               2 = ......... and second derivatives
C      A      - (Input) Function value and derivatives, in the
C               order: function, X,Y, and Z derivatives, and
C               XX, YY, ZZ, XY, XZ, and YZ second derivatives.
C      AI     - (Output) 1/A function value and derivatives.
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
C          1
C     f = ---
C          g
C
C                  (1,0)
C      (1,0)      g
C     f      = - ---
C                g*g
C
C                  (1,0)   (1,0)    (2,0)
C      (2,0)      g     * g        g
C     f      = 2 -----------   -  ---
C                   g*g*g         g*g
C
C                  (1,0)   (0,1)    (1,1)
C      (1,1)      g     * g        g
C     f      = 2 -----------   -  ---
C                   g*g*g         g*g
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ONE=1.D0,TWO=2.D0)

      DIMENSION A(10), AI(10)
C    Function value
      AI(1) = ONE/A(1)
      IF( IMOD.LE.0 ) RETURN
C    First derivatives
      SRG2  = -(AI(1)**2)
      AI(2) = SRG2*A(2)
      AI(3) = SRG2*A(3)
      AI(4) = SRG2*A(4)
      IF( IMOD.LE.1 ) RETURN
C    Second derivatives
      DRG3   = -TWO*SRG2*AI(1)
      AI( 5) = DRG3*A(2)**2 + SRG2*A(5)
      AI( 6) = DRG3*A(3)**2 + SRG2*A(6)
      AI( 7) = DRG3*A(4)**2 + SRG2*A(7)
      AI( 8) = DRG3*A(2)*A(3) + SRG2*A(8)
      AI( 9) = DRG3*A(2)*A(4) + SRG2*A(9)
      AI(10) = DRG3*A(3)*A(4) + SRG2*A(10)
C
      RETURN
      END
