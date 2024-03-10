C     ******************************************************************
C
C     General-purpose matrix handling utilities.
C
C     ******************************************************************
      SUBROUTINE PSPK2U(IROWS,A,B,LDB)
C
C   Unpack matrix held in a packed storage into upper triangle
C   of a general square matrix B.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IROWS  - Number of rows in matrix A
C      A      - Matrix in packed form
C      B      - Output unpacked matrix
C      LDB    - Leading dimension of matrix B
C
C   Accessed common blocks:
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C
C      Probably it would be a good idea to switch to the
C      BLAS1 calls starting at some (large) column length
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION A(*), B(LDB,*)
C
      IIN = 0
      DO 100 I=1,IROWS
          DO 98 J=1,I
              B(J,I) = A(IIN+J)
   98     CONTINUE
          IIN = IIN + I
  100 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSU2PK(IROWS,A,LDA,B)
C
C   Pack upper triangle of symmetric matrix A to the
C   matrix B.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IROWS  - Number of rows in matrix A
C      A      - Matrix in unpacked form, only upper
C               triangle is used
C      LDA    - Leading dimension of A
C      B      - Output packed matrix
C
C   Accessed common blocks:
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C
C      Probably it would be a good idea to switch to the
C      BLAS1 calls starting at some (large) column length
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION A(LDA,*), B(*)
C
      IOUT = 0
      DO 100 I=1,IROWS
          DO 98 J=1,I
              B(IOUT+J) = A(J,I)
   98     CONTINUE
          IOUT = IOUT + I
  100 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSU2PS(IROWS,A,LDA)
C
C   Pack upper triangle of symmetric matrix A in place.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IROWS  - Number of rows in matrix A
C      A      - Input: upper triangle of the matrix in unpacked form
C               Output: matrix A in the packed form.
C      LDA    - Leading dimension of A
C
C   Accessed common blocks:
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C
C      Probably it would be a good idea to switch to the
C      BLAS1 calls starting at some (large) column length
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
C
      DIMENSION A(*)
C
      IF( LDA.LT.IROWS ) THEN
          WRITE(NB6,10000) IROWS, LDA
          STOP 'PSU2PS'
      ENDIF
C
      IOUT = 1
      DO 100 I=2,IROWS
CVD$L NODEPCHK
C*$*ASSERT RELATION (IOUT+I.LT.1+(I-1)*LDA)
          DO 98 J=1,I
              A(IOUT+J) = A(J+(I-1)*LDA)
   98     CONTINUE
          IOUT = IOUT + I
  100 CONTINUE
C
      RETURN
10000 FORMAT(' MATRIX HAS ',I6,' ROWS, BUT LEADING DIMENSION IS ONLY ',
     .       I6,' IN PSU2PS.')
      END
C
      SUBROUTINE PSU2UL(Y,LDY,NORB)
C
C   Complement diagonal upper half-block of the symmetrix
C   matrix by lower counterpart
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      Y      - Symmetric matrix to complement
C      LDY    - Leading dimension of Y
C      NORB   - Number of columns in matrix
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
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(LDY,*)
C
      DO 100 I=2,NORB
          DO 98 J=1,I-1
              Y(I,J) = Y(J,I)
   98     CONTINUE
  100 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSU2AS(Y,LDY,NORB)
C
C   Complement diagonal upper half-block of the antisymmetric
C   matrix by lower counterpart
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      Y      - Antisymmetric matrix to complement
C      LDY    - Leading dimension of Y
C      NORB   - Number of columns in matrix
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
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.D0)
      DIMENSION Y(LDY,*)
C
      DO 100 I=1,NORB
          DO 98 J=1,I-1
              Y(I,J) = -Y(J,I)
   98     CONTINUE
          Y(I,I) = ZERO
  100 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSDSCP(IROWS,ISIZE,A)
C
C   Multiply off-diagonal elements of matrix held in packed storage
C   by 2.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IROWS  - Number of rows in matrix A
C      ISIZE  - Number of elements in linear form of matrix A.
C               It is caller's responsibility to ensure that
C               ISIZE = (IROWS*(IROWS+1))/2
C      A      - Matrix to scale.
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
C      Double whole array, then divide diagonal elements by 2.
C      Trivially oprimise out case ISIZE=1 (which should never
C      occure in practice)
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (TWO =2.0D0)
      PARAMETER (HALF=0.5D0)
      DIMENSION A(ISIZE)
C
      IF( ISIZE.EQ.1 ) RETURN
      CALL DSCAL(ISIZE-1,TWO,A(2),1)
      IA = 3
      DO 100 I=2,IROWS
          A(IA) = A(IA) * HALF
          IA    = IA + I + 1
  100 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSSCTN(NVAR,A)
C
C   Multiply elements of a two-particle density matrix in packed
C   storage by their symmetry number.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      NVAR   - Number of rows in the tensor.
C      A      - Tensor in packed storage.
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
C   Bugs:
C
C      This code is not particularly efficient.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (TWO  =2.D0)
      PARAMETER (FOUR =4.D0)
      PARAMETER (EIGHT=8.D0)
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
C
      DIMENSION A(*)
      IPACK(N) = (N*(N+1))/2
C
      IPOS = 1
      DO 900 L=1,NVAR
          DO 500 K=1,L-1
              DO 300 J=1,L-1
                  DO 200 I=1,J-1
C                    (IJ|KL)
                      A(IPOS) = EIGHT*A(IPOS)
                      IPOS = IPOS + 1
  200             CONTINUE
C                (JJ|KL)
                  A(IPOS) = FOUR*A(IPOS)
                  IPOS = IPOS + 1
  300         CONTINUE
C            J = L
              DO 400 I=1,K-1
C                (IL|KL)
                  A(IPOS) = EIGHT*A(IPOS)
                  IPOS = IPOS + 1
  400         CONTINUE
C            (KL|KL)
              A(IPOS) = FOUR*A(IPOS)
              IPOS = IPOS + 1
  500     CONTINUE
C        K = L
          DO 700 J=1,L-1
              DO 600 I=1,J-1
C                (IJ|LL)
                  A(IPOS) = FOUR*A(IPOS)
                  IPOS = IPOS + 1
  600         CONTINUE
C            (JJ|LL)
              A(IPOS) = TWO*A(IPOS)
              IPOS = IPOS + 1
  700     CONTINUE
C         J = L
          DO 800 I=1,L-1
C            (IL|LL)
              A(IPOS) = FOUR*A(IPOS)
              IPOS = IPOS + 1
  800     CONTINUE
C        (LL|LL)
          IPOS = IPOS + 1
  900 CONTINUE
C
      IF( IPOS-1.NE.IPACK(IPACK(NVAR)) ) THEN
          WRITE(NB6,11010) IPOS-1, IPACK(IPACK(NVAR))
          STOP 'PSSCTN'
      ENDIF
C
      RETURN
11010 FORMAT(' LOGIC FAILURE IN PSSCTN: ACTUAL INDEX = ',I10,
     .       ' EXPECTED INDEX = ',I10)
      END
C
      SUBROUTINE PSZRV(LEN,X)
C
C   Zero out a vector.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      LEN    - Length of the vector
C      X      - (Output) LEN zeros.
C
C   Accessed common blocks:
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.D0)
C
      DIMENSION X(LEN)
C
      DO 100 I=1,LEN 
          X(I) = ZERO
  100 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSZRM(L1,L2,X,LDX)
C
C   Zero out a matrix.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      L1     - Number of rows in X.
C      L2     - Number of columns in X.
C      X      - (Output) zero matrix.
C      LDX    - Leading dimension of X.
C
C   Accessed common blocks:
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.D0)
C
      DIMENSION X(LDX,L2)
C
      DO 200 I2=1,L2
          DO 100 I1=1,L1
              X(I1,I2) = ZERO
  100     CONTINUE
  200 CONTINUE
C
      RETURN
      END
