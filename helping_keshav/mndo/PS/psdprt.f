C     ******************************************************************
C
C     Printing routines for analytical derivative code
C
C     ******************************************************************
      SUBROUTINE PSDPDR(NVARS,DERIV)
C
C   Print first derivatives.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      NVARS  - Number of elements of DERIV to print
C      DERIV  - Derivatives, at least NVARS elements
C
C   Accessed common blocks:
C
C      PSPRT  - Printing unit.
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (IWIDTH=8)
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
      DIMENSION DERIV(*)

      DO 100 I=1,NVARS,IWIDTH
          ITOP = MIN(NVARS,I+IWIDTH-1)
          WRITE(NB6,10000) (II,II=I,ITOP)
          WRITE(NB6,10010) (DERIV(II),II=I,ITOP)
          WRITE(NB6,10020)
  100 CONTINUE

      RETURN
10000 FORMAT( ' I = ', 8(I9,5X) )
10010 FORMAT( ' D = ', 8(F13.6,1X) )
10020 FORMAT()
      END
C
      SUBROUTINE PSDPFC(NVARS,FRC,LDF)
C
C   Print second derivatives matrix
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      NVARS  - Number of elements of DERIV to print
C      FRC    - Derivatives, at least NVARS elements
C      LDF    - Leading dimension of force matrix
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
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION FRC(LDF,*)
      CALL PSDPGM(NVARS,NVARS,FRC,LDF)
      RETURN
      END
C
      SUBROUTINE PSDPPM(NVARS,A)
C
C   Print symmetric matrix held in a packed storage
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      NVARS  - Number of rows of A to print
C      A      - Symmetric matrix in a packed form
C
C   Accessed common blocks:
C
C      PSPRT  - Printing unit.
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (IWIDTH=8)
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
      DIMENSION A(*)
C
      SUM = 0.0D0
      DO 100 J=1,NVARS,IWIDTH
          JTOP = MIN(NVARS,J+IWIDTH-1)
          WRITE(NB6,10000) (JJ,JJ=J,JTOP)
          DO 10 I=J,JTOP
              SUM = SUM + A( (I*(I+1))/2 )
   10     CONTINUE
          DO 90 I=J,NVARS
              J1TOP = MIN(I,JTOP)
              IND   = (I*(I+1))/2 - I
              WRITE(NB6,10010) I, (A(IND+JJ),JJ=J,J1TOP)
   90     CONTINUE
          WRITE(NB6,10020)
  100 CONTINUE
      WRITE(NB6,10030) SUM/DBLE(NVARS)
      RETURN
10000 FORMAT( '        ', 8(I9,5X) )
10010 FORMAT( I7, 1X, 8(F13.6,1X) )
10020 FORMAT()
10030 FORMAT( '        MEAN DIAGONAL VALUE IS ', F13.6 )
      END
C
      SUBROUTINE PSDPGM(NROWS,NCOLS,A,LDA)
C
C   Print general matrix
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      NROWS  - Number of rows to print
C      NCOLS  - Number of columns to print
C      A      - Matrix to print
C      LDA    - Leading dimension of matrix A
C
C   Accessed common blocks:
C
C      PSPRT  - Printing unit.
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (IWIDTH=8)
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
      DIMENSION A(LDA,*)
C
      DO 100 J=1,NCOLS,IWIDTH
          JTOP = MIN(NCOLS,J+IWIDTH-1)
          WRITE(NB6,10000) (JJ,JJ=J,JTOP)
          DO 90 I=1,NROWS
              WRITE(NB6,10010) I, (A(I,JJ),JJ=J,JTOP)
   90     CONTINUE
          WRITE(NB6,10020)
  100 CONTINUE

      RETURN
10000 FORMAT( '        ', 8(I9,5X) )
10010 FORMAT( I7, 1X, 8(F13.6,1X) )
10020 FORMAT()
      END
C
      SUBROUTINE PSDPGI(NROWS,NCOLS,IA,LDA)
C
C   Print general integer matrix
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      NROWS  - Number of rows to print
C      NCOLS  - Number of columns to print
C      IA     - Matrix to print
C      LDA    - Leading dimension of matrix IA
C
C   Accessed common blocks:
C
C      PSPRT  - Printing unit.
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (IWIDTH=8)
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
      DIMENSION IA(LDA,*)
C
      DO 100 J=1,NCOLS,IWIDTH
          JTOP = MIN(NCOLS,J+IWIDTH-1)
          WRITE(NB6,10000) (JJ,JJ=J,JTOP)
          DO 90 I=1,NROWS
              WRITE(NB6,10010) I, (IA(I,JJ),JJ=J,JTOP)
   90     CONTINUE
          WRITE(NB6,10020)
  100 CONTINUE

      RETURN
10000 FORMAT( '        ',  8(I6,':',1X) )
10010 FORMAT( I6, '/', 1X, 8(I7,1X) )
10020 FORMAT()
      END
C
      SUBROUTINE PSDPSU(NVARS,A,LDA)
      ENTRY PSDPAU(NVARS,A,LDA)
C
C   Print upper triangle of a symmetric matrix held in
C   square array.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      NVARS  - Number of rows to print
C      A      - Matrix to print
C      LDA    - Leading dimension of matrix A
C
C   Accessed common blocks:
C
C      PSPRT  - Printing unit.
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (IWIDTH=8)
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
      DIMENSION A(LDA,*)
C
C
      DO 100 J=1,NVARS,IWIDTH
          JTOP = MIN(NVARS,J+IWIDTH-1)
          WRITE(NB6,10000) (JJ,JJ=J,JTOP)
          DO 90 I=J,NVARS
              J1TOP = MIN(I,JTOP)
              WRITE(NB6,10010) I, (A(JJ,I),JJ=J,J1TOP)
   90     CONTINUE
          WRITE(NB6,10020)
  100 CONTINUE


      RETURN
10000 FORMAT( '        ', 8(I9,5X) )
10010 FORMAT( I7, 1X, 8(F13.6,1X) )
10020 FORMAT()
      END
