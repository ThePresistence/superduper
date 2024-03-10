C     ******************************************************************
C
C     Get machine dependent precision through LAPACK.
C
C     ******************************************************************
      DOUBLE PRECISION FUNCTION D1MACH (I)
C     *
C     Return floating point machine dependent constants through LAPACK.
C     Called from DEXINT (SLATEC library).
C     *
C     D1MACH can be used to obtain machine-dependent parameters for the
C     local machine environment.  It is a function subprogram with one
C     (input) argument, and can be referenced as follows:
C
C        D = D1MACH(I)
C
C     where I=1,...,5.  Original conventions:
C
C     D1MACH( 1) = B**(EMIN-1), the smallest positive magnitude.
C     D1MACH( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
C     D1MACH( 3) = B**(-T), the smallest relative spacing.
C     D1MACH( 4) = B**(1-T), the largest relative spacing.
C     D1MACH( 5) = LOG10(B)
C
C     Needed by DEXINT: D1MACH(4), D1MACH(5)
C     Only the needed constants are supplied.
C     The other constants are given as comments.
C     *
      LOGICAL FIRST
      DOUBLE PRECISION DMACH(5), BASE
      DOUBLE PRECISION DLAMCH
C     DOUBLE PRECISION DLAMC2
C     DOUBLE PRECISION EPS, RMIN, RMAX
C     INTEGER BETA, IMAX, IMIN, IT
C     LOGICAL LRND
      SAVE FIRST, DMACH
      DATA FIRST / .TRUE. /
C
      IF( FIRST ) THEN
         FIRST    = .FALSE.
         DMACH(1) = 0.0D0
         DMACH(2) = 0.0D0
         DMACH(3) = 0.0D0
         DMACH(4) = DLAMCH('E')
         BASE     = DLAMCH('B')
         DMACH(5) = DLOG10(BASE)
C        CALL DLAMC2( BETA, IT, LRND, EPS, IMIN, RMIN, IMAX, RMAX )
C        BASE     = BETA
C        DMACH(1) = BASE**( IMIN-1 )
C        DMACH(2) = BASE**IMAX*( 1- BASE**(-IT) )
C        DMACH(3) = BASE**(-IT)
      ENDIF
      IF( I.LT.1 .OR. I.GT.5 ) THEN
         CALL XERMSG ('SLATEC', 'D1MACH', 'I OUT OF BOUNDS', 1, 2)
      ELSE
         D1MACH = DMACH(I)
      ENDIF
      RETURN
      END
C     ******************************************************************
      INTEGER FUNCTION I1MACH (I)
C     *
C     Return integer machine dependent constants through LAPACK.
C     Called from DEXINT (SLATEC library).
C     *
C     I1MACH can be used to obtain machine-dependent parameters for the
C     local machine environment.  It is a function subprogram with one
C     (input) argument and can be referenced as follows:
C
C        K = I1MACH(I)
C
C     where I=1,...,16.  Original conventions:
C
C     I/O unit numbers:
C     I1MACH( 1) = the standard input unit.
C     I1MACH( 2) = the standard output unit.
C     I1MACH( 3) = the standard punch unit.
C     I1MACH( 4) = the standard error message unit.
C
C     Words:
C     I1MACH( 5) = the number of bits per integer storage unit.
C     I1MACH( 6) = the number of characters per integer storage unit.
C
C     Integers:
C     I1MACH( 7) = A, the base.
C     I1MACH( 8) = S, the number of base-A digits.
C     I1MACH( 9) = A**S - 1, the largest magnitude.
C
C     Floating-Point Numbers:
C     I1MACH(10) = B, the base.
C
C     Single-Precision:
C     I1MACH(11) = T, the number of base-B digits.
C     I1MACH(12) = EMIN, the smallest exponent E.
C     I1MACH(13) = EMAX, the largest exponent E.
C
C     Double-Precision:
C     I1MACH(14) = T, the number of base-B digits.
C     I1MACH(15) = EMIN, the smallest exponent E.
C     I1MACH(16) = EMAX, the largest exponent E.
C
C     Needed by DEXINT: I1MACH(4), I1MACH(15)
C     The needed constants are supplied.
C     Some other constants are given as comments.
C     *
      LOGICAL FIRST
      DOUBLE PRECISION DLAMCH
C     REAL SLAMCH
      INTEGER IMACH(16)
      SAVE FIRST, IMACH
      DATA FIRST / .TRUE. /
C
      IF( FIRST ) THEN
         FIRST = .FALSE.
         DO 10 K=1,16
         IMACH(K) = 0
   10    CONTINUE
C        IMACH(1) = 5
C        IMACH(2) = 6
         IMACH(4) = 6
C        IMACH(10)= NINT(DLAMCH('B'))
C        IMACH(11)= NINT(SLAMCH('N'))
C        IMACH(12)= NINT(SLAMCH('M'))
C        IMACH(13)= NINT(SLAMCH('L'))
C        IMACH(14)= NINT(DLAMCH('N'))
         IMACH(15)= NINT(DLAMCH('M'))
C        IMACH(16)= NINT(DLAMCH('L'))
      ENDIF
      IF( I.LT.1 .OR. I.GT.16 ) THEN
         WRITE (6,500)
         STOP 'I1MACH'
      ELSE
         I1MACH = IMACH(I)
      ENDIF
      RETURN
  500 FORMAT ('1ERROR    1 IN I1MACH - I OUT OF BOUNDS')
      END
