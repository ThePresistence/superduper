C     ******************************************************************
C
C     NMR 3-center integrals over Slater AOs: Low-level routines.
C
C     ******************************************************************
C
C
C     Modified Bessel function expansion of Slater AOs with respect
C     to displaced center.
C
C     Evaluation of Loewdin/Sharma expansion using expression of
C     A. Bouferguene, M. Fares, and P.E. Hoggan (Intl. J. Quantum Chem.
C     57, 801-810 (1996))
C
      SUBROUTINE PS3BI(LAMIN,LAMAX,IPMAX,ZETA,NPOINT,ROLESS,ROMORE,
     .                 AINT,LDA1,LDA2)
C
C   Compute batch of the inner integrals in BFH expansion.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      LAMIN  - Minimum lambda parameter value.
C      LAMAX  - Maximum lambda parameter value.
C      IPMAX  - Maximum p parameter value (minimum is zero).
C      ZETA   - Exponential parameter.
C      NPOINT - Number of ro points.
C      ROLESS - Argument values, at least NPOINT entries.
C               Grouping same ro< values together will improve
C               performance.
C      ROMORE - Argument values, at least NPOINT entries.
C               Grouping same ro> values together will improve
C               performance.
C      AINT   - (Output) integrals, three-index array.
C               First index:  LAMIN to LAMAX, lambda parameter values.
C               Second index: 0 to IPMAX, p parameter values.
C               Third index:  1 to NPOINT, ro argumnet values.
C      LDA1   - First leading dimension of AINT.
C      LDA2   - Second leading dimension of AINT.
C
C   Accessed common blocks:
C
C      PSPRT  - Printing unit.
C
C   Local storage:
C
C      MAXORD*(6+4*MAXLAS) + 2*MAXLAS*(MAXPS+1) DOUBLE PRECISION words.
C      With MAXORD=120, MAXLAS=4, and MAXPS=7, it is some 2700 words.
C      Only MAXORD*(2+MAXLAS) + MAXLAS*(MAXPS+1) (some 750 words) are 
C      accessed in the innermost loops, though.
C
C   Module logic:
C
C      The set of integrals:
C          1                      inf
C          /        la  x ro< zeta /        la  -y ro> zeta                p
C          | (1-x^2)   e           | (y^2-1)   e            (x ro< - y ro>)  dy dx
C          /                       /
C         -1                       1
C      is computed, using two-range numerical integration.
C      Integral in y is taken exactly, using the 2*LAMAX+IPMAX+1 Laguerre
C      expression. Integral in x is taken with 2*LAMAX+IPMAX+10 Gauss
C      formula.
C
C   Precision:
C
C      Precision is controlled by the orders of the inner Laguerre (ILO) and
C      outer Gaussian (IGO) integrations. Inner integral can be taken exactly
C      by setting ILO=2*LAMAX+IPMAX+1. Outer integral cannot be taken precisely
C      by Gaussian integration (except for trivial case zeta .eq. 0). Selection
C      IGO=2*LAMAX+IPMAX+4+MAX(ROLESS(I)*ZETA) should provide about nine to
C      ten significant digits. (Slightly lower for large ro< and zeta values).
C
C      This accuracy is excessive for high-order terms, and is reduced gradually
C      (see below), so what the absolute error is maintained about 10**(-12).
C
C   Possible optimizations:
C
C      Writing out more special-case inner loops might improve performance
C      a little bit.
C
C      On vector architectures, moving (a relatively long) IGO loop inside
C      (short) NLAS and IPMAX loops might help. This is a non-trivial amount
C      of work, though.
C
C   Bugs:
C
C      Number of lambda parameters which can be handled at the same time
C      is limited to MAXLAS. IPMAX can't exceed MAXPS. Integratiopn orders 
C      are limited to MAXORD.
C
C      This code is not particularly friendly to vector machines.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXLAS=5)
      PARAMETER (MAXPS=7)
      PARAMETER (MAXORD=120)
      PARAMETER (ZERO= 0.D0)
      PARAMETER (ONE = 1.D0)
      PARAMETER (SONE=-1.D0)
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
C
      DIMENSION ROLESS(1:NPOINT), ROMORE(1:NPOINT)
      DIMENSION AINT(LAMIN:LAMIN+LDA1-1,0:LDA2-1,1:NPOINT)
C    Local temporaries
      DIMENSION ARGG(MAXORD), WGTG(MAXORD), WGTGLA(0:MAXLAS-1,MAXORD)
      DIMENSION ARGL(MAXORD), WGTL(MAXORD), WGTLLA(0:MAXLAS-1,MAXORD)
      DIMENSION WGTGLE(0:MAXLAS-1,MAXORD)
      DIMENSION ARXROL(MAXORD), ARYROM(MAXORD)
      DIMENSION ACCUM(0:MAXLAS-1,0:MAXPS), AINNER(0:MAXLAS-1,0:MAXPS)
C    Range checking
      IF( LDA1.LT.LAMAX-LAMIN+1 ) THEN
          WRITE(NB6,11010) LDA1, LAMAX-LAMIN+1
          STOP 'PS3BI'
      ENDIF
      IF( LDA2.LT.IPMAX+1 ) THEN
          WRITE(NB6,11020) LDA2, IPMAX+1
          STOP 'PS3BI'
      ENDIF
      IF( LAMAX-LAMIN+1.GT.MAXLAS ) THEN
          WRITE(NB6,11030) LAMAX-LAMIN+1, MAXLAS
          STOP 'PS3BI'
      ENDIF
      IF( IPMAX.GT.MAXPS ) THEN
          WRITE(NB6,11050) IPMAX, MAXPS
          STOP 'PS3BI'
      ENDIF
      NLAS = LAMAX - LAMIN
C
C    Decide on the integration orders, and precompute (1-x^2)^la
C    factors at gaussian integration abscissas. Laguerre abscissas
C    depend on the ro> value, and will have to be evaluated later.
C
      IGOINC = 0
      DO 200 I=1,NPOINT
          IGOINC = MAX(IGOINC,NINT(ABS(ROLESS(I)*ZETA)))
  200 CONTINUE
      IGOINC = IGOINC + 4
      ILO = 2*LAMAX+IPMAX+1
      IGO = 2*LAMAX+IPMAX+IGOINC
      IF( ILO.GT.20 ) ILO = 20 + (ILO-20)/3
      IF( IGO.GT.25 ) IGO = 25 + (IGO-25)/3
C
C    Precompute (1-x^2)^la factors at gaussian integration abscissas. 
C    Laguerre abscissas depend on the ro> value, and will have to be 
C    evaluated later.
C
      CALL PSIGA(IGO,MAXORD,SONE,ONE,ARGG,WGTG)
      DO 1000 I=1,IGO
          TERM        = ONE - ARGG(I)**2
          WGTGLA(0,I) = WGTG(I)*TERM**LAMIN
          DO 900 J=1,NLAS
              WGTGLA(J,I) = WGTGLA(J-1,I)*TERM
  900     CONTINUE
 1000 CONTINUE
C    At this point, WGTGLA(J,I) contains a product of gaussian integration
C    weight at the point I by the (1-argg(i)^2)^(lamin+j)
C
      OLDROL = SONE
      OLDROM = SONE
C
C    Go over radial argument values
C
      DO 7000 IRO=1,NPOINT
          ROL = ROLESS(IRO)
          ROM = ROMORE(IRO)
C
C        If ro< changed from the previous point, recompute exponential
C        factors for outer integral.
C
          IF( ROL.NE.OLDROL ) THEN
              OLDROL = ROL
              DO 2000 I=1,IGO
                  ARXROL(I) = ARGG(I)*ROL
                  TERM = EXP(ARGG(I)*ROL*ZETA)
                  DO 1900 J=0,NLAS
                      WGTGLE(J,I) = WGTGLA(J,I)*TERM
 1900             CONTINUE
 2000         CONTINUE
          ENDIF
C
C        If ro> changed from the previous point, recompute Laguerre
C        abscissas and weights, as well as corresponding constant factors.
C
          IF( ROM.NE.OLDROM ) THEN
              OLDROM = ROM
              CALL PSILGB(ILO,MAXORD,ROM*ZETA,ONE,ARGL,WGTL)
              DO 3000 I=1,ILO
                  ARYROM(I)   = ARGL(I)*ROM
*                 WGTL(I)     = WGTL(I)*EXP(-ARGL(I)*ROM*ZETA)
                  TERM        = ARGL(I)**2 - ONE
                  WGTLLA(0,I) = WGTL(I)*TERM**LAMIN
                  DO 2900 J=1,NLAS
                      WGTLLA(J,I) = WGTLLA(J-1,I)*TERM
 2900             CONTINUE
 3000         CONTINUE
          ENDIF
C
C        Integrals for the given radial parameters ro are collected
C        in ACCUM.
C
          DO 3100 IP=0,IPMAX
              DO 3090 LA=0,NLAS
                  ACCUM(LA,IP) = ZERO
 3090         CONTINUE
 3100     CONTINUE
C
C        Actual numerical integration starts here.
C        See common-case loop (far) below at label 4991 
C        for comments.
C
C        NLAS and IPMAX usually take small integer values.
C        The horrow below is intended to eliminate loop
C        overhead in these common cases.
C
          ICASE = MIN(NLAS,4)*6 + MIN(IPMAX,5) + 1
          GOTO (4001,4011,4021,4031,4041,4091,
     .         4001,4111,4121,4131,4141,4191,
     .         4001,4211,4221,4231,4241,4291,
     .         4001,4311,4321,4331,4341,4391,
     .         4001,4911,4921,4931,4941,4991), ICASE
 4001     CONTINUE
          WRITE(NB6,11060) NLAS, IPMAX
          STOP 'PS3BI'
C
C            NLAS = 0, IPMAX = 1
C
 4011     CONTINUE
C
          IF( NLAS.NE.0 .OR. IPMAX.NE.1 ) STOP 'PS3BI'
              DO 4019 IX=1,IGO
                  A00 = ZERO
                  A01 = ZERO
                  DO 4015 IY=1,ILO
                      TOT  = WGTLLA(0,IY)
                      TERM = ARXROL(IX) - ARYROM(IY)
                      A00  = A00 + TOT
                      A01  = A01 + TOT*TERM
 4015             CONTINUE
                  W0         = WGTGLE(0,IX)
                  ACCUM(0,0) = ACCUM(0,0) + W0*A00
                  ACCUM(0,1) = ACCUM(0,1) + W0*A01
 4019         CONTINUE
              GOTO 5400
C
C            NLAS = 0, IPMAX = 2
C
 4021     CONTINUE
          IF( NLAS.NE.0 .OR. IPMAX.NE.2 ) STOP 'PS3BI'
              DO 4029 IX=1,IGO
                  A00 = ZERO
                  A01 = ZERO
                  A02 = ZERO
                  DO 4025 IY=1,ILO
                      TERM = ARXROL(IX) - ARYROM(IY)
                      TOT  = WGTLLA(0,IY)
                      A00  = A00 + TOT
                      A01  = A01 + TOT*TERM
                      A02  = A02 + TOT*TERM*TERM
 4025             CONTINUE
                  W0         = WGTGLE(0,IX)
                  ACCUM(0,0) = ACCUM(0,0) + W0*A00
                  ACCUM(0,1) = ACCUM(0,1) + W0*A01
                  ACCUM(0,2) = ACCUM(0,2) + W0*A02
 4029         CONTINUE
              GOTO 5400
C
C            NLAS = 0, IPMAX = 3
C
 4031     CONTINUE
          IF( NLAS.NE.0 .OR. IPMAX.NE.3 ) STOP 'PS3BI'
              DO 4039 IX=1,IGO
                  A00 = ZERO
                  A01 = ZERO
                  A02 = ZERO
                  A03 = ZERO
                  DO 4035 IY=1,ILO
                      TERM  = ARXROL(IX) - ARYROM(IY)
                      TOT   = WGTLLA(0,IY)
                      TERM2 = TERM*TERM
                      A00   = A00 + TOT
                      A01   = A01 + TOT*TERM
                      A02   = A02 + TOT*TERM2
                      A03   = A03 + TOT*TERM2*TERM
 4035             CONTINUE
                  W0         = WGTGLE(0,IX)
                  ACCUM(0,0) = ACCUM(0,0) + W0*A00
                  ACCUM(0,1) = ACCUM(0,1) + W0*A01
                  ACCUM(0,2) = ACCUM(0,2) + W0*A02
                  ACCUM(0,3) = ACCUM(0,3) + W0*A03
 4039         CONTINUE
              GOTO 5400
C
C            NLAS = 0, IPMAX = 4
C
 4041     CONTINUE
          IF( NLAS.NE.0 .OR. IPMAX.NE.4 ) STOP 'PS3BI'
              DO 4049 IX=1,IGO
                  A00 = ZERO
                  A01 = ZERO
                  A02 = ZERO
                  A03 = ZERO
                  A04 = ZERO
                  DO 4045 IY=1,ILO
                      TERM  = ARXROL(IX) - ARYROM(IY)
                      TOT   = WGTLLA(0,IY)
                      TERM2 = TERM*TERM
                      A00   = A00 + TOT
                      A01   = A01 + TOT*TERM
                      A02   = A02 + TOT*TERM2
                      A03   = A03 + TOT*TERM2*TERM
                      A04   = A04 + TOT*TERM2*TERM2
 4045             CONTINUE
                  W0         = WGTGLE(0,IX)
                  ACCUM(0,0) = ACCUM(0,0) + W0*A00
                  ACCUM(0,1) = ACCUM(0,1) + W0*A01
                  ACCUM(0,2) = ACCUM(0,2) + W0*A02
                  ACCUM(0,3) = ACCUM(0,3) + W0*A03
                  ACCUM(0,4) = ACCUM(0,4) + W0*A04
 4049         CONTINUE
              GOTO 5400
C
C            NLAS = 0
C
 4091     CONTINUE
          IF( NLAS.NE.0 ) STOP 'PS3BI'
              DO 4099 IX=1,IGO
                  DO 4093 IP=0,IPMAX
                      AINNER(0,IP) = ZERO
 4093             CONTINUE
                  DO 4095 IY=1,ILO
                      TERM = ARXROL(IX) - ARYROM(IY)
                      TOT  = WGTLLA(0,IY)
                      AINNER(0,0) = AINNER(0,0) + TOT
                      DO 4094 IP=1,IPMAX
                          TOT          = TOT * TERM
                          AINNER(0,IP) = AINNER(0,IP) + TOT
 4094                 CONTINUE 
 4095             CONTINUE
                  W0         = WGTGLE(0,IX)
                  DO 4098 IP=0,IPMAX
                      ACCUM(0,IP) = ACCUM(0,IP) + W0*AINNER(0,IP)
 4098             CONTINUE
 4099         CONTINUE
              GOTO 5400
C
C            NLAS = 1, IPMAX = 1
C
 4111     CONTINUE
          IF( NLAS.NE.1 .OR. IPMAX.NE.1 ) STOP 'PS3BI'
              DO 4119 IX=1,IGO
                  A00 = ZERO
                  A10 = ZERO
                  A01 = ZERO
                  A11 = ZERO
                  DO 4116 IY=1,ILO
                      TERM = ARXROL(IX) - ARYROM(IY)
                      TOT0 = WGTLLA(0,IY)
                      TOT1 = WGTLLA(1,IY)
                      A00  = A00 + TOT0
                      A10  = A10 + TOT1
                      A01  = A01 + TOT0*TERM
                      A11  = A11 + TOT1*TERM
 4116             CONTINUE
                  W0         = WGTGLE(0,IX)
                  W1         = WGTGLE(1,IX)
                  ACCUM(0,0) = ACCUM(0,0) + W0*A00
                  ACCUM(1,0) = ACCUM(1,0) + W1*A10
                  ACCUM(0,1) = ACCUM(0,1) + W0*A01
                  ACCUM(1,1) = ACCUM(1,1) + W1*A11
 4119         CONTINUE
              GOTO 5400
C
C            NLAS = 1, IPMAX = 2
C
 4121     CONTINUE
          IF( NLAS.NE.1 .OR. IPMAX.NE.2 ) STOP 'PS3BI'
              DO 4129 IX=1,IGO
                  A00 = ZERO
                  A10 = ZERO
                  A01 = ZERO
                  A11 = ZERO
                  A02 = ZERO
                  A12 = ZERO
                  DO 4126 IY=1,ILO
                      TERM  = ARXROL(IX) - ARYROM(IY)
                      TOT0  = WGTLLA(0,IY)
                      TOT1  = WGTLLA(1,IY)
                      TERM2 = TERM*TERM
                      A00   = A00 + TOT0
                      A10   = A10 + TOT1
                      A01   = A01 + TOT0*TERM
                      A11   = A11 + TOT1*TERM
                      A02   = A02 + TOT0*TERM2
                      A12   = A12 + TOT1*TERM2
 4126             CONTINUE
                  W0         = WGTGLE(0,IX)
                  W1         = WGTGLE(1,IX)
                  ACCUM(0,0) = ACCUM(0,0) + W0*A00
                  ACCUM(1,0) = ACCUM(1,0) + W1*A10
                  ACCUM(0,1) = ACCUM(0,1) + W0*A01
                  ACCUM(1,1) = ACCUM(1,1) + W1*A11
                  ACCUM(0,2) = ACCUM(0,2) + W0*A02
                  ACCUM(1,2) = ACCUM(1,2) + W1*A12
 4129         CONTINUE
              GOTO 5400
C
C            NLAS = 1, IPMAX = 3
C
 4131     CONTINUE
          IF( NLAS.NE.1 .OR. IPMAX.NE.3 ) STOP 'PS3BI'
              DO 4139 IX=1,IGO
                  A00 = ZERO
                  A10 = ZERO
                  A01 = ZERO
                  A11 = ZERO
                  A02 = ZERO
                  A12 = ZERO
                  A03 = ZERO
                  A13 = ZERO
                  DO 4136 IY=1,ILO
                      TERM  = ARXROL(IX) - ARYROM(IY)
                      TOT0  = WGTLLA(0,IY)
                      TOT1  = WGTLLA(1,IY)
                      TERM2 = TERM*TERM
                      A00   = A00 + TOT0
                      A10   = A10 + TOT1
                      TERM3 = TERM*TERM2
                      A01   = A01 + TOT0*TERM
                      A11   = A11 + TOT1*TERM
                      A02   = A02 + TOT0*TERM2
                      A12   = A12 + TOT1*TERM2
                      A03   = A03 + TOT0*TERM3
                      A13   = A13 + TOT1*TERM3
 4136             CONTINUE
                  W0         = WGTGLE(0,IX)
                  W1         = WGTGLE(1,IX)
                  ACCUM(0,0) = ACCUM(0,0) + W0*A00
                  ACCUM(1,0) = ACCUM(1,0) + W1*A10
                  ACCUM(0,1) = ACCUM(0,1) + W0*A01
                  ACCUM(1,1) = ACCUM(1,1) + W1*A11
                  ACCUM(0,2) = ACCUM(0,2) + W0*A02
                  ACCUM(1,2) = ACCUM(1,2) + W1*A12
                  ACCUM(0,3) = ACCUM(0,3) + W0*A03
                  ACCUM(1,3) = ACCUM(1,3) + W1*A13
 4139         CONTINUE
              GOTO 5400
C
C            NLAS = 1, IPMAX = 4
C
 4141     CONTINUE
          IF( NLAS.NE.1 .OR. IPMAX.NE.4 ) STOP 'PS3BI'
              DO 4149 IX=1,IGO
                  A00 = ZERO
                  A10 = ZERO
                  A01 = ZERO
                  A11 = ZERO
                  A02 = ZERO
                  A12 = ZERO
                  A03 = ZERO
                  A13 = ZERO
                  A04 = ZERO
                  A14 = ZERO
                  DO 4146 IY=1,ILO
                      TERM  = ARXROL(IX) - ARYROM(IY)
                      TOT0  = WGTLLA(0,IY)
                      TOT1  = WGTLLA(1,IY)
                      TERM2 = TERM*TERM
                      A00   = A00 + TOT0
                      A10   = A10 + TOT1
                      TERM3 = TERM*TERM2
                      TERM4 = TERM2*TERM2
                      A01   = A01 + TOT0*TERM
                      A11   = A11 + TOT1*TERM
                      A02   = A02 + TOT0*TERM2
                      A12   = A12 + TOT1*TERM2
                      A03   = A03 + TOT0*TERM3
                      A13   = A13 + TOT1*TERM3
                      A04   = A04 + TOT0*TERM4
                      A14   = A14 + TOT1*TERM4
 4146             CONTINUE
                  W0         = WGTGLE(0,IX)
                  W1         = WGTGLE(1,IX)
                  ACCUM(0,0) = ACCUM(0,0) + W0*A00
                  ACCUM(1,0) = ACCUM(1,0) + W1*A10
                  ACCUM(0,1) = ACCUM(0,1) + W0*A01
                  ACCUM(1,1) = ACCUM(1,1) + W1*A11
                  ACCUM(0,2) = ACCUM(0,2) + W0*A02
                  ACCUM(1,2) = ACCUM(1,2) + W1*A12
                  ACCUM(0,3) = ACCUM(0,3) + W0*A03
                  ACCUM(1,3) = ACCUM(1,3) + W1*A13
                  ACCUM(0,4) = ACCUM(0,4) + W0*A04
                  ACCUM(1,4) = ACCUM(1,4) + W1*A14
 4149         CONTINUE
              GOTO 5400
C
C            NLAS = 1
C
 4191     CONTINUE
          IF( NLAS.NE.1 ) STOP 'PS3BI'
              DO 4199 IX=1,IGO
                  DO 4193 IP=0,IPMAX
                      AINNER(0,IP) = ZERO
                      AINNER(1,IP) = ZERO
 4193             CONTINUE
                  DO 4196 IY=1,ILO
                      TERM        = ARXROL(IX) - ARYROM(IY)
                      TOT0        = WGTLLA(0,IY)
                      TOT1        = WGTLLA(1,IY)
                      AINNER(0,0) = AINNER(0,0) + TOT0
                      AINNER(1,0) = AINNER(1,0) + TOT1
                      DO 4194 IP=1,IPMAX
                          TOT0         = TOT0 * TERM
                          TOT1         = TOT1 * TERM
                          AINNER(0,IP) = AINNER(0,IP) + TOT0
                          AINNER(1,IP) = AINNER(1,IP) + TOT1
 4194                 CONTINUE 
 4196             CONTINUE
                  W0         = WGTGLE(0,IX)
                  W1         = WGTGLE(1,IX)
                  DO 4198 IP=0,IPMAX
                      ACCUM(0,IP) = ACCUM(0,IP) + W0*AINNER(0,IP)
                      ACCUM(1,IP) = ACCUM(1,IP) + W1*AINNER(1,IP)
 4198             CONTINUE
 4199         CONTINUE
              GOTO 5400
C
C            NLAS = 2, IPMAX = 1
C
 4211     CONTINUE
          IF( NLAS.NE.2 .OR. IPMAX.NE.1 ) STOP 'PS3BI'
              DO 4219 IX=1,IGO
                  A00 = ZERO
                  A10 = ZERO
                  A20 = ZERO
                  A01 = ZERO
                  A11 = ZERO
                  A21 = ZERO
                  DO 4216 IY=1,ILO
                      TERM = ARXROL(IX) - ARYROM(IY)
                      TOT0 = WGTLLA(0,IY)
                      TOT1 = WGTLLA(1,IY)
                      TOT2 = WGTLLA(2,IY)
                      A00  = A00 + TOT0
                      A10  = A10 + TOT1
                      A20  = A20 + TOT2
                      A01  = A01 + TOT0*TERM
                      A11  = A11 + TOT1*TERM
                      A21  = A21 + TOT2*TERM
 4216             CONTINUE
                  W0         = WGTGLE(0,IX)
                  W1         = WGTGLE(1,IX)
                  W2         = WGTGLE(2,IX)
                  ACCUM(0,0) = ACCUM(0,0) + W0*A00
                  ACCUM(1,0) = ACCUM(1,0) + W1*A10
                  ACCUM(2,0) = ACCUM(2,0) + W2*A20
                  ACCUM(0,1) = ACCUM(0,1) + W0*A01
                  ACCUM(1,1) = ACCUM(1,1) + W1*A11
                  ACCUM(2,1) = ACCUM(2,1) + W2*A21
 4219         CONTINUE
              GOTO 5400
C
C            NLAS = 2, IPMAX = 2
C
 4221     CONTINUE
          IF( NLAS.NE.2 .OR. IPMAX.NE.2 ) STOP 'PS3BI'
              DO 4229 IX=1,IGO
                  A00 = ZERO
                  A10 = ZERO
                  A20 = ZERO
                  A01 = ZERO
                  A11 = ZERO
                  A21 = ZERO
                  A02 = ZERO
                  A12 = ZERO
                  A22 = ZERO
                  DO 4226 IY=1,ILO
                      TERM  = ARXROL(IX) - ARYROM(IY)
                      TOT0  = WGTLLA(0,IY)
                      TOT1  = WGTLLA(1,IY)
                      TOT2  = WGTLLA(2,IY)
                      TERM2 = TERM*TERM
                      A00   = A00 + TOT0
                      A10   = A10 + TOT1
                      A20   = A20 + TOT2
                      A01   = A01 + TOT0*TERM
                      A11   = A11 + TOT1*TERM
                      A21   = A21 + TOT2*TERM
                      A02   = A02 + TOT0*TERM2
                      A12   = A12 + TOT1*TERM2
                      A22   = A22 + TOT2*TERM2
 4226             CONTINUE
                  W0         = WGTGLE(0,IX)
                  W1         = WGTGLE(1,IX)
                  W2         = WGTGLE(2,IX)
                  ACCUM(0,0) = ACCUM(0,0) + W0*A00
                  ACCUM(1,0) = ACCUM(1,0) + W1*A10
                  ACCUM(2,0) = ACCUM(2,0) + W2*A20
                  ACCUM(0,1) = ACCUM(0,1) + W0*A01
                  ACCUM(1,1) = ACCUM(1,1) + W1*A11
                  ACCUM(2,1) = ACCUM(2,1) + W2*A21
                  ACCUM(0,2) = ACCUM(0,2) + W0*A02
                  ACCUM(1,2) = ACCUM(1,2) + W1*A12
                  ACCUM(2,2) = ACCUM(2,2) + W2*A22
 4229         CONTINUE
              GOTO 5400
C
C            NLAS = 2, IPMAX = 3
C
 4231     CONTINUE
          IF( NLAS.NE.2 .OR. IPMAX.NE.3 ) STOP 'PS3BI'
              DO 4239 IX=1,IGO
                  A00 = ZERO
                  A10 = ZERO
                  A20 = ZERO
                  A01 = ZERO
                  A11 = ZERO
                  A21 = ZERO
                  A02 = ZERO
                  A12 = ZERO
                  A22 = ZERO
                  A03 = ZERO
                  A13 = ZERO
                  A23 = ZERO
                  DO 4236 IY=1,ILO
                      TERM  = ARXROL(IX) - ARYROM(IY)
                      TOT0  = WGTLLA(0,IY)
                      TOT1  = WGTLLA(1,IY)
                      TOT2  = WGTLLA(2,IY)
                      TERM2 = TERM*TERM
                      A00   = A00 + TOT0
                      A10   = A10 + TOT1
                      A20   = A20 + TOT2
                      TERM3 = TERM2*TERM
                      A01   = A01 + TOT0*TERM
                      A11   = A11 + TOT1*TERM
                      A21   = A21 + TOT2*TERM
                      A02   = A02 + TOT0*TERM2
                      A12   = A12 + TOT1*TERM2
                      A22   = A22 + TOT2*TERM2
                      A03   = A03 + TOT0*TERM3
                      A13   = A13 + TOT1*TERM3
                      A23   = A23 + TOT2*TERM3
 4236             CONTINUE
                  W0         = WGTGLE(0,IX)
                  W1         = WGTGLE(1,IX)
                  W2         = WGTGLE(2,IX)
                  ACCUM(0,0) = ACCUM(0,0) + W0*A00
                  ACCUM(1,0) = ACCUM(1,0) + W1*A10
                  ACCUM(2,0) = ACCUM(2,0) + W2*A20
                  ACCUM(0,1) = ACCUM(0,1) + W0*A01
                  ACCUM(1,1) = ACCUM(1,1) + W1*A11
                  ACCUM(2,1) = ACCUM(2,1) + W2*A21
                  ACCUM(0,2) = ACCUM(0,2) + W0*A02
                  ACCUM(1,2) = ACCUM(1,2) + W1*A12
                  ACCUM(2,2) = ACCUM(2,2) + W2*A22
                  ACCUM(0,3) = ACCUM(0,3) + W0*A03
                  ACCUM(1,3) = ACCUM(1,3) + W1*A13
                  ACCUM(2,3) = ACCUM(2,3) + W2*A23
 4239         CONTINUE
              GOTO 5400
C
C            NLAS = 2, IPMAX = 4
C
 4241     CONTINUE
          IF( NLAS.NE.2 .OR. IPMAX.NE.4 ) STOP 'PS3BI'
              DO 4249 IX=1,IGO
                  A00 = ZERO
                  A10 = ZERO
                  A20 = ZERO
                  A01 = ZERO
                  A11 = ZERO
                  A21 = ZERO
                  A02 = ZERO
                  A12 = ZERO
                  A22 = ZERO
                  A03 = ZERO
                  A13 = ZERO
                  A23 = ZERO
                  A04 = ZERO
                  A14 = ZERO
                  A24 = ZERO
                  DO 4246 IY=1,ILO
                      TERM  = ARXROL(IX) - ARYROM(IY)
                      TOT0  = WGTLLA(0,IY)
                      TOT1  = WGTLLA(1,IY)
                      TOT2  = WGTLLA(2,IY)
                      TERM2 = TERM*TERM
                      A00   = A00 + TOT0
                      A10   = A10 + TOT1
                      A20   = A20 + TOT2
                      TERM3 = TERM2*TERM
                      TERM4 = TERM2*TERM2
                      A01   = A01 + TOT0*TERM
                      A11   = A11 + TOT1*TERM
                      A21   = A21 + TOT2*TERM
                      A02   = A02 + TOT0*TERM2
                      A12   = A12 + TOT1*TERM2
                      A22   = A22 + TOT2*TERM2
                      A03   = A03 + TOT0*TERM3
                      A13   = A13 + TOT1*TERM3
                      A23   = A23 + TOT2*TERM3
                      A04   = A04 + TOT0*TERM4
                      A14   = A14 + TOT1*TERM4
                      A24   = A24 + TOT2*TERM4
 4246             CONTINUE
                  W0         = WGTGLE(0,IX)
                  W1         = WGTGLE(1,IX)
                  W2         = WGTGLE(2,IX)
                  ACCUM(0,0) = ACCUM(0,0) + W0*A00
                  ACCUM(1,0) = ACCUM(1,0) + W1*A10
                  ACCUM(2,0) = ACCUM(2,0) + W2*A20
                  ACCUM(0,1) = ACCUM(0,1) + W0*A01
                  ACCUM(1,1) = ACCUM(1,1) + W1*A11
                  ACCUM(2,1) = ACCUM(2,1) + W2*A21
                  ACCUM(0,2) = ACCUM(0,2) + W0*A02
                  ACCUM(1,2) = ACCUM(1,2) + W1*A12
                  ACCUM(2,2) = ACCUM(2,2) + W2*A22
                  ACCUM(0,3) = ACCUM(0,3) + W0*A03
                  ACCUM(1,3) = ACCUM(1,3) + W1*A13
                  ACCUM(2,3) = ACCUM(2,3) + W2*A23
                  ACCUM(0,4) = ACCUM(0,4) + W0*A04
                  ACCUM(1,4) = ACCUM(1,4) + W1*A14
                  ACCUM(2,4) = ACCUM(2,4) + W2*A24
 4249         CONTINUE
              GOTO 5400
C
C            NLAS = 2
C
 4291     CONTINUE
          IF( NLAS.NE.2 ) STOP 'PS3BI'
              DO 4299 IX=1,IGO
                  DO 4293 IP=0,IPMAX
                      AINNER(0,IP) = ZERO
                      AINNER(1,IP) = ZERO
                      AINNER(2,IP) = ZERO
 4293             CONTINUE
                  DO 4296 IY=1,ILO
                      TERM        = ARXROL(IX) - ARYROM(IY)
                      TOT0        = WGTLLA(0,IY)
                      TOT1        = WGTLLA(1,IY)
                      TOT2        = WGTLLA(2,IY)
                      AINNER(0,0) = AINNER(0,0) + TOT0
                      AINNER(1,0) = AINNER(1,0) + TOT1
                      AINNER(2,0) = AINNER(2,0) + TOT2
                      DO 4294 IP=1,IPMAX
                          TOT0         = TOT0 * TERM
                          TOT1         = TOT1 * TERM
                          TOT2         = TOT2 * TERM
                          AINNER(0,IP) = AINNER(0,IP) + TOT0
                          AINNER(1,IP) = AINNER(1,IP) + TOT1
                          AINNER(2,IP) = AINNER(2,IP) + TOT2
 4294                 CONTINUE 
 4296             CONTINUE
                  W0         = WGTGLE(0,IX)
                  W1         = WGTGLE(1,IX)
                  W2         = WGTGLE(2,IX)
                  DO 4298 IP=0,IPMAX
                      ACCUM(0,IP) = ACCUM(0,IP) + W0*AINNER(0,IP)
                      ACCUM(1,IP) = ACCUM(1,IP) + W1*AINNER(1,IP)
                      ACCUM(2,IP) = ACCUM(2,IP) + W2*AINNER(2,IP)
 4298             CONTINUE
 4299         CONTINUE
              GOTO 5400
C
C            NLAS = 3, IPMAX = 1
C
 4311     CONTINUE
          IF( NLAS.NE.3 .OR. IPMAX.NE.1 ) STOP 'PS3BI'
              DO 4319 IX=1,IGO
                  A00 = ZERO
                  A10 = ZERO
                  A20 = ZERO
                  A30 = ZERO
                  A01 = ZERO
                  A11 = ZERO
                  A21 = ZERO
                  A31 = ZERO
                  DO 4316 IY=1,ILO
                      TERM = ARXROL(IX) - ARYROM(IY)
                      TOT0 = WGTLLA(0,IY)
                      TOT1 = WGTLLA(1,IY)
                      TOT2 = WGTLLA(2,IY)
                      TOT3 = WGTLLA(3,IY)
                      A00  = A00 + TOT0
                      A10  = A10 + TOT1
                      A20  = A20 + TOT2
                      A30  = A30 + TOT3
                      A01  = A01 + TOT0*TERM
                      A11  = A11 + TOT1*TERM
                      A21  = A21 + TOT2*TERM
                      A31  = A31 + TOT3*TERM
 4316             CONTINUE
                  W0         = WGTGLE(0,IX)
                  W1         = WGTGLE(1,IX)
                  W2         = WGTGLE(2,IX)
                  W3         = WGTGLE(3,IX)
                  ACCUM(0,0) = ACCUM(0,0) + W0*A00
                  ACCUM(1,0) = ACCUM(1,0) + W1*A10
                  ACCUM(2,0) = ACCUM(2,0) + W2*A20
                  ACCUM(3,0) = ACCUM(3,0) + W3*A30
                  ACCUM(0,1) = ACCUM(0,1) + W0*A01
                  ACCUM(1,1) = ACCUM(1,1) + W1*A11
                  ACCUM(2,1) = ACCUM(2,1) + W2*A21
                  ACCUM(3,1) = ACCUM(3,1) + W3*A31
 4319         CONTINUE
              GOTO 5400
C
C            NLAS = 3, IPMAX = 2
C
 4321     CONTINUE
          IF( NLAS.NE.3 .OR. IPMAX.NE.2 ) STOP 'PS3BI'
              DO 4329 IX=1,IGO
                  A00 = ZERO
                  A10 = ZERO
                  A20 = ZERO
                  A30 = ZERO
                  A01 = ZERO
                  A11 = ZERO
                  A21 = ZERO
                  A31 = ZERO
                  A02 = ZERO
                  A12 = ZERO
                  A22 = ZERO
                  A32 = ZERO
                  DO 4326 IY=1,ILO
                      TERM  = ARXROL(IX) - ARYROM(IY)
                      TOT0  = WGTLLA(0,IY)
                      TOT1  = WGTLLA(1,IY)
                      TOT2  = WGTLLA(2,IY)
                      TOT3  = WGTLLA(3,IY)
                      TERM2 = TERM*TERM
                      A00   = A00 + TOT0
                      A10   = A10 + TOT1
                      A20   = A20 + TOT2
                      A30   = A30 + TOT3
                      A01   = A01 + TOT0*TERM
                      A11   = A11 + TOT1*TERM
                      A21   = A21 + TOT2*TERM
                      A31   = A31 + TOT3*TERM
                      A02   = A02 + TOT0*TERM2
                      A12   = A12 + TOT1*TERM2
                      A22   = A22 + TOT2*TERM2
                      A32   = A32 + TOT3*TERM2
 4326             CONTINUE
                  W0         = WGTGLE(0,IX)
                  W1         = WGTGLE(1,IX)
                  W2         = WGTGLE(2,IX)
                  W3         = WGTGLE(3,IX)
                  ACCUM(0,0) = ACCUM(0,0) + W0*A00
                  ACCUM(1,0) = ACCUM(1,0) + W1*A10
                  ACCUM(2,0) = ACCUM(2,0) + W2*A20
                  ACCUM(3,0) = ACCUM(3,0) + W3*A30
                  ACCUM(0,1) = ACCUM(0,1) + W0*A01
                  ACCUM(1,1) = ACCUM(1,1) + W1*A11
                  ACCUM(2,1) = ACCUM(2,1) + W2*A21
                  ACCUM(3,1) = ACCUM(3,1) + W3*A31
                  ACCUM(0,2) = ACCUM(0,2) + W0*A02
                  ACCUM(1,2) = ACCUM(1,2) + W1*A12
                  ACCUM(2,2) = ACCUM(2,2) + W2*A22
                  ACCUM(3,2) = ACCUM(3,2) + W3*A32
 4329         CONTINUE
              GOTO 5400
C
C            NLAS = 3, IPMAX = 3
C
 4331     CONTINUE
          IF( NLAS.NE.3 .OR. IPMAX.NE.3 ) STOP 'PS3BI'
              DO 4339 IX=1,IGO
                  A00 = ZERO
                  A10 = ZERO
                  A20 = ZERO
                  A30 = ZERO
                  A01 = ZERO
                  A11 = ZERO
                  A21 = ZERO
                  A31 = ZERO
                  A02 = ZERO
                  A12 = ZERO
                  A22 = ZERO
                  A32 = ZERO
                  A03 = ZERO
                  A13 = ZERO
                  A23 = ZERO
                  A33 = ZERO
                  DO 4336 IY=1,ILO
                      TERM  = ARXROL(IX) - ARYROM(IY)
                      TOT0  = WGTLLA(0,IY)
                      TOT1  = WGTLLA(1,IY)
                      TOT2  = WGTLLA(2,IY)
                      TOT3  = WGTLLA(3,IY)
                      TERM2 = TERM*TERM
                      A00   = A00 + TOT0
                      A10   = A10 + TOT1
                      A20   = A20 + TOT2
                      A30   = A30 + TOT3
                      TERM3 = TERM2*TERM
                      A01   = A01 + TOT0*TERM
                      A11   = A11 + TOT1*TERM
                      A21   = A21 + TOT2*TERM
                      A31   = A31 + TOT3*TERM
                      A02   = A02 + TOT0*TERM2
                      A12   = A12 + TOT1*TERM2
                      A22   = A22 + TOT2*TERM2
                      A32   = A32 + TOT3*TERM2
                      A03   = A03 + TOT0*TERM3
                      A13   = A13 + TOT1*TERM3
                      A23   = A23 + TOT2*TERM3
                      A33   = A33 + TOT3*TERM3
 4336             CONTINUE
                  W0         = WGTGLE(0,IX)
                  W1         = WGTGLE(1,IX)
                  W2         = WGTGLE(2,IX)
                  W3         = WGTGLE(3,IX)
                  ACCUM(0,0) = ACCUM(0,0) + W0*A00
                  ACCUM(1,0) = ACCUM(1,0) + W1*A10
                  ACCUM(2,0) = ACCUM(2,0) + W2*A20
                  ACCUM(3,0) = ACCUM(3,0) + W3*A30
                  ACCUM(0,1) = ACCUM(0,1) + W0*A01
                  ACCUM(1,1) = ACCUM(1,1) + W1*A11
                  ACCUM(2,1) = ACCUM(2,1) + W2*A21
                  ACCUM(3,1) = ACCUM(3,1) + W3*A31
                  ACCUM(0,2) = ACCUM(0,2) + W0*A02
                  ACCUM(1,2) = ACCUM(1,2) + W1*A12
                  ACCUM(2,2) = ACCUM(2,2) + W2*A22
                  ACCUM(3,2) = ACCUM(3,2) + W3*A32
                  ACCUM(0,3) = ACCUM(0,3) + W0*A03
                  ACCUM(1,3) = ACCUM(1,3) + W1*A13
                  ACCUM(2,3) = ACCUM(2,3) + W2*A23
                  ACCUM(3,3) = ACCUM(3,3) + W3*A33
 4339         CONTINUE
              GOTO 5400
C
C            NLAS = 3, IPMAX = 4
C
 4341     CONTINUE
          IF( NLAS.NE.3 .OR. IPMAX.NE.4 ) STOP 'PS3BI'
              DO 4349 IX=1,IGO
                  A00 = ZERO
                  A10 = ZERO
                  A20 = ZERO
                  A30 = ZERO
                  A01 = ZERO
                  A11 = ZERO
                  A21 = ZERO
                  A31 = ZERO
                  A02 = ZERO
                  A12 = ZERO
                  A22 = ZERO
                  A32 = ZERO
                  A03 = ZERO
                  A13 = ZERO
                  A23 = ZERO
                  A33 = ZERO
                  A04 = ZERO
                  A14 = ZERO
                  A24 = ZERO
                  A34 = ZERO
                  DO 4346 IY=1,ILO
                      TERM  = ARXROL(IX) - ARYROM(IY)
                      TOT0  = WGTLLA(0,IY)
                      TOT1  = WGTLLA(1,IY)
                      TOT2  = WGTLLA(2,IY)
                      TOT3  = WGTLLA(3,IY)
                      TERM2 = TERM*TERM
                      A00   = A00 + TOT0
                      A10   = A10 + TOT1
                      A20   = A20 + TOT2
                      A30   = A30 + TOT3
                      TERM3 = TERM2*TERM
                      TERM4 = TERM2*TERM2
                      A01   = A01 + TOT0*TERM
                      A11   = A11 + TOT1*TERM
                      A21   = A21 + TOT2*TERM
                      A31   = A31 + TOT3*TERM
                      A02   = A02 + TOT0*TERM2
                      A12   = A12 + TOT1*TERM2
                      A22   = A22 + TOT2*TERM2
                      A32   = A32 + TOT3*TERM2
                      A03   = A03 + TOT0*TERM3
                      A13   = A13 + TOT1*TERM3
                      A23   = A23 + TOT2*TERM3
                      A33   = A33 + TOT3*TERM3
                      A04   = A04 + TOT0*TERM4
                      A14   = A14 + TOT1*TERM4
                      A24   = A24 + TOT2*TERM4
                      A34   = A34 + TOT3*TERM4
 4346             CONTINUE
                  W0         = WGTGLE(0,IX)
                  W1         = WGTGLE(1,IX)
                  W2         = WGTGLE(2,IX)
                  W3         = WGTGLE(3,IX)
                  ACCUM(0,0) = ACCUM(0,0) + W0*A00
                  ACCUM(1,0) = ACCUM(1,0) + W1*A10
                  ACCUM(2,0) = ACCUM(2,0) + W2*A20
                  ACCUM(3,0) = ACCUM(3,0) + W3*A30
                  ACCUM(0,1) = ACCUM(0,1) + W0*A01
                  ACCUM(1,1) = ACCUM(1,1) + W1*A11
                  ACCUM(2,1) = ACCUM(2,1) + W2*A21
                  ACCUM(3,1) = ACCUM(3,1) + W3*A31
                  ACCUM(0,2) = ACCUM(0,2) + W0*A02
                  ACCUM(1,2) = ACCUM(1,2) + W1*A12
                  ACCUM(2,2) = ACCUM(2,2) + W2*A22
                  ACCUM(3,2) = ACCUM(3,2) + W3*A32
                  ACCUM(0,3) = ACCUM(0,3) + W0*A03
                  ACCUM(1,3) = ACCUM(1,3) + W1*A13
                  ACCUM(2,3) = ACCUM(2,3) + W2*A23
                  ACCUM(3,3) = ACCUM(3,3) + W3*A33
                  ACCUM(0,4) = ACCUM(0,4) + W0*A04
                  ACCUM(1,4) = ACCUM(1,4) + W1*A14
                  ACCUM(2,4) = ACCUM(2,4) + W2*A24
                  ACCUM(3,4) = ACCUM(3,4) + W3*A34
 4349         CONTINUE
              GOTO 5400
C
C            NLAS = 3
C
 4391     CONTINUE
          IF( NLAS.NE.3 ) STOP 'PS3BI'
              DO 4399 IX=1,IGO
                  DO 4393 IP=0,IPMAX
                      AINNER(0,IP) = ZERO
                      AINNER(1,IP) = ZERO
                      AINNER(2,IP) = ZERO
                      AINNER(3,IP) = ZERO
 4393             CONTINUE
                  DO 4396 IY=1,ILO
                      TERM        = ARXROL(IX) - ARYROM(IY)
                      TOT0        = WGTLLA(0,IY)
                      TOT1        = WGTLLA(1,IY)
                      TOT2        = WGTLLA(2,IY)
                      TOT3        = WGTLLA(3,IY)
                      AINNER(0,0) = AINNER(0,0) + TOT0
                      AINNER(1,0) = AINNER(1,0) + TOT1
                      AINNER(2,0) = AINNER(2,0) + TOT2
                      AINNER(3,0) = AINNER(3,0) + TOT3
                      DO 4394 IP=1,IPMAX
                          TOT0         = TOT0 * TERM
                          TOT1         = TOT1 * TERM
                          TOT2         = TOT2 * TERM
                          TOT3         = TOT3 * TERM
                          AINNER(0,IP) = AINNER(0,IP) + TOT0
                          AINNER(1,IP) = AINNER(1,IP) + TOT1
                          AINNER(2,IP) = AINNER(2,IP) + TOT2
                          AINNER(3,IP) = AINNER(3,IP) + TOT3
 4394                 CONTINUE 
 4396             CONTINUE
                  W0         = WGTGLE(0,IX)
                  W1         = WGTGLE(1,IX)
                  W2         = WGTGLE(2,IX)
                  W3         = WGTGLE(3,IX)
                  DO 4398 IP=0,IPMAX
                      ACCUM(0,IP) = ACCUM(0,IP) + W0*AINNER(0,IP)
                      ACCUM(1,IP) = ACCUM(1,IP) + W1*AINNER(1,IP)
                      ACCUM(2,IP) = ACCUM(2,IP) + W2*AINNER(2,IP)
                      ACCUM(3,IP) = ACCUM(3,IP) + W3*AINNER(3,IP)
 4398             CONTINUE
 4399         CONTINUE
              GOTO 5400
C
C            IPMAX = 1
C
 4911     CONTINUE
          IF( IPMAX.NE.1 ) STOP 'PS3BI'
              DO 4919 IX=1,IGO
                  DO 4912 LA=0,NLAS
                      AINNER(LA,0) = ZERO
                      AINNER(LA,1) = ZERO
 4912             CONTINUE
                  DO 4916 IY=1,ILO
                      TERM = ARXROL(IX) - ARYROM(IY)
                      DO 4915 LA=0,NLAS
                          TOT          = WGTLLA(LA,IY)
                          AINNER(LA,0) = AINNER(LA,0) + TOT
                          AINNER(LA,1) = AINNER(LA,1) + TOT*TERM
 4915                 CONTINUE
 4916             CONTINUE
                  DO 4917 LA=0,NLAS
                      W           = WGTGLE(LA,IX)
                      ACCUM(LA,0) = ACCUM(LA,0) + W*AINNER(LA,0)
                      ACCUM(LA,1) = ACCUM(LA,1) + W*AINNER(LA,1)
 4917             CONTINUE
 4919         CONTINUE
              GOTO 5400
C
C            IPMAX = 2
C
 4921     CONTINUE
          IF( IPMAX.NE.2 ) STOP 'PS3BI'
              DO 4929 IX=1,IGO
                  DO 4922 LA=0,NLAS
                      AINNER(LA,0) = ZERO
                      AINNER(LA,1) = ZERO
                      AINNER(LA,2) = ZERO
 4922             CONTINUE
                  DO 4926 IY=1,ILO
                      TERM  = ARXROL(IX) - ARYROM(IY)
                      TERM2 = TERM*TERM
                      DO 4925 LA=0,NLAS
                          TOT          = WGTLLA(LA,IY)
                          AINNER(LA,0) = AINNER(LA,0) + TOT
                          AINNER(LA,1) = AINNER(LA,1) + TOT*TERM
                          AINNER(LA,2) = AINNER(LA,2) + TOT*TERM2
 4925                 CONTINUE
 4926             CONTINUE
                  DO 4927 LA=0,NLAS
                      W           = WGTGLE(LA,IX)
                      ACCUM(LA,0) = ACCUM(LA,0) + W*AINNER(LA,0)
                      ACCUM(LA,1) = ACCUM(LA,1) + W*AINNER(LA,1)
                      ACCUM(LA,2) = ACCUM(LA,2) + W*AINNER(LA,2)
 4927             CONTINUE
 4929         CONTINUE
              GOTO 5400
C
C            IPMAX = 3
C
 4931     CONTINUE
          IF( IPMAX.NE.3 ) STOP 'PS3BI'
              DO 4939 IX=1,IGO
                  DO 4932 LA=0,NLAS
                      AINNER(LA,0) = ZERO
                      AINNER(LA,1) = ZERO
                      AINNER(LA,2) = ZERO
                      AINNER(LA,3) = ZERO
 4932             CONTINUE
                  DO 4936 IY=1,ILO
                      TERM  = ARXROL(IX) - ARYROM(IY)
                      TERM2 = TERM*TERM
                      TERM3 = TERM2*TERM
                      DO 4935 LA=0,NLAS
                          TOT          = WGTLLA(LA,IY)
                          AINNER(LA,0) = AINNER(LA,0) + TOT
                          AINNER(LA,1) = AINNER(LA,1) + TOT*TERM
                          AINNER(LA,2) = AINNER(LA,2) + TOT*TERM2
                          AINNER(LA,3) = AINNER(LA,3) + TOT*TERM3
 4935                 CONTINUE
 4936             CONTINUE
                  DO 4937 LA=0,NLAS
                      W           = WGTGLE(LA,IX)
                      ACCUM(LA,0) = ACCUM(LA,0) + W*AINNER(LA,0)
                      ACCUM(LA,1) = ACCUM(LA,1) + W*AINNER(LA,1)
                      ACCUM(LA,2) = ACCUM(LA,2) + W*AINNER(LA,2)
                      ACCUM(LA,3) = ACCUM(LA,3) + W*AINNER(LA,3)
 4937             CONTINUE
 4939         CONTINUE
              GOTO 5400
C
C            IPMAX = 4
C
 4941     CONTINUE
          IF( IPMAX.NE.4 ) STOP 'PS3BI'
              DO 4949 IX=1,IGO
                  DO 4942 LA=0,NLAS
                      AINNER(LA,0) = ZERO
                      AINNER(LA,1) = ZERO
                      AINNER(LA,2) = ZERO
                      AINNER(LA,3) = ZERO
                      AINNER(LA,4) = ZERO
 4942             CONTINUE
                  DO 4946 IY=1,ILO
                      TERM  = ARXROL(IX) - ARYROM(IY)
                      TERM2 = TERM*TERM
                      TERM3 = TERM2*TERM
                      TERM4 = TERM2*TERM2
                      DO 4945 LA=0,NLAS
                          TOT          = WGTLLA(LA,IY)
                          AINNER(LA,0) = AINNER(LA,0) + TOT
                          AINNER(LA,1) = AINNER(LA,1) + TOT*TERM
                          AINNER(LA,2) = AINNER(LA,2) + TOT*TERM2
                          AINNER(LA,3) = AINNER(LA,3) + TOT*TERM3
                          AINNER(LA,4) = AINNER(LA,4) + TOT*TERM4
 4945                 CONTINUE
 4946             CONTINUE
                  DO 4947 LA=0,NLAS
                      W           = WGTGLE(LA,IX)
                      ACCUM(LA,0) = ACCUM(LA,0) + W*AINNER(LA,0)
                      ACCUM(LA,1) = ACCUM(LA,1) + W*AINNER(LA,1)
                      ACCUM(LA,2) = ACCUM(LA,2) + W*AINNER(LA,2)
                      ACCUM(LA,3) = ACCUM(LA,3) + W*AINNER(LA,3)
                      ACCUM(LA,4) = ACCUM(LA,4) + W*AINNER(LA,4)
 4947             CONTINUE
 4949         CONTINUE
              GOTO 5400
C
C                This is the common case, which should be reasonably
C                efficient for these longer loops.
C
 4991     CONTINUE
C
C            Outer integration loop (d x). Inner integrals (d y) are
C            accumulated in AINNER.
C
              DO 4999 IX=1,IGO
                  DO 4993 IP=0,IPMAX
                      DO 4992 LA=0,NLAS
                          AINNER(LA,IP) = ZERO
 4992                 CONTINUE
 4993             CONTINUE
                  DO 4996 IY=1,ILO
                      TERM = ARXROL(IX) - ARYROM(IY)
                      DO 4995 LA=0,NLAS
                          TOT          = WGTLLA(LA,IY)
                          AINNER(LA,0) = AINNER(LA,0) + TOT
                          DO 4994 IP=1,IPMAX
                              TOT           = TOT * TERM
                              AINNER(LA,IP) = AINNER(LA,IP) + TOT
 4994                     CONTINUE 
 4995                 CONTINUE
 4996             CONTINUE
C
C                Inner integral is done at this point. If we had not yet
C                reduced Laguerre integration order, it is exact, too.
C
                  DO 4998 IP=0,IPMAX
                      DO 4997 LA=0,NLAS
                          ACCUM(LA,IP) = ACCUM(LA,IP) 
     .                                 + WGTGLE(LA,IX)*AINNER(LA,IP)
 4997                 CONTINUE
 4998             CONTINUE
 4999         CONTINUE
              GOTO 5400
C
C            End of horror code
C
 5400         CONTINUE
C
C        Move integrals to their final resting place.
C
          DO 6500 IP=0,IPMAX
              DO 6490 LA=0,NLAS
                  AINT(LAMIN+LA,IP,IRO) = ACCUM(LA,IP)
 6490         CONTINUE
 6500     CONTINUE
 7000 CONTINUE
C
      RETURN
11010 FORMAT(' FIRST LEADING DIMENSION OF AINT (',I4,
     .       ') IS TOO SMALL IN PS3BI, NEED AT LEAST ',I4)
11020 FORMAT(' SECOND LEADING DIMENSION OF AINT (',I4,
     .       ') IS TOO SMALL IN PS3BI, NEED AT LEAST ',I4)
11030 FORMAT(' NUMBER OF LAMBDA VALUES (',I4,
     .       ') IS TOO BIG IN PS3BI, LIMIT IS ',I4)
11050 FORMAT(' IPMAX (',I4,') IS TOO LARGE IN PS3BI, LIMIT IS ',I4)
11060 FORMAT(' UNSUPPORTED NLAS,IPMAX (',I4,',',I4,
     .       ') ENCOUNTERED IN THE LINEARIZED LOOP IN PS3BI.' )
      END
C
      SUBROUTINE PS3BJ(LAMIN,LAMAX,IQ,ZETA,NPOINT,ROLESS,ROMORE,
     .                 AINT,LDA)
C
C   Compute batch of the outer integrals in BFH expansion.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      LAMIN  - Minimum lambda parameter value.
C      LAMAX  - Maximum lambda parameter value.
C      IQ     - q parameter value.
C      ZETA   - Exponential parameter.
C      NPOINT - Number of ro points.
C      ROLESS - Argument values, at least NPOINT entries.
C               Grouping same ro< values together will improve
C               performance.
C      ROMORE - Argument values, at least NPOINT entries.
C               Grouping same ro> values together will improve
C               performance.
C      AINT   - (Output) integrals, two-index array.
C               First index:  LAMIN to LAMAX, lambda parameter values.
C               Second index: 1 to NPOINT, ro argumnet values.
C      LDA1   - Leading dimension of AINT.
C
C   Accessed common blocks:
C
C      PSPRT  - Printing unit.
C
C   Local storage:
C
C      (MAXQ+1)*(MAXPT*MAXLAS+2) + 2*(MAXFAC+1) DOUBLE PRECISION words.
C      With MAXQ=7, MAXLAS=4, MAXPT=80, and MAXFAC=100, it is some 2800
C      DP words.
C
C   Module logic:
C
C           min(q,2*la+1) / q \   k!   / 2*la+1 \
C    I     =   Sum        |   | -----  |        |  I
C     q,la     k=0        \ k / zeta^k \   k    /   la,q-k
C
C   Precision:
C
C      No significant loss of precision should occur in PS3BJ.
C
C   Possible optimizations:
C
C      Don't bother. This is a minor time consumer compared to
C      PS3BI.
C
C   Bugs:
C
C      Number of lambda values which can be handled at the same time
C      is limited to MAXLAS. IQ value is limited to MAXQ. Number of
C      radial arguments is limited to MAXPT. LAMAX is limited to
C      (MAXFAC-1)/2 by the factorial and binomial tables.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXLAS=5)
      PARAMETER (MAXQ  =7)
      PARAMETER (MAXPT =80)
      PARAMETER (MAXFAC=100)
      PARAMETER (ZERO=0.D0)
      PARAMETER (ONE =1.D0)
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
C
      DIMENSION ROLESS(NPOINT), ROMORE(NPOINT)
      DIMENSION AINT(LAMIN:LAMIN+LDA-1,NPOINT)
C    Local temporaries
      DIMENSION AI(0:MAXLAS-1,0:MAXQ,MAXPT)
      DIMENSION FACT(0:MAXFAC), RFACT(0:MAXFAC)
      DIMENSION QSCALE(0:MAXQ), SCALE(0:MAXQ)
C    Range checking
      IF( LAMAX-LAMIN+1.GT.LDA ) THEN
          WRITE(NB6,11010) LDA, LAMAX-LAMIN+1
          STOP 'PS3BJ'
      ENDIF
      IF( LAMAX-LAMIN+1.GT.MAXLAS ) THEN
          WRITE(NB6,11020) LAMAX-LAMIN+1, MAXLAS
          STOP 'PS3BJ'
      ENDIF
      IF( IQ.GT.MAXQ ) THEN
          WRITE(NB6,11030) IQ, MAXQ
          STOP 'PS3BJ'
      ENDIF
      IF( NPOINT.GT.MAXPT ) THEN
          WRITE(NB6,11040) NPOINT, MAXPT
          STOP 'PS3BJ'
      ENDIF
      NFACT = MAX(IQ,2*LAMAX+1)
      IF( NFACT.GT.MAXFAC ) THEN
          WRITE(NB6,11050) NFACT, MAXFAC
          STOP 'PS3BJ'
      ENDIF
C    Evaluate inner integrals
      CALL PS3BI(LAMIN,LAMAX,IQ,ZETA,NPOINT,ROLESS,ROMORE,
     .           AI,MAXLAS,MAXQ+1)
C    Evaluate factorial values.
      FACT (0) = ONE
      RFACT(0) = ONE
      DO 1000 I=1,NFACT
          FACT (I) = FACT(I-1)*I
          RFACT(I) = ONE/FACT(I)
 1000 CONTINUE
C    Evaluate parts of the coefficients which is independent of lambda
      RZETA     = ONE/ZETA
      POWZ      = ONE
      QSCALE(0) = ONE
      DO 2000 K=1,IQ
          POWZ = POWZ*RZETA
          QSCALE(K) = FACT(IQ)*RFACT(IQ-K)*RFACT(K)*POWZ
 2000 CONTINUE
C    Run over possible lambda values
      DO 7000 LA=0,LAMAX-LAMIN
          ILA  = LAMIN+LA
          KTOP = MIN(IQ,2*ILA+1)
C        Evaluate parts of the coefficients which depend on lambda
          DO 3000 K=0,KTOP
              SCALE(K) = QSCALE(K)*FACT(2*ILA+1)*RFACT(2*ILA+1-K)
 3000     CONTINUE
C        Run over radial argument values, computing total integral in each case
          DO 5000 IPT=1,NPOINT
              SUM = ZERO
              DO 4000 K=0,KTOP
                  SUM = SUM + SCALE(K)*AI(LA,IQ-K,IPT)
 4000         CONTINUE
              AINT(ILA,IPT) = SUM
 5000     CONTINUE
 7000 CONTINUE
C
      RETURN
11010 FORMAT(' LEADING DIMENSION OF AINT (',I4,
     .       ') IS TOO SMALL IN PS3BJ. NEED AT LEAST ',I4)
11020 FORMAT(' TOO MUCH LABMDA VALUES IN PS3BJ (',I4,'). LIMIT IS ',I4)
11030 FORMAT(' IQ PARAMETER VALUE (',I4,
     .       ') TOO LARGE IN PS3BJ. LIMIT IS ',I4)
11040 FORMAT(' NUMBER OF RADIAL ARGUMENTS (',I4,
     .       ') TOO LARGE IN PS3BJ. LIMIT IS ',I4)
11050 FORMAT(' TOO MUCH FACTORIALS (',I4,
     .       ') ARE NEEDED IN PS3BJ. LIMIT IS ',I4)
      END
C
      SUBROUTINE PS3BA(N,L,ZETA,LS,A,NR,RS,AJNT,LDG,KHAVE,FUN,LDF)
C
C   Compute values of the Loewdin/Sharma expansion of a given order
C   at a set of points using BFH expression. This function is intended
C   to serve as a direct replacement for PS3SHQ.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      N,L    - Major and orbital quantum number of the orbital
C               being expanded. L should be smaller than N.
C      ZETA   - Orbital exponent.
C      LS     - Order of the expansion term
C      A      - Separation of the centers. A should be non-zero.
C      NR     - Number of radial points where function should be
C               computed
C      RS     - Array of radial arguments. Arguments should be
C               in an increasing order for the best performance.
C      AJNT   - Input/Output array of the outer auxiliary integrals.
C               First index: lambda parameter value, 0-based.
C               Second index: argument value, one-to one correspondence
C                   with arguments.
C      LDG    - Leading dimension of the AJNTS array.
C      KHAVE  - Input/Output: maximum order of the auxiliary integral in
C                   the AJNTS plus 1, should not exceed LDG.
C      FUN    - Two-dimensional array for computed expansion
C               function values.
C               First dimension:
C                  1 = value of the expansion function
C                  2 = bound for the error in the expansion function
C               Second dimension: absolute value of the magnetic 
C                  quantum number M, in the order 0 to L (i.e., at 
C                  least L + 1 entries). 
C               Third dimension: radial argument values (i.e.,
C                  at least NR entries)
C      LDF    - Second leading dimension of the FUN matrix.
C      
C   Accessed common blocks:
C
C      PSPRT  - Printing unit.
C
C   Local storage:
C
C      MAXNR*(3*MAXL+4) + (MAXL+1)*(MAXL+4) + 2*(MAXLS+MAXL) DOUBLE PRECISION
C      words. With MAXNR=80, MAXL=4, and MAXLS=40, it is some 1400 DP words.
C
C   Module logic:
C
C      This one is a fairly hairy bit:
C
C          N+M  L
C      (-1)    A  Sqrt[2*L+1] Sqrt[2*LS+1] Sqrt[(L+M)!] Sqrt[(L-M)!]       (A)
C
C                          
C           L                    1                          1      lap
C          Sum   --------------------------------------- -------  r        (B)
C        lap=|M|  (L-lap)! Sqrt[(lap+M)!] Sqrt[(lap-M)!]     lap
C                                                        (-A) 
C                               2*la+1  la
C         LS+lap   (2*la+1) zeta       A      la+1
C          Sum    -------------------------  r                             (C)
C      la=|LS-lap|            2  la
C                      2 (la!)  4
C
C        /  lap  la LS \  /  lap  la  LS \
C        |             |  |              |  J      (zeta,roless,romore)
C        \   M   0  -M /  \   0   0   0  /   N-L,la
C
C       Factors (A), (B), and (C), which are independent of radial parameter r,
C       are precomputed outside of the inner loop.
C
C   Precision:
C
C      EPS1 and EPS2 reflect the idea of the precision of the results
C      rather than control it. Precision is ultimately controlled by
C      the selection of the integration orders in PS3BI above.
C
C      Note that if the precision is estimated by comparison with PS3SHQ,
C      expansion Sharma's Bnu coefficients should be computed in full
C      quadruple precision.
C
C   Possible optimizations:
C
C      It might be a good idea to play with INCLA a little bit.
C      Running with larger INCLA will result in larger bodies of
C      inner loops in PS3BI. However, unrolled loops are supplied
C      for INCLA up to 3, so larger increments will degrade performance.
C      Additionally, non-zero INCLA value will result in unnecessary
C      evaluation of some integrals, also decreasing performance.
C
C      Going from INCLA=0 to INCLA=3 decreases execution times by a
C      factor of 2 on R8k, and by 1.5 on R4k. R10k and Power1 are
C      less sensitive (about 20% improvement). YMMV.
C
C   Bugs:
C
C      L is limited by MAXL. LS is limited by LS. NR is limited by MAXNR.
C      There is no reason why any of these cannot be increased, if desired.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXL  =4)
      PARAMETER (MAXLS =40)
      PARAMETER (MAXNR =80)
      PARAMETER (INCLA =3)
      PARAMETER (EPS1  =1.D-7)
      PARAMETER (EPS2  =1.D-12)
      PARAMETER (ZERO  =0.D0)
      PARAMETER (ONE   =1.D0)
      PARAMETER (HALF  =0.5D0)
      PARAMETER (FOURTH=0.25D0)
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
C
      EXTERNAL PS3J
      DIMENSION RS(NR), FUN(2,0:LDF-1,NR), AJNT(0:LDG-1,NR)
C    Local temporaries.
      DIMENSION FACT(0:MAXLS+MAXL), RFACT(0:MAXLS+MAXL)
      DIMENSION ROLESS(MAXNR), ROMORE(MAXNR)
      DIMENSION SCALEA(0:MAXL), SCALEB(0:MAXL,0:MAXL)
      DIMENSION SCALEC(-MAXL:MAXL)
      DIMENSION RPOWL(0:MAXL,MAXNR), RPOWH(-MAXL+1:MAXL+1,MAXNR)
C    Range checking.
C
      IF( N.LE.0 .OR. L.LT.0 .OR. L.GE.N ) THEN
          WRITE(NB6,11010) N, L
          STOP 'PS3BA'
      ENDIF
      IF( L.GT.MAXL ) THEN
          WRITE(NB6,11020) L, MAXL
          STOP 'PS3BA'
      ENDIF
      IF( ZETA.LE.ZERO ) THEN
          WRITE(NB6,11030) ZETA
          STOP 'PS3BA'
      ENDIF
      IF( LS.LT.0 .OR. LS.GT.MAXLS ) THEN
          WRITE(NB6,11035) LS, MAXLS
          STOP 'PS3BA'
      ENDIF
      IF( A.EQ.ZERO ) THEN
          WRITE(NB6,11040)
          STOP 'PS3BA'
      ENDIF
      IF( NR.LE.0 .OR. NR.GT.MAXNR ) THEN
          WRITE(NB6,11050) NR, MAXNR
          STOP 'PS3BA'
      ENDIF
      IF( LDF.LE.L ) THEN
          WRITE(NB6,11060) LDF, L+1
          STOP 'PS3BA'
      ENDIF
C    Evaluate factorials
      FACT (0) = ONE
      RFACT(0) = ONE
      DO 500 I=1,L+MAX(L,LS)
          FACT (I) = FACT(I-1)*I
          RFACT(I) = ONE/FACT(I)
  500 CONTINUE
C    Evaluate necesary J integrals
      MAXLA = LS + L
      IF( KHAVE-1.LT.MAXLA ) THEN
          IF( MAXLA+1.GT.LDG ) THEN
              WRITE(NB6,11070) MAXLA+1, LDG
              STOP 'PS3BA'
          ENDIF
C        First, we'd need roless and romore, rather than A's and R's
          ABSA = ABS(A)
          DO 600 I=1,NR
              ROLESS(I) = MIN(ABSA,RS(I))
              ROMORE(I) = MAX(ABSA,RS(I))
  600     CONTINUE
          MAXLAT = MAX(KHAVE+INCLA,MAXLA)
          MAXLAT = MIN(MAXLAT,LDG-1)
C        Compute integrals and update cache's upper watermark
          CALL PS3BJ(KHAVE,MAXLAT,N-L,ZETA,NR,ROLESS,ROMORE,
     .               AJNT(KHAVE,1),LDG)
          KHAVE = MAXLAT + 1
      ENDIF
C    Compute common pre-factors (A)
      SCALE = (-1)**N * A**L * SQRT(DBLE((2*L+1)*(2*LS+1)))
      DO 1000 M=0,L
          SCALEA(M) = SCALE * (-1)**M * SQRT(FACT(L+M)*FACT(L-M))
 1000 CONTINUE
C    Compute lambda-prim dependent pre-factors (B)
      DO 1500 LAP=0,L
          SCALE = RFACT(L-LAP) / (-A)**LAP
          DO 1400 M=0,MIN(L,LAP)
              SCALEB(M,LAP) = SCALE*SQRT(RFACT(LAP+M)*RFACT(LAP-M))
 1400     CONTINUE
 1500 CONTINUE
C    Compute lambda dependent per-factors (C)
      LAMIN = MAX(0,LS-L)
      DO 2000 LA=LAMIN,LS+L
          SCALEC(LA-LS) = HALF * (2*LA+1) * ZETA**(2*LA+1) 
     .                         * (FOURTH*A)**LA * RFACT(LA)**2 
 2000 CONTINUE
C    Compute low powers of radial parameters
      DO 2500 IPT=1,NR
          RPOWL(0,IPT) = ONE
          DO 2400 LAP=1,L
              RPOWL(LAP,IPT) = RPOWL(LAP-1,IPT)*RS(IPT)
 2400     CONTINUE
 2500 CONTINUE
C    Compute high powers of radial parameters
      DO 3000 IPT=1,NR
          RPOWH(LAMIN+1-LS,IPT) = RS(IPT)**(LAMIN+1)
          DO 2900 LA=LAMIN+1,LS+L
              RPOWH(LA+1-LS,IPT) = RPOWH(LA-LS,IPT)*RS(IPT)
 2900     CONTINUE
 3000 CONTINUE
C    First element of each function value pair will accumulate
C    actual function value, while the second element will get 
C    sum of the absolute values of terms.
      DO 4000 IPT=1,NR
          DO 3900 M=0,L
              FUN(1,M,IPT) = ZERO
              FUN(2,M,IPT) = ZERO
 3900     CONTINUE
 4000 CONTINUE
C    Run over all possible values of lambda and lambda prim, and 
C    acumulate contributions to the approapriate function values.
C    Loop over LA goes by 2 because C3JA is zero otherwise.
C    There is a specialized version of PS3J (PSVXFI) geared for
C    3j coefficients with zero Ms, but it was never tested for 
C    large Ls.
      DO 8000 LAP=0,L
          DO 7000 LA=ABS(LS-LAP),LS+LAP,2
              C3JA = PS3J(LAP,LA,LS,0,0,0)
              DO 6000 M=0,MIN(LAP,L,LS)
                  IF( M.EQ.0 ) THEN
                      C3JB = C3JA
                  ELSE
                      C3JB = PS3J(LAP,LA,LS,M,0,-M)
                  ENDIF
                  SCALE = SCALEA(M)*SCALEB(M,LAP)*SCALEC(LA-LS)
     .                             *C3JA*C3JB
                  DO 5000 IPT=1,NR
                      TERM = SCALE*RPOWL(LAP,IPT)
     .                            *RPOWH(LA+1-LS,IPT)*AJNT(LA,IPT)
                      FUN(1,M,IPT) = FUN(1,M,IPT) + TERM
                      FUN(2,M,IPT) = FUN(2,M,IPT) + ABS(TERM)
 5000             CONTINUE
 6000         CONTINUE
 7000     CONTINUE
 8000 CONTINUE
C    Convert absolute sums into error estimations, and away we go
      DO 9000 IPT=1,NR
          DO 8900 M=0,L
              FUN(2,M,IPT) = MAX(EPS1*FUN(2,M,IPT),EPS2)
 8900     CONTINUE
 9000 CONTINUE
      RETURN
11010 FORMAT(' QUANTUM NUMBERS N (',I4,') AND L (',I4,
     .       ') ARE NON-PHYSICAL IN PS3BA.')
11020 FORMAT(' QUANTUM NUMBER L (',I4,
     .       ') IS TOO LARGE IN PS3BA. LIMIT IS ',I4)
11030 FORMAT(' ORBITAL EXPONENT ZETA (',G15.7,
     .       ') IS NON-POSITIVE IN PS3BA.')
11035 FORMAT(' EXPANSION ORDER LS (',I4,
     .       ') IS TOO LARGE IN PS3BA. LIMIT IS ',I4)
11040 FORMAT(' CENTERS SEPARATION A IS ZERO IN PS3BA.')
11050 FORMAT(' NUMBER OF RADIAL POINTS NR (',I4,
     .       ') IS TOO LARGE IN PS3BA. LIMIT IS ',I4)
11060 FORMAT(' LEADING DIMENSION OF OUTPUT ARRAY LDF (',I4,
     .       ') IS TOO SMALL. AT LEAST ',I4,' IS NEEDED.')
11070 FORMAT(' AJNT IS TOO SMALL TO HOLD ALL LAMBDAS. NEED ',I4,
     .       ', BUT ONLY ',I4,' IS AVAILABLE.' )
      END
