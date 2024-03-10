      FUNCTION RSC (K,NA,EA,NB,EB,NC,EC,ND,ED)
C     *
C     CALCULATE THE RADIAL PART OF ONE-CENTER TWO-ELECTRON INTEGRALS
C     (SLATER-CONDON PARAMETER).
c     K     - TYPE OF INTEGRAL, CAN BE EQUAL TO 0,1,2,3,4 IN SPD-BASIS
C     NA,NB - PRINCIPLE QUANTUM NUMBER OF AO, ELECTRON 1
C     EA,EB - EXPONENTS OF AO, ELECTRON 1
C     NC,ND - PRINCIPLE QUANTUM NUMBER OF AO, ELECTRON 2
C     EC,ED - EXPONENTS OF AO, ELECTRON 2
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./PSC   / F(30),B(30,30)
      AEA    = LOG(EA)
      AEB    = LOG(EB)
      AEC    = LOG(EC)
      AED    = LOG(ED)
      NAB    = NA+NB
      NCD    = NC+ND
      ECD    = EC+ED
      EAB    = EA+EB
      E      = ECD+EAB
      N      = NAB+NCD
      AE     = LOG(E)
      A2     = LOG(TWO)
      ACD    = LOG(ECD)
      AAB    = LOG(EAB)
      C      = EXP(F(N)+NA*AEA+NB*AEB+NC*AEC+ND*AED
     1             +PT5*(AEA+AEB+AEC+AED)+A2*(N+2)
     2             -PT5*(F(2*NA+1)+F(2*NB+1)
     3             +F(2*NC+1)+F(2*ND+1))-AE*N)
      C      = C*EV
      S0     = ONE/E
      S1     = ZERO
      S2     = ZERO
      M      = NCD-K
      DO 10 I=1,M
      S0     = S0*E/ECD
      S1     = S1+S0*(B(NCD-K,I)-B(NCD+K+1,I))/B(N,I)
   10 CONTINUE
      M1     = M+1
      M2     = NCD+K+1
      DO 20 I=M1,M2
      S0     = S0*E/ECD
      S2     = S2+S0*B(M2,I)/B(N,I)
   20 CONTINUE
      S3     = EXP(AE*N-ACD*M2-AAB*(NAB-K))/B(N,M2)
      RSC    = C*(S1-S2+S3)
      RETURN
      END
