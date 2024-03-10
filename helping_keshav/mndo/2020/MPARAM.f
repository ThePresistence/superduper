      SUBROUTINE MPARAM
C     *
C     MINDO/3 PARAMETERS.
C     *
      USE LIMIT, ONLY: LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./DIPOL1/ HYF(LMZ),HYFPD(LMZ)
     ./MNDOBL/ USSM(LMZ),UPPM(LMZ),ZSM(LMZ),ZPM(LMZ),
     .         BETASM(LMZ),BETAPM(LMZ),ALPM(LMZ),EISOLM(LMZ),
     .         DDM(LMZ),QQM(LMZ),AMM(LMZ),ADM(LMZ),AQM(LMZ),
     .         GSSM(LMZ),GPPM(LMZ),GSPM(LMZ),GP2M(LMZ),HSPM(LMZ)
     ./MNDOEL/ IMM(LMZ)
     ./PARAVL/ IMPAR(LMZ)
     ./PARDER/ CORE(LMZ),EHEAT(LMZ),EISOL(LMZ)
     ./PARMIN/ F0(18),VS(18),VP(18),BETA(55),ALP(55)
     ./PAROPT/ USS(LMZ),UPP(LMZ),ZS(LMZ),ZP(LMZ),BETAS(LMZ*3)
     ./REP   / GSS(LMZ),GPP(LMZ),GSP(LMZ),GP2(LMZ),HSP(LMZ),HPP(LMZ)
C     OPTIMIZED PARAMETERS.
      IMPAR(1)  = 1
      IMPAR(5)  = 1
      IMPAR(6)  = 1
      IMPAR(7)  = 1
      IMPAR(8)  = 1
      IMPAR(9)  = 1
      IMPAR(14) = 1
      IMPAR(15) = 1
      IMPAR(16) = 1
      IMPAR(17) = 1
      USS(1)    = -12.505D0
      USS(5)    = -33.61D0
      USS(6)    = -51.79D0
      USS(7)    = -66.06D0
      USS(8)    = -91.73D0
      USS(9)    =-129.86D0
      USS(14)   = -39.82D0
      USS(15)   = -56.23D0
      USS(16)   = -73.39D0
      USS(17)   = -98.99D0
      UPP(1)    = ZERO
      UPP(5)    = -25.11D0
      UPP(6)    = -39.18D0
      UPP(7)    = -56.40D0
      UPP(8)    = -78.80D0
      UPP(9)    =-105.93D0
      UPP(14)   = -29.15D0
      UPP(15)   = -42.31D0
      UPP(16)   = -57.25D0
      UPP(17)   = -76.43D0
      ZS(1)     = 1.3D0
      ZS(5)     = 1.211156D0
      ZS(6)     = 1.739391D0
      ZS(7)     = 2.704546D0
      ZS(8)     = 3.640575D0
      ZS(9)     = 3.111270D0
      ZS(14)    = 1.629173D0
      ZS(15)    = 1.926108D0
      ZS(16)    = 1.719480D0
      ZS(17)    = 3.430887D0
      ZP(1)     = ZERO
      ZP(5)     = 0.972826D0
      ZP(6)     = 1.709645D0
      ZP(7)     = 1.870839D0
      ZP(8)     = 2.168448D0
      ZP(9)     = 1.41986D0
      ZP(14)    = 1.381721D0
      ZP(15)    = 1.590665D0
      ZP(16)    = 1.403205D0
      ZP(17)    = 1.627017D0
      BETA(1)   = 0.244770D0
      BETA(2)   = 0.185347D0
      BETA(3)   = 0.151324D0
      BETA(4)   = 0.315011D0
      BETA(5)   = 0.250031D0
      BETA(6)   = 0.419907D0
      BETA(7)   = 0.360776D0
      BETA(8)   = 0.310959D0
      BETA(9)   = 0.410886D0
      BETA(10)  = 0.377342D0
      BETA(11)  = 0.417759D0
      BETA(12)  = 0.349745D0
      BETA(13)  = 0.464514D0
      BETA(14)  = 0.458110D0
      BETA(15)  = 0.659407D0
      BETA(16)  = 0.195242D0
      BETA(17)  = 0.219591D0
      BETA(18)  = 0.247494D0
      BETA(19)  = 0.205347D0
      BETA(20)  = 0.334044D0
      BETA(21)  = 0.197464D0
      BETA(22)  = 0.289647D0
      BETA(24)  = 0.411377D0
      BETA(26)  = 0.530779D0
      BETA(28)  = 0.291703D0
      BETA(29)  = 0.320118D0
      BETA(31)  = 0.457816D0
      BETA(36)  = 0.311790D0
      BETA(37)  = 0.220654D0
      BETA(39)  = 0.284620D0
      BETA(40)  = 0.313170D0
      BETA(41)  = 0.422890D0
      BETA(45)  = 0.202489D0
      BETA(46)  = 0.231653D0
      BETA(48)  = 0.315480D0
      BETA(49)  = 0.462366D0
      BETA(50)  = 0.474860D0
      BETA(53)  = 0.277322D0
      BETA(54)  = 0.221764D0
      BETA(55)  = 0.258969D0
      ALP(1)    = 1.489450D0
      ALP(2)    = 2.090352D0
      ALP(3)    = 2.280544D0
      ALP(4)    = 1.475836D0
      ALP(5)    = 2.138291D0
      ALP(6)    = 1.371208D0
      ALP(7)    = 0.589380D0
      ALP(8)    = 1.909763D0
      ALP(9)    = 1.635259D0
      ALP(10)   = 2.029618D0
      ALP(11)   = 0.478901D0
      ALP(12)   = 2.484827D0
      ALP(13)   = 1.820975D0
      ALP(14)   = 1.873859D0
      ALP(15)   = 1.537190D0
      ALP(16)   = 3.771362D0
      ALP(17)   = 2.862183D0
      ALP(18)   = 2.725913D0
      ALP(19)   = 2.861667D0
      ALP(20)   = 2.266949D0
      ALP(21)   = 3.864997D0
      ALP(22)   = 0.940789D0
      ALP(24)   = 1.101382D0
      ALP(26)   = 1.361185D0
      ALP(28)   = 0.918432D0
      ALP(29)   = 0.923170D0
      ALP(31)   = 1.029693D0
      ALP(36)   = 1.186652D0
      ALP(37)   = 1.700698D0
      ALP(39)   = 1.761370D0
      ALP(40)   = 1.878176D0
      ALP(41)   = 2.077240D0
      ALP(45)   = 1.751617D0
      ALP(46)   = 2.089404D0
      ALP(48)   = 1.676222D0
      ALP(49)   = 1.155635D0
      ALP(50)   = 1.534841D0
      ALP(53)   = 1.543720D0
      ALP(54)   = 1.950318D0
      ALP(55)   = 1.792125D0
C     DERIVED PARAMETERS.
      EHEAT(1)  =   52.102D0
      EHEAT(5)  =  135.7D0
      EHEAT(6)  =  170.89D0
      EHEAT(7)  =  113.0D0
      EHEAT(8)  =   59.559D0
      EHEAT(9)  =   18.86D0
      EHEAT(14) =  106.0D0
      EHEAT(15) =   79.8D0
      EHEAT(16) =   65.65D0
      EHEAT(17) =   28.95D0
      EISOL(1)  =  -12.505D0
      EISOL(5)  =  -61.70D0
      EISOL(6)  = -119.47D0
      EISOL(7)  = -187.51D0
      EISOL(8)  = -307.07D0
      EISOL(9)  = -475.00D0
      EISOL(14) =  -90.98D0
      EISOL(15) = -150.81D0
      EISOL(16) = -229.15D0
      EISOL(17) = -345.93D0
      F0(1)     = 12.848D0
      F0(5)     =  8.958D0
      F0(6)     = 10.833D0
      F0(7)     = 12.377D0
      F0(8)     = 13.985D0
      F0(9)     = 16.250D0
      F0(14)    =  7.57D0
      F0(15)    =  9.00D0
      F0(16)    = 10.20D0
      F0(17)    = 11.73D0
      VS(1)     =-13.605D0
      VS(5)     =-15.16D0
      VS(6)     =-21.34D0
      VS(7)     =-27.51D0
      VS(8)     =-35.30D0
      VS(9)     =-43.70D0
      VS(14)    =-17.82D0
      VS(15)    =-21.10D0
      VS(16)    =-23.84D0
      VS(17)    =-25.26D0
      VP(1)     =  0.0D0
      VP(5)     = -8.52D0
      VP(6)     =-11.54D0
      VP(7)     =-14.34D0
      VP(8)     =-17.91D0
      VP(9)     =-20.89D0
      VP(14)    = -8.51D0
      VP(15)    =-10.29D0
      VP(16)    =-12.41D0
      VP(17)    =-15.09D0
      HYF(5)    = 6.520587D0
      HYF(6)    = 4.253676D0
      HYF(7)    = 2.947501D0
      HYF(8)    = 2.139793D0
      HYF(9)    = 2.2210719D0
      HYF(14)   = 6.663059D0
      HYF(15)   = 5.657623D0
      HYF(16)   = 6.345552D0
      HYF(17)   = 2.522964D0
      BETAS(1)  = BETA(1)*VS(1)*TWO
      CORE(1)   = ONE
      CORE(2)   = TWO
      DO 10 I=3,10
      CORE(I) = I-2
   10 CONTINUE
      DO 20 I=11,18
      CORE(I) = I-10
   20 CONTINUE
      DO 30 I=1,17
      IF(IMM(I).EQ.1) THEN
         GSS(I)   = GSSM(I)
         GSP(I)   = GSPM(I)
         GPP(I)   = GPPM(I)
         GP2(I)   = GP2M(I)
         HSP(I)   = HSPM(I)
         HPP(I)   = PT5*(GPPM(I)-GP2M(I))
         HYFPD(I) = ZERO
      ENDIF
   30 CONTINUE
      RETURN
      END
