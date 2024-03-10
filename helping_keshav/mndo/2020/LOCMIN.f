      SUBROUTINE LOCMIN (X,G,P,XX,GG,N,T,F,MAXLIN,IPRINT,
     1                   ARRAY,LM5,ICALL,SCFCAL)
C     *
C     QUADRATIC LINE SEARCH.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL SCFCAL
      COMMON
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./NBFILE/ NBF(20)
     ./OPTTOL/ TOLEND(4),TOLX,TOLF,TREL,TABS,XMAXST
     ./OVERLY/ IOV,JOV,KOV,LOV
      DIMENSION X(N),G(N),P(N),XX(N),GG(N)
      DIMENSION ARRAY(LM5)
      DIMENSION PHI(3),VT(3)
C *** FILE NUMBERS.
      NB6    = NBF(6)
C *** INITIALIZATION.
      PMAX   = ZERO
      DO 10 I=1,N
      PMAX   = MAX(PMAX,ABS(P(I)))
   10 CONTINUE
      TMAX   = XMAXST/PMAX
      FIN    = F
      PHI(1) = F
      VT(1)  = ZERO
      CALL EXCHNG(F,FF,VT(1),TT,X,XX,N)
C *** SECOND POINT IN LINE SEARCH.
      VT(2)  = T*PT25
      IF(VT(2).GT.TMAX) VT(2)=TMAX
      T      = VT(2)
      DO 20 I=1,N
      X(I)   = XX(I)+T*P(I)
   20 CONTINUE
      LOV    = 1
      CALL COMPFG(N,X,PHI(2),FF,GG,ARRAY,LM5,ICALL,SCFCAL)
      IF(ICALL.EQ.-1) RETURN
      IF(PHI(1).LE.PHI(2)) THEN
         VT(3)  = -VT(2)
         LEFT   = 3
         LCENTR = 1
         LRIGHT = 2
      ELSE
         VT(3)  = TWO*VT(2)
         LEFT   = 1
         LCENTR = 2
         LRIGHT = 3
         CALL EXCHNG(PHI(2),FF,T,TT,GG,G,N)
      ENDIF
C *** THIRD POINT IN LINE SEARCH.
      TLAST  = VT(3)
      T      = VT(3)
      DO 40 I=1,N
      X(I)   = XX(I)+T*P(I)
   40 CONTINUE
      LOV    = 2
      CALL COMPFG(N,X,F,FF,GG,ARRAY,LM5,ICALL,SCFCAL)
      IF(ICALL.EQ.-1) RETURN
      IF(F.LT.FF) CALL EXCHNG(F,FF,T,TT,GG,G,N)
      PHI(3) = F
      IF(IPRINT.GT.0) WRITE(NB6,500) VT(1),PHI(1),VT(2),PHI(2),VT(3),
     1                             PHI(3)
C *** LOOP OVER NEW POINTS IN LINE SEARCH.
      FLAST  = FF
      ICTR   = 2
   60 ALPHA  = VT(2)-VT(3)
      BETA   = VT(3)-VT(1)
      GAMMA  = VT(1)-VT(2)
      ABG    = ALPHA*BETA*GAMMA
      IF(ABS(ABG).GT.1.0D-20) THEN
         ALPHA = -(PHI(1)*ALPHA+PHI(2)*BETA+PHI(3)*GAMMA)/ABG
      ELSE
         GO TO 200
      ENDIF
      BETA   = ((PHI(1)-PHI(2))/GAMMA)-ALPHA*(VT(1)+VT(2))
      IF(ALPHA.LE.ZERO) THEN
         IF(PHI(LRIGHT).LE.PHI(LEFT)) THEN
            T = THREE*VT(LRIGHT)-TWO*VT(LCENTR)
         ELSE
            T = THREE*VT(LEFT)-TWO*VT(LCENTR)
         ENDIF
         S   = T-TLAST
         IF(ABS(S).GT.TMAX) S=SIGN(TMAX,S)
         T   = S+TLAST
      ELSE
         T   = -BETA/(TWO*ALPHA)
         S   = T-TLAST
         XXM = TWO*TMAX
         IF(ABS(S).GT.XXM) S=SIGN(XXM,S)
         T   = S+TLAST
      ENDIF
C     CHECK WHETHER NEW POINT SHOULD STILL BE CALCULATED.
      IF(FF.LT.FIN) THEN
         IF(ABS(S).LT.(TREL*ABS(T+TLAST)+TABS)) GO TO 200
         IF(ABS(S*PMAX).LT.TOLX .AND. ICTR.GT.2) GO TO 200
      ENDIF
C     CALCULATION FOR NEW POINT.
      DO 70 I=1,N
      X(I)   = XX(I)+T*P(I)
   70 CONTINUE
      LOV    = 3
      CALL COMPFG(N,X,F,FF,GG,ARRAY,LM5,ICALL,SCFCAL)
      IF(ICALL.EQ.-1) RETURN
      IF(F.LT.FF) CALL EXCHNG(F,FF,T,TT,GG,G,N)
      ICTR   = ICTR+1
      IF(IPRINT.GT.0) WRITE(NB6,501) VT(LEFT),PHI(LEFT),VT(LCENTR),
     1  PHI(LCENTR),VT(LRIGHT),PHI(LRIGHT),T,F
C     TEST FOR CONVERGENCE.
      IF(FF.LT.FIN) THEN
         IF(ABS(T-TLAST).LT.(ABS(T+TLAST)*TREL+TABS)) GO TO 200
         IF(ABS(F-FLAST).LT.TOLF) GO TO 200
      ENDIF
      FLAST  = F
      TLAST  = T
C     SORT CURRENT POINTS.
      IF((T.GT.VT(LRIGHT)) .OR. (T.GT.VT(LCENTR) .AND. F.LT.
     1  PHI(LCENTR)) .OR. (T.GT.VT(LEFT) .AND. T.LT.VT(LCENTR) .AND.
     2  F.GT.PHI(LCENTR))) GO TO 90
      VT(LRIGHT) = T
      PHI(LRIGHT) = F
      GO TO 100
   90 VT(LEFT) = T
      PHI(LEFT) = F
  100 IF(VT(LCENTR).LT.VT(LRIGHT))  GO TO 110
      I      = LCENTR
      LCENTR = LRIGHT
      LRIGHT = I
  110 IF(VT(LEFT).LT.VT(LCENTR))  GO TO 120
      I      = LEFT
      LEFT   = LCENTR
      LCENTR = I
  120 IF(VT(LCENTR).LT.VT(LRIGHT))  GO TO 130
      I      = LCENTR
      LCENTR = LRIGHT
      LRIGHT = I
  130 IF(ICTR.LT.MAXLIN) GO TO 60
C *** EXIT AFTER FINISHING THE LINE SEARCH.
  200 IF(T.NE.TT) THEN
         DO 210 I=1,N
         X(I) = XX(I)+TT*P(I)
  210    CONTINUE
      ENDIF
      F      = FF
      T      = TT
      IF(F.LT.FIN .AND. T.LT.ZERO) THEN
         T   = -T
         DO 220 I=1,N
         P(I) = -P(I)
  220    CONTINUE
      ENDIF
      RETURN
  500 FORMAT(/1X,'////LOCMIN',11X,'ALPHA',13X,'F(X)'// 5X,'LEFT   ...',
     1           2E17.8/5X,'CENTER ...',2E17.8/5X,'RIGHT  ...',2E17.8/)
  501 FORMAT(5X,'LEFT   ...',2E17.8/5X,'CENTER ...',2E17.8/5X,
     1          'RIGHT  ...',2E17.8/5X,'NEW    ...',2E17.8/)
      END
