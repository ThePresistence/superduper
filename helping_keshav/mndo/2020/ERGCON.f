      SUBROUTINE ERGCON(EA,EB,EIG,IMOCI,LM3,II,IJ,IK,IL,IU,IV,IW,IX,IOP,
     1                  G,LM10)
C     *
C     ENERGY OF AN EXCITED CONFIGURATION DEFINED BY EXCITATION INDICES
C     II->IU, IJ->IV, IK->IW, IL->IX.
C     *
C     NOTE. THE ARRAY G(I,J) CONTAINS THE COULOMB (I.GE.J) AND
C     EXCHANGE (I.LT.J) INTEGRALS. FOR COMPATIBILITY WITH THIS
C     STORAGE, THE EXCITATION INDICES MUST FULFIL.
C     II.GE.IJ.GE.IK.GE.IL, AND IU.GE.IV.GE.IW.GE.IX.
C     ALL INDICES FOR OCCUPIED MOS MUST BE SMALLER THAN THOSE
C     FOR UNOCCUPIED MOS, I.E. II.LT.IU ETC.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./NBFILE/ NBF(20)
      DIMENSION EIG(LM3),IMOCI(LM3)
      DIMENSION G(LM10,LM10)
      DATA ONE5/1.5D0/
      IF(II.LT.IJ .OR. IU.LT.IV) GO TO 500
C     INITIALIZATION.
      III    = IMOCI(II)
      IIJ    = IMOCI(IJ)
      IIU    = IMOCI(IU)
      IIV    = IMOCI(IV)
      IF(IL.GT.0 .AND. IX.GT.0) GO TO 300
      IF(IK.GT.0 .AND. IW.GT.0) GO TO 200
C *** DOUBLE EXCITATIONS.
      EA     = EIG(IIU)+EIG(IIV)-EIG(III)-EIG(IIJ)
      EB     = EA
      IF(IOP.EQ.1 .OR. IOP.EQ.3) RETURN
      EA     = EA+G(II,IJ)+G(IU,IV)-G(IU,II)-G(IV,II)-G(IU,IJ)-G(IV,IJ)
  100 EE     = EA
      IF(II.EQ.IJ) THEN
         IF(IU.EQ.IV) THEN
C           CASE II,II-IU,IU (TYPE 2).
            EA  = EA+TWO*G(II,IU)
         ELSE
C           CASE II,II-IU,IV (TYPE 3).
            EA  = EA+G(II,IU)+G(II,IV)+G(IV,IU)
         ENDIF
      ELSE
         IF(IU.EQ.IV) THEN
C           CASE II,IJ-IU,IU (TYPE 4).
            EA  = EA+G(II,IU)+G(IJ,IU)+G(IJ,II)
         ELSE
C           CASE II,IJ-IU,IV (TYPE 5,6).
            GS1 = G(II,IU)+G(II,IV)+G(IJ,IU)+G(IJ,IV)
            GS2 = G(IJ,II)+G(IV,IU)
            EA  = EE+PT5*GS1+GS2
            EB  = EE+ONE5*GS1-GS2
         ENDIF
      ENDIF
      RETURN
C *** TRIPLE EXCITATIONS.
  200 IF(IJ.LT.IK .OR. IV.LT.IW) GO TO 500
      IIK    = IMOCI(IK)
      IIW    = IMOCI(IW)
      EA     = EIG(IIU)+EIG(IIV)+EIG(IIW)-EIG(III)-EIG(IIJ)-EIG(IIK)
      EB     = EA
      IF(IOP.EQ.1 .OR. IOP.EQ.3) RETURN
      EA     = EA+G(II,IJ)+G(II,IK)+G(IJ,IK)+G(IU,IV)+G(IU,IW)+G(IV,IW)
      EA     = EA-G(IU,II)-G(IU,IJ)-G(IU,IK)-G(IV,II)-G(IV,IJ)-G(IV,IK)
      EA     = EA-G(IW,II)-G(IW,IJ)-G(IW,IK)
      EE     = EA
      IF(IJ.EQ.IK) THEN
         IF(IV.EQ.IW) THEN
C           CASE II,IJ,IJ-IU,IV,IV (TYPE 8C).
            EA  = EA-G(IJ,II)-G(IV,IU)+TWO*G(II,IU)+G(II,IV)+G(IJ,IU)
     1            +TWO*G(IJ,IV)
         ELSE
C           CASE II,IJ,IJ-IU,IV,IW (TYPE 9,10).
            IF(IU.NE.IV) THEN
               GS1 = -G(IJ,II)+PT5*(G(II,IU)+G(II,IV))+G(IJ,IU)
     1               +G(IJ,IV)+G(IJ,IW)
               GS2 = -PT5*(G(IW,IV)+G(IW,IU))+G(IV,IU)
               EA  = EE+GS1+GS2+TWO*G(II,IW)
               EB  = EE+GS1-GS2+G(II,IV)+G(II,IU)
            ELSE
C              CASE II,IJ,IJ-IU,IU,IW (TYPE 8B).
               EA  = EA-G(IJ,II)-G(IW,IU)+TWO*G(II,IW)+G(II,IU)
     1               +G(IJ,IW)+TWO*G(IJ,IU)
            ENDIF
         ENDIF
      ELSE
         IF(IV.NE.IW) GO TO 500
         IF(II.NE.IJ) THEN
C           CASE II,IJ,IK-IU,IV,IV (TYPE 12,13).
            GS1 = -G(IV,IU)+PT5*(G(IK,IU)+G(IJ,IU))+G(II,IV)
     1            +G(IJ,IV)+G(IK,IV)
            GS2 = -PT5*(G(IK,II)+G(IJ,II))+G(IK,IJ)
            EA  = EE+GS1+GS2+TWO*G(II,IU)
            EB  = EE+GS1-GS2+G(IK,IU)+G(IJ,IU)
         ELSE
C           CASE II,II,IK-IU,IV,IV (TYPE 8B).
            EA  = EA-G(IK,II)-G(IV,IU)+TWO*G(IK,IU)
     1            +G(IK,IV)+G(II,IU)+TWO*G(II,IV)
         ENDIF
      ENDIF
      RETURN
C *** QUADRUPLE EXCITATIONS.
  300 IF(IK.NE.IL .OR. IW.NE.IX) GO TO 500
      IIK    = IMOCI(IK)
      IIL    = IMOCI(IL)
      IIW    = IMOCI(IW)
      IIX    = IMOCI(IX)
      EA     = EIG(IIU)+EIG(IIV)+EIG(IIW)+EIG(IIX)-EIG(III)-EIG(IIJ)
     1        -EIG(IIK)-EIG(IIL)
      EB     = EA
      IF(IOP.EQ.1 .OR. IOP.EQ.3) RETURN
      EA     = EA+G(II,IJ)+G(II,IK)+G(II,IL)+G(IJ,IK)+G(IJ,IL)+G(IK,IL)
      EA     = EA+G(IU,IV)+G(IU,IW)+G(IU,IX)+G(IV,IW)+G(IV,IX)+G(IW,IX)
      EA     = EA-G(IU,II)-G(IU,IJ)-G(IU,IK)-G(IU,IL)
      EA     = EA-G(IV,II)-G(IV,IJ)-G(IV,IK)-G(IV,IL)
      EA     = EA-G(IW,II)-G(IW,IJ)-G(IW,IK)-G(IW,IL)
      EA     = EA-G(IX,II)-G(IX,IJ)-G(IX,IK)-G(IX,IL)
      EA     = EA-G(IK,II)-G(IK,IJ)-G(IW,IU)-G(IW,IV)
     1           +G(II,IW)+G(IJ,IW)
     2           +G(IK,IU)+G(IK,IV)+TWO*G(IK,IW)
      GO TO 100
C *** ERROR EXIT.
  500 CONTINUE
      NB6 = NBF(6)
      WRITE(NB6,510) II,IJ,IK,IL,IU,IV,IW,IX
      STOP 'ERGCON'
  510 FORMAT(//1X,'THE PROGRAM CANNOT HANDLE THE FOLLOWING EXCITATION',
     1            ' INDICES.'/1X,8I5//)
      END
