      SUBROUTINE WSTORE (W,LM6,MODE)
C     *
C     COMPLETE DEFINITION OF SQUARE MATRIX OF MNDO TWO-ELECTRON
C     INTEGRALS BY INCLUDING THE ONE-CENTER TERMS AND THE TERMS
C     WITH TRANSPOSED INDICES.
C     *
C     MODE=-1     INCLUDE RAW ONE-CENTER INTEGRALS.
C     MODE= 0     INCLUDE RHF ONE-CENTER INTEGRALS AND TRANSPOSE.
C     MODE= 1     INCLUDE RAW ONE-CENTER INTEGRALS AND TRANSPOSE (UHF).
C     *
      USE LIMIT, ONLY: LM1, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL UHF
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./DPARM4/ REPD(52,LMZ)
     ./DPARM5/ INTIJ(243),INTKL(243),INTREP(243),INTRF1(243),INTRF2(243)
     ./INDEXW/ NW(LM1)
     ./REP   / GSS(LMZ),GPP(LMZ),GSP(LMZ),GP2(LMZ),HSP(LMZ),HPP(LMZ)
     ./UHF   / UHF
      DIMENSION W(LM6,LM6)
C *** INITIALIZE ONE-CENTER INTEGRALS (LOWER TRIANGLE) TO ZERO.
      MODW   = MODE
      IF(MODE.EQ.0 .AND. UHF) MODW=1
      IF(MODW.GE.0) THEN
         DO 20 II=1,NUMAT
         IORBS  = NLAST(II)-NFIRST(II)+1
         IF(IORBS.GE.4) THEN
            IW  = IORBS*(IORBS+1)/2
            IP  = NW(II)
            IPM = IP+IW-1
            DO 15 J=IP,IPM-1
            DO 10 I=J+1,IPM
            W(I,J) = ZERO
   10       CONTINUE
   15       CONTINUE
         ENDIF
   20    CONTINUE
      ENDIF
C *** INCLUDE NONZERO ONE-CENTER TERMS.
      IF(MODW.NE.0) THEN
CDIR$ NOBL
         DO 40 II=1,NUMAT
         IP  = NW(II)
         NI  = NAT(II)
         W(IP,IP) = GSS(NI)
         IORBS    = NLAST(II)-NFIRST(II)+1
         IF(IORBS.GE.4) THEN
            IPX = IP+2
            IPY = IP+5
            IPZ = IP+9
            W(IPX ,IP  ) = GSP(NI)
            W(IPY ,IP  ) = GSP(NI)
            W(IPZ ,IP  ) = GSP(NI)
            W(IP  ,IPX ) = GSP(NI)
            W(IP  ,IPY ) = GSP(NI)
            W(IP  ,IPZ ) = GSP(NI)
            W(IPX ,IPX ) = GPP(NI)
            W(IPY ,IPY ) = GPP(NI)
            W(IPZ ,IPZ ) = GPP(NI)
            W(IPY ,IPX ) = GP2(NI)
            W(IPZ ,IPX ) = GP2(NI)
            W(IPZ ,IPY ) = GP2(NI)
            W(IPX ,IPY ) = GP2(NI)
            W(IPX ,IPZ ) = GP2(NI)
            W(IPY ,IPZ ) = GP2(NI)
            W(IP+1,IP+1) = HSP(NI)
            W(IP+3,IP+3) = HSP(NI)
            W(IP+6,IP+6) = HSP(NI)
            W(IP+4,IP+4) = HPP(NI)
            W(IP+7,IP+7) = HPP(NI)
            W(IP+8,IP+8) = HPP(NI)
            IF(IORBS.GE.9) THEN
               IJ0    = IP-1
               DO 30 I=1,243
               IJ     = INTIJ(I)
               KL     = INTKL(I)
               INT    = INTREP(I)
               W(IJ+IJ0,KL+IJ0) = REPD(INT,NI)
   30          CONTINUE
            ENDIF
         ENDIF
   40    CONTINUE
      ELSE
         DO 60 II=1,NUMAT
         IP  = NW(II)
         NI  = NAT(II)
         W(IP,IP) = GSS(NI)*PT5
         IORBS    = NLAST(II)-NFIRST(II)+1
         IF(IORBS.GE.4) THEN
            IPX = IP+2
            IPY = IP+5
            IPZ = IP+9
            GSPNI = GSP(NI)-HSP(NI)*PT5
            GPPNI = GPP(NI)*PT5
            GP2NI = GP2(NI)-HPP(NI)*PT5
            HSPNI = HSP(NI)*0.75D0-GSP(NI)*PT25
            HPPNI = HPP(NI)*0.75D0-GP2(NI)*PT25
            W(IPX ,IP  ) = GSPNI
            W(IPY ,IP  ) = GSPNI
            W(IPZ ,IP  ) = GSPNI
            W(IP  ,IPX ) = GSPNI
            W(IP  ,IPY ) = GSPNI
            W(IP  ,IPZ ) = GSPNI
            W(IPX ,IPX ) = GPPNI
            W(IPY ,IPY ) = GPPNI
            W(IPZ ,IPZ ) = GPPNI
            W(IPY ,IPX ) = GP2NI
            W(IPZ ,IPX ) = GP2NI
            W(IPZ ,IPY ) = GP2NI
            W(IPX ,IPY ) = GP2NI
            W(IPX ,IPZ ) = GP2NI
            W(IPY ,IPZ ) = GP2NI
            W(IP+1,IP+1) = HSPNI
            W(IP+3,IP+3) = HSPNI
            W(IP+6,IP+6) = HSPNI
            W(IP+4,IP+4) = HPPNI
            W(IP+7,IP+7) = HPPNI
            W(IP+8,IP+8) = HPPNI
            IF(IORBS.GE.9) THEN
               IJ0    = IP-1
               DO 50 I=1,243
               IJ     = INTIJ(I)
               KL     = INTKL(I)
               INT    = INTREP(I)
               INT1   = INTRF1(I)
               INT2   = INTRF2(I)
               RF     = REPD(INT,NI)
               IF(INT1.GT.0) RF = RF-PT25*REPD(INT1,NI)
               IF(INT2.GT.0) RF = RF-PT25*REPD(INT2,NI)
               W(IJ+IJ0,KL+IJ0) = RF
   50          CONTINUE
            ENDIF
         ENDIF
   60    CONTINUE
      ENDIF
C *** INCLUDE TERMS WITH TRANSPOSED INDICES.
      IF(MODW.GE.0) THEN
         DO 80 I=2,LM6
CDIR$ IVDEP
C$DIR NO_RECURRENCE
*VDIR NODEP
*VOCL LOOP,NOVREC
         DO 70 J=1,I-1
         W(J,I) = W(I,J)
   70    CONTINUE
   80    CONTINUE
      ENDIF
      RETURN
      END
