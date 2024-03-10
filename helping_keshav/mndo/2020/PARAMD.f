      SUBROUTINE PARAMD (IOP,IPAROK)
C     *
C     PUT MNDO/d OR AM1/d PARAMETERS INTO THE WORKING ARRAYS.
C     *
      USE LIMIT, ONLY: LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (DIPFAC=5.0832D0)
C *** COMMON BLOCKS USED IN THE ACTUAL CALCULATIONS.
      COMMON
     ./AMPGAU/ GUESS1(LMZ,4),GUESS2(LMZ,4),GUESS3(LMZ,4),GSCAL(LMZ),
     .         IMP(LMZ)
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./DIPOL1/ HYF(LMZ),HYFPD(LMZ)
     ./DPARM / LORBS(LMZ)
     ./DPARM1/ UDD(LMZ),ZD(LMZ),BETAD(LMZ)
     ./DPARM2/ ZSN(LMZ),ZPN(LMZ),ZDN(LMZ)
     ./DPARM3/ F0DD(LMZ),F2DD(LMZ),F4DD(LMZ),F0SD(LMZ),G2SD(LMZ),
     .         F0PD(LMZ),F2PD(LMZ),G1PD(LMZ),G3PD(LMZ)
     ./DPARM7/ IF0SD(LMZ),IG2SD(LMZ),IGSP(LMZ),IHSP(LMZ)
     ./DPARM8/ JF0SD(LMZ),JG2SD(LMZ),JGSP(LMZ),JHSP(LMZ)
     ./MULTIP/ DD(6,0:LMZ),PO(9,0:LMZ)
     ./PARAVL/ IMPAR(LMZ)
     ./PARDER/ CORE(LMZ),EHEAT(LMZ),EISOL(LMZ)
     ./PAROPT/ USS(LMZ),UPP(LMZ),ZS(LMZ),ZP(LMZ),
     .         BETAS(LMZ),BETAP(LMZ),ALP(LMZ)
     ./REP   / GSS(LMZ),GPP(LMZ),GSP(LMZ),GP2(LMZ),HSP(LMZ),
     .         HPP(LMZ)
C *** COMMON BLOCK FOR EXPERIMENTAL HEATS OF FORMATION OF ATOMS.
      COMMON
     ./EXPHAT/ EXHEAT(LMZ)
C *** COMMON BLOCKS FOR MNDO/d
      COMMON
     ./MNDOD / USSD(LMZ),UPPD(LMZ),ZSD(LMZ),ZPD(LMZ),
     .         BETASD(LMZ),BETAPD(LMZ),ALPD(LMZ),
     .         GSSD(LMZ),GPPD(LMZ),GSPD(LMZ),GP2D(LMZ),
     .         HSPD(LMZ),UDDD(LMZ),ZDD(LMZ),BETADD(LMZ),
     .         ZSND(LMZ),ZPND(LMZ),ZDND(LMZ),POCORD(LMZ),
     .         F0SDD(LMZ),G2SDD(LMZ),IMMD(LMZ)
C *** COMMON BLOCKS FOR AM1/d
      COMMON
     ./AM1D  / USSAMD(LMZ),UPPAMD(LMZ),ZSAMD(LMZ),ZPAMD(LMZ),
     .         BSAMD(LMZ),BPAMD(LMZ),ALPAMD(LMZ),
     .         GSSAMD(LMZ),GSPAMD(LMZ),GPPAMD(LMZ),GP2AMD(LMZ),
     .         HSPAMD(LMZ),UDDAMD(LMZ),ZDAMD(LMZ),BDAMD(LMZ),
     .         ZSNAMD(LMZ),ZPNAMD(LMZ),ZDNAMD(LMZ),POCAMD(LMZ),
     .         F0SDAM(LMZ),G2SDAM(LMZ),IMAMD(LMZ)
     ./AMDGAU/ GUESD1(LMZ,4),GUESD2(LMZ,4),GUESD3(LMZ,4),IM1D(LMZ)
C *** COMMON BLOCKS FOR MNDO
      COMMON
     ./MNDOBL/ USSM(LMZ),UPPM(LMZ),ZSM(LMZ),ZPM(LMZ),
     .         BETASM(LMZ),BETAPM(LMZ),ALPM(LMZ),EISOLM(LMZ),
     .         DDM(LMZ),QQM(LMZ),AMM(LMZ),ADM(LMZ),AQM(LMZ),
     .         GSSM(LMZ),GPPM(LMZ),GSPM(LMZ),GP2M(LMZ),HSPM(LMZ)
     ./MNDOEL/ IMM(LMZ)
C *** COMMON BLOCKS FOR AM1
      COMMON
     ./AM1   / USSAM1(LMZ),UPPAM1(LMZ),ZSAM1(LMZ),ZPAM1(LMZ),
     .         BETASA(LMZ),BETAPA(LMZ),ALPAM1(LMZ),EISOLA(LMZ),
     .         DDAM1(LMZ),QQAM1(LMZ),AMAM1(LMZ),ADAM1(LMZ),AQAM1(LMZ),
     .         GSSAM1(LMZ),GSPAM1(LMZ),GPPAM1(LMZ),GP2AM1(LMZ),
     .         HSPAM1(LMZ)
     ./AM1GAU/ GUESA1(LMZ,4),GUESA2(LMZ,4),GUESA3(LMZ,4),IM1(LMZ)
C     *
C *** AM1/d PARAMETERS.
C     *
C     FIRST SECTION.
C     USE STANDARD AM1 PARAMETERS BY DEFAULT.
      IF(IOP.EQ.-12) THEN
         DO 10 I=1,LMZ
         IMPAR(I) = IM1(I)
         IF(IMPAR(I).EQ.0) GO TO 10
         IMX      = IM1(I)
         IMP(I)   = IM1(I)
         USS(I)   = USSAM1(I)
         UPP(I)   = UPPAM1(I)
         ZS(I)    = ZSAM1(I)
         ZP(I)    = ZPAM1(I)
         BETAS(I) = BETASA(I)
         BETAP(I) = BETAPA(I)
         ALP(I)   = ALPAM1(I)
         GSS(I)   = GSSAM1(I)
         GSP(I)   = GSPAM1(I)
         GPP(I)   = GPPAM1(I)
         GP2(I)   = GP2AM1(I)
         HSP(I)   = HSPAM1(I)
         HPP(I)   = PT5*(GPPAM1(I)-GP2AM1(I))
         DO 11 J=1,IMX
         GUESS1(I,J) = GUESA1(I,J)
         GUESS2(I,J) = GUESA2(I,J)
         GUESS3(I,J) = GUESA3(I,J)
   11    CONTINUE
         GSCAL(I) = 1.D0
         IF(IPAROK.NE.7) THEN
            DD(2,I)  = DDAM1(I)
            DD(3,I)  = QQAM1(I)
            PO(1,I)  = PT5/AMAM1(I)
            PO(2,I)  = PT5/ADAM1(I)
            PO(3,I)  = PT5/AQAM1(I)
            PO(7,I)  = PO(1,I)
            PO(9,I)  = PO(1,I)
            HYF(I)   = DIPFAC*DD(2,I)
            HYFPD(I) = ZERO
            EISOL(I) = EISOLA(I)
         ELSE
            CALL DDPOHY(I)
            PO(9,I)  = PO(1,I)
            EISOL(I) = EATOM(I,0,0,0)
         ENDIF
   10    CONTINUE
C     SECOND SECTION.
C     USE AVAILABLE AM1/d PARAMETERS.
C     IMAMD(I)=2: DEFINE A FULL SET OF SP  PARAMETERS FOR I.GT.10.
C     IMAMD(I)=3: DEFINE A FULL SET OF SPD PARAMETERS FOR I.GT.10.
         DO 20 I=11,LMZ
         IF(IMAMD(I).EQ.0) GO TO 20
         IMPAR(I) = IMAMD(I)
         IMX      = IM1D(I)
         IMP(I)   = IM1D(I)
         USS(I)   = USSAMD(I)
         UPP(I)   = UPPAMD(I)
         ZS(I)    = ZSAMD(I)
         ZP(I)    = ZPAMD(I)
         BETAS(I) = BSAMD(I)
         BETAP(I) = BPAMD(I)
         ALP(I)   = ALPAMD(I)
         PO(9,I)  = POCAMD(I)
         IF(IMAMD(I).EQ.2) THEN
            GSS(I)   = GSSAMD(I)
            GSP(I)   = GSPAMD(I)
            GPP(I)   = GPPAMD(I)
            GP2(I)   = GP2AMD(I)
            HSP(I)   = HSPAMD(I)
            HPP(I)   = PT5*(GPPAMD(I)-GP2AMD(I))
         ELSE IF(IMAMD(I).EQ.3) THEN
            UDD(I)   = UDDAMD(I)
            ZD(I)    = ZDAMD(I)
            BETAD(I) = BDAMD(I)
            ZSN(I)   = ZSNAMD(I)
            ZPN(I)   = ZPNAMD(I)
            ZDN(I)   = ZDNAMD(I)
            IF(IF0SD(I).GT.0) F0SD(I) = F0SDAM(I)
            IF(IG2SD(I).GT.0) G2SD(I) = G2SDAM(I)
            IF(IGSP(I) .GT.0) GSP(I)  = GSPAMD(I)
            IF(IHSP(I) .GT.0) HSP(I)  = HSPAMD(I)
            JF0SD(I) = IF0SD(I)
            JG2SD(I) = IG2SD(I)
            JGSP(I)  = IGSP(I)
            JHSP(I)  = IHSP(I)
            LORBS(I) = 9
            CALL INIGHD(I)
         ENDIF
         CALL DDPOHY(I)
         EISOL(I) = EATOM(I,0,0,0)
         DO 21 J=1,IMX
         GUESS1(I,J) = GUESD1(I,J)
         GUESS2(I,J) = GUESD2(I,J)
         GUESS3(I,J) = GUESD3(I,J)
   21    CONTINUE
         GSCAL(I) = 1.D0
   20    CONTINUE
      ENDIF
C     *
C *** MNDO/d PARAMETERS.
C     *
C     FIRST SECTION.
C     USE STANDARD MNDO PARAMETERS FOR I.LE.10 BY DEFAULT.
C     USE STANDARD MNDO PARAMETERS FOR I.EQ.86 BY DEFAULT (LINK ATOM).
      IF(IOP.EQ.-10 .OR. IOP.EQ.-13) THEN
         DO 30 K=1,11
         I = K
         IF(K.EQ.11) I=86
         IMPAR(I) = IMM(I)
         IF(IMPAR(I).EQ.0) GO TO 30
         USS(I)   = USSM(I)
         UPP(I)   = UPPM(I)
         ZS(I)    = ZSM(I)
         ZP(I)    = ZPM(I)
         BETAS(I) = BETASM(I)
         BETAP(I) = BETAPM(I)
         ALP(I)   = ALPM(I)
         GSS(I)   = GSSM(I)
         GSP(I)   = GSPM(I)
         GPP(I)   = GPPM(I)
         GP2(I)   = GP2M(I)
         HSP(I)   = HSPM(I)
         HPP(I)   = PT5*(GPPM(I)-GP2M(I))
         IF(IPAROK.NE.7) THEN
            DD(2,I)  = DDM(I)
            DD(3,I)  = QQM(I)
            PO(1,I)  = PT5/AMM(I)
            PO(2,I)  = PT5/ADM(I)
            PO(3,I)  = PT5/AQM(I)
            PO(7,I)  = PO(1,I)
            PO(9,I)  = PO(1,I)
            HYF(I)   = DIPFAC*DD(2,I)
            HYFPD(I) = ZERO
            EISOL(I) = EISOLM(I)
         ELSE
            CALL DDPOHY(I)
            PO(9,I)  = PO(1,I)
            EISOL(I) = EATOM(I,0,0,0)
         ENDIF
   30    CONTINUE
C     SECOND SECTION.
C     IMMD(I)=2, NI.LE.10: DEFINE A FULL SET OF NEW SP PARAMETERS.
         DO 40 I=3,10
         IF(IMMD(I).EQ.2) THEN
            IMPAR(I) = IMMD(I)
            USS(I)   = USSD(I)
            ZS(I)    = ZSD(I)
            ALP(I)   = ALPD(I)
            BETAS(I) = BETASD(I)
            PO(9,I)  = POCORD(I)
            GSS(I)   = GSSD(I)
            UPP(I)   = UPPD(I)
            ZP(I)    = ZPD(I)
            BETAP(I) = BETAPD(I)
            GSP(I)   = GSPD(I)
            GPP(I)   = GPPD(I)
            GP2(I)   = GP2D(I)
            HSP(I)   = HSPD(I)
            HPP(I)   = PT5*(GPPD(I)-GP2D(I))
            CALL DDPOHY(I)
            EISOL(I) = EATOM(I,0,0,0)
         ENDIF
   40    CONTINUE
C     THIRD SECTION.
C     IMMD(I)=2, NI.GT.10: DEFINE A FULL SET OF NEW SP  PARAMETERS.
C     IMMD(I)=3, NI.GT.10: DEFINE A FULL SET OF NEW SPD PARAMETERS.
         DO 50 I=11,LMZ
         IMPAR(I) = IMMD(I)
         IF(IMMD(I).LE.1) GO TO 50
         USS(I)   = USSD(I)
         UPP(I)   = UPPD(I)
         ZS(I)    = ZSD(I)
         ZP(I)    = ZPD(I)
         BETAS(I) = BETASD(I)
         BETAP(I) = BETAPD(I)
         ALP(I)   = ALPD(I)
         PO(9,I)  = POCORD(I)
         IF(IMMD(I).EQ.2) THEN
            GSS(I)   = GSSD(I)
            GSP(I)   = GSPD(I)
            GPP(I)   = GPPD(I)
            GP2(I)   = GP2D(I)
            HSP(I)   = HSPD(I)
            HPP(I)   = PT5*(GPPD(I)-GP2D(I))
         ELSE IF(IMMD(I).EQ.3) THEN
            UDD(I)   = UDDD(I)
            ZD(I)    = ZDD(I)
            BETAD(I) = BETADD(I)
            ZSN(I)   = ZSND(I)
            ZPN(I)   = ZPND(I)
            ZDN(I)   = ZDND(I)
            IF(IF0SD(I).GT.0) F0SD(I) = F0SDD(I)
            IF(IG2SD(I).GT.0) G2SD(I) = G2SDD(I)
            IF(IGSP(I) .GT.0) GSP(I)  = GSPD(I)
            IF(IHSP(I) .GT.0) HSP(I)  = HSPD(I)
            JF0SD(I) = IF0SD(I)
            JG2SD(I) = IG2SD(I)
            JGSP(I)  = IGSP(I)
            JHSP(I)  = IHSP(I)
            LORBS(I) = 9
            CALL INIGHD(I)
         ENDIF
         CALL DDPOHY(I)
         EISOL(I) = EATOM(I,0,0,0)
   50    CONTINUE
      ENDIF
C *** DEFINE EXPERIMENTAL HEATS OF FORMATION OF THE ATOMS.
      DO 200 I=1,LMZ
      EHEAT(I) = EXHEAT(I)
  200 CONTINUE
      RETURN
      END
