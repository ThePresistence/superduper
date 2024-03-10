      SUBROUTINE PARAM (IOP,IPAROK)
C     *
C     PUT MNDO, AM1, PM3, MNDOC, OM1, OM2, OM3, OM4, ODM2 OR ODM3 PARAMETERS
C     INTO THE WORKING ARRAYS.
C     *
      USE LIMIT, ONLY: LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (DIPFAC=5.0832D0)
      COMMON /INOPT2/ IN2(300)
C *** COMMON BLOCKS USED IN THE ACTUAL CALCULATIONS.
      COMMON
     ./AMPGAU/ GUESS1(LMZ,4),GUESS2(LMZ,4),GUESS3(LMZ,4),GSCAL(LMZ),
     .         IMP(LMZ)
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./DIPOL1/ HYF(LMZ),HYFPD(LMZ)
     ./MULTIP/ DD(6,0:LMZ),PO(9,0:LMZ)
     ./PARAVL/ IMPAR(LMZ)
     ./PARDER/ CORE(LMZ),EHEAT(LMZ),EISOL(LMZ)
     ./PAROMM/ BETPI(LMZ,16)
     ./PAROPT/ USS(LMZ),UPP(LMZ),ZS(LMZ),ZP(LMZ),
     .         BETAS(LMZ),BETAP(LMZ),ALP(LMZ)
     ./PDDG1 / DA(LMZ,2),PA(LMZ,2)
     ./REP   / GSS(LMZ),GPP(LMZ),GSP(LMZ),GP2(LMZ),HSP(LMZ),
     .         HPP(LMZ)
     ./MMDPP / C6(LMZ),R0(LMZ),D3S6,D3S8,D3A1,D3A2
C *** COMMON BLOCK FOR EXPERIMENTAL HEATS OF FORMATION OF ATOMS.
      COMMON
     ./EXPHAT/ EXHEAT(LMZ)
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
     ./AM1DBL/ AM1C6(LMZ),AM1R0(LMZ),USSAMD(LMZ),UPPAMD(LMZ),BSAMD(LMZ),
     .         BPAMD(LMZ),ALPAMD(LMZ)
C *** COMMON BLOCKS FOR PM3
      COMMON
     ./PM3   / USSPM3(LMZ),UPPPM3(LMZ),ZSPM3(LMZ),ZPPM3(LMZ),
     .         BETASP(LMZ),BETAPP(LMZ),ALPPM3(LMZ),EISOLP(LMZ),
     .         DDPM3(LMZ),QQPM3(LMZ),AMPM3(LMZ),ADPM3(LMZ),AQPM3(LMZ),
     .         GSSPM3(LMZ),GSPPM3(LMZ),GPPPM3(LMZ),GP2PM3(LMZ),
     .         HSPPM3(LMZ)
     ./PM3GAU/ GUESP1(LMZ,2),GUESP2(LMZ,2),GUESP3(LMZ,2),IM3(LMZ)
     ./PM3DBL/ PM3C6(LMZ),PM3R0(LMZ),USSPMD(LMZ),UPPPMD(LMZ),BSPMD(LMZ),
     .         BPPMD(LMZ),ALPPMD(LMZ)
C *** COMMON BLOCKS FOR OM1
      COMMON
     ./OM1BL / USS1(LMZ),UPP1(LMZ),ZS1(LMZ),ZP1(LMZ),BETAS1(LMZ),
     .         BETAP1(LMZ),HYFOM1(LMZ),BETPI1(LMZ,10),
     .         DDOM1(LMZ),QQOM1(LMZ),AMOM1(LMZ),ADOM1(LMZ),AQOM1(LMZ),
     .         GSSOM1(LMZ),GSPOM1(LMZ),GPPOM1(LMZ),GP2OM1(LMZ),
     .         HSPOM1(LMZ),EISOL1(LMZ),IMOM1(LMZ)
C *** COMMON BLOCKS FOR OM2
      COMMON
     ./OM2BL / USS2(LMZ),UPP2(LMZ),ZS2(LMZ),ZP2(LMZ),BETAS2(LMZ),
     .         BETAP2(LMZ),HYFOM2(LMZ),BETPI2(LMZ,12),ZSCOR2(LMZ,4),
     .         DDOM2(LMZ),QQOM2(LMZ),AMOM2(LMZ),ADOM2(LMZ),AQOM2(LMZ),
     .         GSSOM2(LMZ),GSPOM2(LMZ),GPPOM2(LMZ),GP2OM2(LMZ),
     .         HSPOM2(LMZ),EISOL2(LMZ),IMOM2(LMZ)
     ./OM2DBL/ C62(LMZ,2),R02(LMZ,2)
C *** COMMON BLOCKS FOR OM3
      COMMON
     ./OM3BL/  USS3(LMZ),UPP3(LMZ),ZS3(LMZ),ZP3(LMZ),BETAS3(LMZ),
     .         BETAP3(LMZ),HYFOM3(LMZ),BETPI3(LMZ,12),ZSCOR3(LMZ,4),
     .         DDOM3(LMZ),QQOM3(LMZ),AMOM3(LMZ),ADOM3(LMZ),AQOM3(LMZ),
     .         GSSOM3(LMZ),GSPOM3(LMZ),GPPOM3(LMZ),GP2OM3(LMZ),
     .         HSPOM3(LMZ),EISOL3(LMZ),IMOM3(LMZ)
     ./OM3DBL/ C63(LMZ,2),R03(LMZ,2)
C *** COMMON BLOCKS FOR OM4
      COMMON
     ./OM4BL / USS4(LMZ),UPP4(LMZ),ZS4(LMZ),ZP4(LMZ),BETAS4(LMZ),
     .         BETAP4(LMZ),HYFOM4(LMZ),BETPI4(LMZ,8),ZSCOR4(LMZ,4),
     .         DDOM4(LMZ),QQOM4(LMZ),AMOM4(LMZ),ADOM4(LMZ),AQOM4(LMZ),
     .         GSSOM4(LMZ),GSPOM4(LMZ),GPPOM4(LMZ),GP2OM4(LMZ),
     .         HSPOM4(LMZ),EISOL4(LMZ),IMOM4(LMZ)
C *** COMMON BLOCKS FOR ODM2
      COMMON
     ./ODM2BL/ USSOD2(LMZ),UPPOD2(LMZ),ZSOMD2(LMZ),
     .         ZPOMD2(LMZ),BETSD2(LMZ),BETPD2(LMZ),HYFOD2(LMZ),
     .         BTPID2(LMZ,12),ZSCOD2(LMZ,4),
     .         DDODM2(LMZ),QQODM2(LMZ),AMODM2(LMZ),ADODM2(LMZ),
     .         AQODM2(LMZ),GSSOD2(LMZ),GSPOD2(LMZ),GPPOD2(LMZ),
     .         GP2OD2(LMZ),HSPOD2(LMZ),EISLD2(LMZ),IMODM2(LMZ)
C *** COMMON BLOCKS FOR ODM3
      COMMON
     ./ODM3BL/ USSOD3(LMZ),UPPOD3(LMZ),ZSOMD3(LMZ),
     .         ZPOMD3(LMZ),BETSD3(LMZ),BETPD3(LMZ),HYFOD3(LMZ),
     .         BTPID3(LMZ,12),ZSCOD3(LMZ,4),
     .         DDODM3(LMZ),QQODM3(LMZ),AMODM3(LMZ),ADODM3(LMZ),
     .         AQODM3(LMZ),GSSOD3(LMZ),GSPOD3(LMZ),GPPOD3(LMZ),
     .         GP2OD3(LMZ),HSPOD3(LMZ),EISLD3(LMZ),IMODM3(LMZ)
C *** COMMON BLOCKS FOR PDDG/MNDO
      COMMON
     ./PDDGM / USSPDM(LMZ),UPPPDM(LMZ),ZSPDM(LMZ),ZPPDM(LMZ),
     .         BSPDM(LMZ),BPPDM(LMZ),ALPPDM(LMZ),EISPDM(LMZ),
     .         DDPDM(LMZ),QQPDM(LMZ),AMPDM(LMZ),ADPDM(LMZ),AQPDM(LMZ),
     .         GSSPDM(LMZ),GPPPDM(LMZ),GSPPDM(LMZ),GP2PDM(LMZ),
     .         HSPPDM(LMZ)
     ./PDDGM1/ PAPDM(LMZ,2),DAPDM(LMZ,2),IMPDM(LMZ)
C *** COMMON BLOCKS FOR PDDG/PM3
      COMMON
     ./PDDGP / USSPDP(LMZ),UPPPDP(LMZ),ZSPDP(LMZ),ZPPDP(LMZ),
     .         BSPDP(LMZ),BPPDP(LMZ),ALPPDP(LMZ),EISPDP(LMZ),
     .         DDPDP(LMZ),QQPDP(LMZ),AMPDP(LMZ),ADPDP(LMZ),AQPDP(LMZ),
     .         GSSPDP(LMZ),GPPPDP(LMZ),GSPPDP(LMZ),GP2PDP(LMZ),
     .         HSPPDP(LMZ)
     ./PDDGP1/ GUPDP1(LMZ,2),GUPDP2(LMZ,2),GUPDP3(LMZ,2),
     .         PAPDP(LMZ,2),DAPDP(LMZ,2),IMPDP(LMZ)
      IMMDP=IN2(30)
C *** PM3 PARAMETERS.
      IF(IOP.EQ.-7) THEN
         DO 20 I=1,LMZ
         IMPAR(I) = IM3(I)
         IF(IMPAR(I).EQ.0) GO TO 20
         IMX      = IM3(I)
         IMP(I)   = IM3(I)
         USS(I)   = USSPM3(I)
         UPP(I)   = UPPPM3(I)
         ZS(I)    = ZSPM3(I)
         ZP(I)    = ZPPM3(I)
         BETAS(I) = BETASP(I)
         BETAP(I) = BETAPP(I)
         ALP(I)   = ALPPM3(I)
         GSS(I)   = GSSPM3(I)
         GSP(I)   = GSPPM3(I)
         GPP(I)   = GPPPM3(I)
         GP2(I)   = GP2PM3(I)
         HSP(I)   = HSPPM3(I)
         HPP(I)   = PT5*(GPPPM3(I)-GP2PM3(I))
         DO 10 J=1,IMX
         GUESS1(I,J) = GUESP1(I,J)
         GUESS2(I,J) = GUESP2(I,J)
         GUESS3(I,J) = GUESP3(I,J)
   10    CONTINUE
         GSCAL(I) = 1.D0
         IF(IPAROK.NE.7) THEN
            DD(2,I)  = DDPM3(I)
            DD(3,I)  = QQPM3(I)
            PO(1,I)  = PT5/AMPM3(I)
            PO(2,I)  = PT5/ADPM3(I)
            PO(3,I)  = PT5/AQPM3(I)
            PO(7,I)  = PO(1,I)
            PO(9,I)  = PO(1,I)
            HYF(I)   = DIPFAC*DD(2,I)
            HYFPD(I) = ZERO
            EISOL(I) = EISOLP(I)
         ELSE
            CALL DDPOHY(I)
            PO(9,I)  = PO(1,I)
            EISOL(I) = EATOM(I,0,0,0)
         ENDIF
         IF(IMMDP.GT.0) THEN
           SELECT CASE (I)
           CASE (1, 6, 7, 8)
             USS(I)=USSPMD(I)
             UPP(I)=UPPPMD(I)
             BETAS(I)=BSPMD(I)
             BETAP(I)=BPPMD(I)
             ALP(I)=ALPPMD(I)
             C6(I)=PM3C6(I)
             R0(I)=PM3R0(I)
           CASE DEFAULT
             C6(I)=0.0D0
             R0(I)=1.0D0
           END SELECT
         ENDIF
   20    CONTINUE
      ENDIF
C *** AM1 PARAMETERS.
      IF(IOP.EQ.-2) THEN
         DO 40 I=1,LMZ
         IMPAR(I) = IM1(I)
         IF(IMPAR(I).EQ.0) GO TO 40
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
         DO 30 J=1,IMX
         GUESS1(I,J) = GUESA1(I,J)
         GUESS2(I,J) = GUESA2(I,J)
         GUESS3(I,J) = GUESA3(I,J)
   30    CONTINUE
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
         IF(IMMDP.GT.0) THEN
           SELECT CASE (I)
           CASE (1, 6, 7, 8)
             USS(I)=USSAMD(I)
             UPP(I)=UPPAMD(I)
             BETAS(I)=BSAMD(I)
             BETAP(I)=BPAMD(I)
             ALP(I)=ALPAMD(I)
             C6(I)=AM1C6(I)
             R0(I)=AM1R0(I)
           CASE DEFAULT
             C6(I)=0.0D0
             R0(I)=1.0D0
           END SELECT
         ENDIF
   40    CONTINUE
      ENDIF
C *** MNDO PARAMETERS.
      IF(IOP.EQ.0 .OR. IOP.EQ.-3 .OR. IOP.EQ.-4) THEN
         DO 50 I=1,LMZ
         IMPAR(I) = IMM(I)
         IF(IMPAR(I).EQ.0) GO TO 50
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
   50    CONTINUE
      ENDIF
C     OLD HPP VALUES FOR Be,N (NOT SYMMETRY-ADAPTED).
C     OLD HPP VALUES FOR Si,S,Cl ARE NOT LISTED SINCE ALMOST
C     ALL PUBLISHED CALCULATIONS USE SYMMETRY-ADAPTED VALUES.
C     THE OLD CODE IS KEPT FOR DOCUMENTATION (COMMENT).
C     THIS OPTION IS NO LONGER SUPPORTED.
C     IF(IOP.EQ.-4) THEN
C        HPP(4)   = 0.380D0
C        HPP(7)   = 0.700D0
C        EISOL(7) =-202.581201D0
C        PO(3,4)  = 1.29435131D0
C        PO(3,7)  = 0.61389465D0
C     ENDIF
C     OLD MNDO PARAMETERS FOR SILICON AND SULPHUR, SEE:
C     M.J.S.DEWAR, M.L.MCKEE, H.S.RZEPA, J.AM.CHEM.SOC. 100,3607(1978).
C     NOTE THAT THE SYMMETRY-ADAPTED HPP VALUES FOR SILICON AND
C     SULPHUR ARE USED PRESENTLY (SAME CONVENTION AS IN MOPAC)
C     WHICH EXPLAINS THE DEVIATIONS OF THE EISOL VALUES FROM THE
C     ORIGINAL PAPER (SEE ABOVE).
      IF((IOP.EQ.0 .OR. IOP.EQ.-4) .AND. IPAROK.EQ.4) THEN
         USS(14)  = -40.5682920D0
         UPP(14)  = -28.0891870D0
         ZS(14)   =   1.4353060D0
         ZP(14)   = ZP(85)
         BETAS(14)=  -4.2562180D0
         BETAP(14)= BETAP(85)
         ALP(14)  =   2.1961078D0
         EISOL(14)= -90.5399580D0
         DD(2,14) =   1.4078712D0
         DD(3,14) =   1.1658281D0
         PO(1,14) = PT5/0.3608967D0
         PO(2,14) = PT5/0.3441817D0
         PO(3,14) = PT5/0.3999442D0
         PO(7,14) = PO(1,85)
         PO(9,14) = PO(1,85)
         HYF(14)  = DD(2,85)*DIPFAC
         HYFPD(14)= ZERO
         USS(16)  = -75.2391520D0
         UPP(16)  = -57.8320130D0
         ZS(16)   =   2.6135910D0
         ZP(16)   =   2.0343930D0
         BETAS(16)= -11.1422310D0
         BETAP(16)= BETAP(86)
         ALP(16)  =   2.4916445D0
         EISOL(16)=-235.4413560D0
         DD(2,16) =   0.8231596D0
         DD(3,16) =   0.8225156D0
         PO(1,16) = PT5/0.4733554D0
         PO(2,16) = PT5/0.5889395D0
         PO(3,16) = PT5/0.5632724D0
         PO(7,16) = PO(1,86)
         PO(9,16) = PO(1,86)
         HYF(16)  = DD(2,86)*DIPFAC
         HYFPD(16)= ZERO
      ENDIF
C     THE DEPENDENT MNDO PARAMETERS FOR H,Li,Be,B,C,N,O,F ARE
C     DEFINED MORE PRECISELY THAN IN MOPAC(6.0).
C     OPTION IOP=-4 COMBINED WITH IPAROK=5.
      IF(IOP.EQ.-4 .AND. IPAROK.EQ.5) THEN
         ALP(1)   = 2.54413400D0
         DD(2,1)  = ZERO
         DD(2,3)  = 2.05497832D0
         DD(2,4)  = 1.43732508D0
         DD(2,5)  = 0.95790723D0
         DD(2,6)  = 0.80746618D0
         DD(2,7)  = 0.63990367D0
         DD(2,8)  = 0.53460239D0
         DD(2,9)  = 0.50671661D0
         DD(3,1)  = ZERO
         DD(3,3)  = 1.74370693D0
         DD(3,4)  = 1.21961120D0
         DD(3,5)  = 0.81281124D0
         DD(3,6)  = 0.68515777D0
         DD(3,7)  = 0.54297627D0
         DD(3,8)  = 0.45362517D0
         DD(3,9)  = 0.42996330D0
         PO(1,1)  = 1.05891968D0
         PO(2,1)  = 1.05891968D0
         PO(3,1)  = 1.05891968D0
         PO(7,1)  = PO(1,1)
         PO(9,1)  = PO(1,1)
         PO(1,3)  = 1.86369863D0
         PO(2,3)  = 2.20284359D0
         PO(3,3)  = 1.91235216D0
         PO(1,4)  = 1.51166645D0
         PO(2,4)  = 1.48980568D0
         PO(3,4)  = 1.29992593D0
         PO(1,5)  = 1.28470255D0
         PO(2,5)  = 1.01942409D0
         PO(3,5)  = 0.89976948D0
         PO(1,6)  = 1.11242846D0
         PO(2,6)  = 0.81307772D0
         PO(3,6)  = 0.74784277D0
         PO(1,7)  = 1.00110375D0
         PO(2,7)  = 0.63745891D0
         PO(3,7)  = 0.61527519D0
         PO(1,8)  = 0.88229572D0
         PO(2,8)  = 0.52123718D0
         PO(3,8)  = 0.52654116D0
         PO(1,9)  = 0.80407801D0
         PO(2,9)  = 0.46081672D0
         PO(3,9)  = 0.48338867D0
         DO 60 I=3,9
         PO(7,I)  = PO(1,I)
         PO(9,I)  = PO(1,I)
         HYF(I)   = DIPFAC*DD(2,I)
         HYFPD(I) = ZERO
   60    CONTINUE
      ENDIF
C     SPECIAL MNDO PARAMETERS FOR CARBON ATOMS IN FULLERENES,
C     SEE THEOR.CHIM.ACTA 92, 269 (1995), TABLE I, SET B.
C     OPTION IOP=-4 COMBINED WITH IPAROK=6.
      IF(IOP.EQ.-4 .AND. IPAROK.EQ.6) THEN
         USS(6)   = -52.253460D0
         UPP(6)   = -39.215668D0
         ZS(6)    =   1.806560D0
         ZP(6)    =   1.780856D0
         BETAS(6) = -18.977848D0
         BETAP(6) = - 7.849880D0
         ALP(6)   =   2.566805D0
         EISOL(6) = -120.46825600D0
         DD(2,6)  = 0.80458493D0
         DD(3,6)  = 0.68772819D0
         PO(1,6)  = 1.11242845D0
         PO(2,6)  = 0.81148576D0
         PO(3,6)  = 0.74974989D0
         PO(7,6)  = PO(1,6)
         PO(9,6)  = PO(1,6)
         HYF(6)   = DD(2,6)*DIPFAC
         HYFPD(6) = ZERO
      ENDIF
C *** MNDOC PARAMETERS.
      IF(IOP.EQ.-1) THEN
         IMPAR(1) = 1
         IMPAR(6) = 1
         IMPAR(7) = 1
         IMPAR(8) = 1
         USS(1)   = -12.113619D0
         USS(6)   = -51.856594D0
         USS(7)   = -70.348077D0
         USS(8)   = -99.400519D0
         UPP(1)   = ZERO
         UPP(6)   = -39.306788D0
         UPP(7)   = -57.350327D0
         UPP(8)   = -77.268163D0
         ZS(1)    = 1.359938D0
         ZS(6)    = 1.828090D0
         ZS(7)    = 2.295131D0
         ZS(8)    = 2.733991D0
         ZP(1)    = ZS(1)
         ZP(6)    = ZS(6)
         ZP(7)    = ZS(7)
         ZP(8)    = ZS(8)
         BETAS(1) = - 6.929584D0
         BETAS(6) = -15.235343D0
         BETAS(7) = -20.129705D0
         BETAS(8) = -32.722231D0
         BETAP(1) = ZERO
         BETAP(6) = -10.297988D0
         BETAP(7) = BETAS(7)
         BETAP(8) = BETAS(8)
         ALP(1)   = 2.592386D0
         ALP(6)   = 2.550031D0
         ALP(7)   = 2.857123D0
         ALP(8)   = 3.168422D0
         EISOL(1) = USS(1)
         EISOL(6) = -120.133645D0
         EISOL(7) = -199.932135D0
         EISOL(8) = -315.263690D0
         DD(2,1)  = ZERO
         DD(2,6)  = 0.78955395D0
         DD(2,7)  = 0.62888596D0
         DD(2,8)  = 0.52793724D0
         DD(3,1)  = ZERO
         DD(3,6)  = 0.66995874D0
         DD(3,7)  = 0.53362744D0
         DD(3,8)  = 0.44796961D0
         PO(1,1)  = 1.05891968D0
         PO(2,1)  = 1.05891968D0
         PO(3,1)  = 1.05891968D0
         PO(1,6)  = 1.11242846D0
         PO(2,6)  = 0.80312540D0
         PO(3,6)  = 0.73650635D0
         PO(1,7)  = 1.00110375D0
         PO(2,7)  = 0.63143234D0
         PO(3,7)  = 0.60797443D0
         PO(1,8)  = 0.88229572D0
         PO(2,8)  = 0.51769876D0
         PO(3,8)  = 0.52199472D0
         DO 70 K=5,8
         I        = K
         IF(K.EQ.5) I=1
         GSS(I)   = GSSM(I)
         GSP(I)   = GSPM(I)
         GPP(I)   = GPPM(I)
         GP2(I)   = GP2M(I)
         HSP(I)   = HSPM(I)
         HPP(I)   = PT5*(GPPM(I)-GP2M(I))
         PO(7,I)  = PO(1,I)
         PO(9,I)  = PO(1,I)
         HYF(I)   = DIPFAC*DD(2,I)
         HYFPD(I) = ZERO
   70    CONTINUE
      ENDIF
C *** OM1 PARAMETERS.
      IF(IOP.EQ.-5) THEN
         LMOM1  = 86
         DO 90 I=1,LMOM1
         IMPAR(I) = IMOM1(I)
         IF(IMPAR(I).EQ.0) GO TO 90
         USS(I)   = USS1(I)
         UPP(I)   = UPP1(I)
         ZS(I)    = ZS1(I)
         ZP(I)    = ZP1(I)
         BETAS(I) = BETAS1(I)
         BETAP(I) = BETAP1(I)
         ALP(I)   = ZERO
         DO 80 J=1,10
         BETPI(I,J) = BETPI1(I,J)
   80    CONTINUE
         EISOL(I) = EISOL1(I)
         DD(2,I)  = DDOM1(I)
         DD(3,I)  = QQOM1(I)
         PO(1,I)  = PT5/AMOM1(I)
         PO(2,I)  = PT5/ADOM1(I)
         PO(3,I)  = PT5/AQOM1(I)
         PO(7,I)  = PO(1,I)
         PO(9,I)  = PO(1,I)
         GSS(I)   = GSSOM1(I)
         GSP(I)   = GSPOM1(I)
         GPP(I)   = GPPOM1(I)
         GP2(I)   = GP2OM1(I)
         HSP(I)   = HSPOM1(I)
         HPP(I)   = PT5*(GPPOM1(I)-GP2OM1(I))
         HYF(I)   = HYFOM1(I)
         HYFPD(I) = ZERO
   90    CONTINUE
      ENDIF
C *** OM2 PARAMETERS.
      IF(IOP.EQ.-6) THEN
         LMOM2   = 86
         DO 120 I=1,LMOM2
         IMPAR(I) = IMOM2(I)
         IF(IMPAR(I).EQ.0) GO TO 120
         USS(I)   = USS2(I)
         UPP(I)   = UPP2(I)
         ZS(I)    = ZS2(I)
         ZP(I)    = ZP2(I)
         BETAS(I) = BETAS2(I)
         BETAP(I) = BETAP2(I)
         ALP(I)   = ZERO
         DO 100 J=1,12
         BETPI(I,J) = BETPI2(I,J)
  100    CONTINUE
         DO 110 J=13,16
         BETPI(I,J) = ZSCOR2(I,J-12)
  110    CONTINUE
         EISOL(I) = EISOL2(I)
         DD(2,I)  = DDOM2(I)
         DD(3,I)  = QQOM2(I)
         PO(1,I)  = PT5/AMOM2(I)
         PO(2,I)  = PT5/ADOM2(I)
         PO(3,I)  = PT5/AQOM2(I)
         PO(7,I)  = PO(1,I)
         PO(9,I)  = PO(1,I)
         GSS(I)   = GSSOM2(I)
         GSP(I)   = GSPOM2(I)
         GPP(I)   = GPPOM2(I)
         GP2(I)   = GP2OM2(I)
         HSP(I)   = HSPOM2(I)
         HPP(I)   = PT5*(GPPOM2(I)-GP2OM2(I))
         HYF(I)   = HYFOM2(I)
         HYFPD(I) = ZERO
         C6(I)    = 0.0D0
         R0(I)    = 0.0D0
         D3S6     = 1.000d0
         D3S8     = 0.531d0
         D3A1     = 0.690d0
         D3A2     = 3.446d0
         IF(IMMDP.EQ.1.OR.IMMDP.EQ.2) THEN
           C6(I)    = C62(I,IMMDP)
           R0(I)    = R02(I,IMMDP)
         ENDIF
  120    CONTINUE
      ENDIF
C *** OM3 PARAMETERS.
      IF(IOP.EQ.-8) THEN
         LMOM3   = 9
         DO 150 I=1,LMOM3
         IMPAR(I) = IMOM3(I)
         IF(IMPAR(I).EQ.0) GO TO 150
         USS(I)   = USS3(I)
         UPP(I)   = UPP3(I)
         ZS(I)    = ZS3(I)
         ZP(I)    = ZP3(I)
         BETAS(I) = BETAS3(I)
         BETAP(I) = BETAP3(I)
         ALP(I)   = ZERO
         DO 130 J=1,12
         BETPI(I,J) = BETPI3(I,J)
  130    CONTINUE
         DO 140 J=13,16
         BETPI(I,J) = ZSCOR3(I,J-12)
  140    CONTINUE
         EISOL(I) = EISOL3(I)
         DD(2,I)  = DDOM3(I)
         DD(3,I)  = QQOM3(I)
         PO(1,I)  = PT5/AMOM3(I)
         PO(2,I)  = PT5/ADOM3(I)
         PO(3,I)  = PT5/AQOM3(I)
         PO(7,I)  = PO(1,I)
         PO(9,I)  = PO(1,I)
         GSS(I)   = GSSOM3(I)
         GSP(I)   = GSPOM3(I)
         GPP(I)   = GPPOM3(I)
         GP2(I)   = GP2OM3(I)
         HSP(I)   = HSPOM3(I)
         HPP(I)   = PT5*(GPPOM3(I)-GP2OM3(I))
         HYF(I)   = HYFOM3(I)
         HYFPD(I) = ZERO
         C6(I)    = 0.0D0
         R0(I)    = 0.0D0
         D3S6     = 1.000d0
         D3S8     = 0.501d0
         D3A1     = 0.613d0
         D3A2     = 3.258d0
         IF(IMMDP.EQ.1.OR.IMMDP.EQ.2) THEN
           C6(I)    = C63(I,IMMDP)
           R0(I)    = R03(I,IMMDP)
         ENDIF
  150    CONTINUE
      ENDIF
C *** OM4 PARAMETERS.
      IF(IOP.EQ.-9) THEN
         LMOM4   = 9
         DO 180 I=1,LMOM4
         IMPAR(I) = IMOM4(I)
         IF(IMPAR(I).EQ.0) GO TO 180
         USS(I)   = USS4(I)
         UPP(I)   = UPP4(I)
         ZS(I)    = ZS4(I)
         ZP(I)    = ZP4(I)
         BETAS(I) = BETAS4(I)
         BETAP(I) = BETAP4(I)
         ALP(I)   = ZERO
         DO 160 J=1,8
         BETPI(I,J) = BETPI4(I,J)
  160    CONTINUE
         DO 170 J=13,16
         BETPI(I,J) = ZSCOR4(I,J-12)
  170    CONTINUE
         EISOL(I) = EISOL4(I)
         DD(2,I)  = DDOM4(I)
         DD(3,I)  = QQOM4(I)
         PO(1,I)  = PT5/AMOM4(I)
         PO(2,I)  = PT5/ADOM4(I)
         PO(3,I)  = PT5/AQOM4(I)
         PO(7,I)  = PO(1,I)
         PO(9,I)  = PO(1,I)
         GSS(I)   = GSSOM4(I)
         GSP(I)   = GSPOM4(I)
         GPP(I)   = GPPOM4(I)
         GP2(I)   = GP2OM4(I)
         HSP(I)   = HSPOM4(I)
         HPP(I)   = PT5*(GPPOM4(I)-GP2OM4(I))
         HYF(I)   = HYFOM4(I)
         HYFPD(I) = ZERO
  180    CONTINUE
      ENDIF
C *** ODM2 PARAMETERS.
      IF(IOP.EQ.-22) THEN
         LMODM2  = 86
         DO 210 I=1,LMODM2
         IMPAR(I) = IMODM2(I)
         IF(IMPAR(I).EQ.0) GO TO 210
         USS(I)   = USSOD2(I)
         UPP(I)   = UPPOD2(I)
         ZS(I)    = ZSOMD2(I)
         ZP(I)    = ZPOMD2(I)
         BETAS(I) = BETSD2(I)
         BETAP(I) = BETPD2(I)
         ALP(I)   = ZERO
         DO 190 J=1,12
         BETPI(I,J) = BTPID2(I,J)
  190    CONTINUE
         DO 200 J=13,16
         BETPI(I,J) = ZSCOD2(I,J-12)
  200    CONTINUE
         EISOL(I) = EISLD2(I)
         DD(2,I)  = DDODM2(I)
         DD(3,I)  = QQODM2(I)
         PO(1,I)  = PT5/AMODM2(I)
         PO(2,I)  = PT5/ADODM2(I)
         PO(3,I)  = PT5/AQODM2(I)
         PO(7,I)  = PO(1,I)
         PO(9,I)  = PO(1,I)
         GSS(I)   = GSSOD2(I)
         GSP(I)   = GSPOD2(I)
         GPP(I)   = GPPOD2(I)
         GP2(I)   = GP2OD2(I)
         HSP(I)   = HSPOD2(I)
         HPP(I)   = PT5*(GPPOD2(I)-GP2OD2(I))
         HYF(I)   = HYFOD2(I)
         HYFPD(I) = ZERO
         C6(I)    = 0.0D0
         R0(I)    = 0.0D0
         D3S6     = 1.000D0
         D3S8     = 0.56467295D0
         D3A1     = 0.66492067D0
         D3A2     = 3.34510640D0
  210    CONTINUE
      ENDIF
C *** ODM3 PARAMETERS.
      IF(IOP.EQ.-23) THEN
         LMODM3  = 86
         DO 240 I=1,LMODM3
         IMPAR(I) = IMODM3(I)
         IF(IMPAR(I).EQ.0) GO TO 240
         USS(I)   = USSOD3(I)
         UPP(I)   = UPPOD3(I)
         ZS(I)    = ZSOMD3(I)
         ZP(I)    = ZPOMD3(I)
         BETAS(I) = BETSD3(I)
         BETAP(I) = BETPD3(I)
         ALP(I)   = ZERO
         DO 220 J=1,12
         BETPI(I,J) = BTPID3(I,J)
  220    CONTINUE
         DO 230 J=13,16
         BETPI(I,J) = ZSCOD3(I,J-12)
  230    CONTINUE
         EISOL(I) = EISLD3(I)
         DD(2,I)  = DDODM3(I)
         DD(3,I)  = QQODM3(I)
         PO(1,I)  = PT5/AMODM3(I)
         PO(2,I)  = PT5/ADODM3(I)
         PO(3,I)  = PT5/AQODM3(I)
         PO(7,I)  = PO(1,I)
         PO(9,I)  = PO(1,I)
         GSS(I)   = GSSOD3(I)
         GSP(I)   = GSPOD3(I)
         GPP(I)   = GPPOD3(I)
         GP2(I)   = GP2OD3(I)
         HSP(I)   = HSPOD3(I)
         HPP(I)   = PT5*(GPPOD3(I)-GP2OD3(I))
         HYF(I)   = HYFOD3(I)
         HYFPD(I) = ZERO
         C6(I)    = 0.0D0
         R0(I)    = 0.0D0
         D3S6     = 1.000D0
         D3S8     = 0.78131768D0
         D3A1     = 0.69959980D0
         D3A2     = 3.05404380D0
  240    CONTINUE
      ENDIF
C *** PCCG/PM3 PARAMETERS.
      IF(IOP.EQ.-7 .AND. IPAROK.EQ.8) THEN
         DO 270 I=1,LMZ
         IMPAR(I) = IMPDP(I)
         IF(IMPAR(I).EQ.0) GO TO 270
         IMX      = IMPDP(I)
         IMP(I)   = IMPDP(I)
         USS(I)   = USSPDP(I)
         UPP(I)   = UPPPDP(I)
         ZS(I)    = ZSPDP(I)
         ZP(I)    = ZPPDP(I)
         BETAS(I) = BSPDP(I)
         BETAP(I) = BPPDP(I)
         ALP(I)   = ALPPDP(I)
         GSS(I)   = GSSPDP(I)
         GSP(I)   = GSPPDP(I)
         GPP(I)   = GPPPDP(I)
         GP2(I)   = GP2PDP(I)
         HSP(I)   = HSPPDP(I)
         HPP(I)   = PT5*(GPPPDP(I)-GP2PDP(I))
         DO 250 J=1,IMX
         GUESS1(I,J) = GUPDP1(I,J)
         GUESS2(I,J) = GUPDP2(I,J)
         GUESS3(I,J) = GUPDP3(I,J)
  250    CONTINUE
         GSCAL(I) = 1.D0
         DO 260 J=1,2
         DA(I,J)  = DAPDP(I,J)
         PA(I,J)  = PAPDP(I,J)
  260    CONTINUE
         DD(2,I)  = DDPDP(I)
         DD(3,I)  = QQPDP(I)
         PO(1,I)  = PT5/AMPDP(I)
         PO(2,I)  = PT5/ADPDP(I)
         PO(3,I)  = PT5/AQPDP(I)
         PO(7,I)  = PO(1,I)
         PO(9,I)  = PO(1,I)
         HYF(I)   = DIPFAC*DD(2,I)
         HYFPD(I) = ZERO
         EISOL(I) = EISPDP(I)
  270    CONTINUE
      ENDIF
C *** PCCG/MNDO PARAMETERS.
      IF((IOP.EQ.0 .OR. IOP.EQ.-4) .AND. IPAROK.EQ.8) THEN
         DO 290 I=1,LMZ
         IMPAR(I) = IMPDM(I)
         IF(IMPAR(I).EQ.0) GO TO 290
         IMX      = IMPDM(I)
         IMP(I)   = IMPDM(I)
         USS(I)   = USSPDM(I)
         UPP(I)   = UPPPDM(I)
         ZS(I)    = ZSPDM(I)
         ZP(I)    = ZPPDM(I)
         BETAS(I) = BSPDM(I)
         BETAP(I) = BPPDM(I)
         ALP(I)   = ALPPDM(I)
         GSS(I)   = GSSPDM(I)
         GSP(I)   = GSPPDM(I)
         GPP(I)   = GPPPDM(I)
         GP2(I)   = GP2PDM(I)
         HSP(I)   = HSPPDM(I)
         HPP(I)   = PT5*(GPPPDM(I)-GP2PDM(I))
         DO 280 J=1,2
         DA(I,J)  = DAPDM(I,J)
         PA(I,J)  = PAPDM(I,J)
  280    CONTINUE
         DD(2,I)  = DDPDM(I)
         DD(3,I)  = QQPDM(I)
         PO(1,I)  = PT5/AMPDM(I)
         PO(2,I)  = PT5/ADPDM(I)
         PO(3,I)  = PT5/AQPDM(I)
         PO(7,I)  = PO(1,I)
         PO(9,I)  = PO(1,I)
         HYF(I)   = DIPFAC*DD(2,I)
         HYFPD(I) = ZERO
         EISOL(I) = EISPDM(I)
  290    CONTINUE
      ENDIF
C *** SPECIAL CONVENTIONS.
C     ZP(1) = ZS(1)
C *** DEFINE EXPERIMENTAL HEATS OF FORMATION OF THE ATOMS.
      DO 300 I=1,LMZ
      EHEAT(I) = EXHEAT(I)
  300 CONTINUE
      RETURN
      END
