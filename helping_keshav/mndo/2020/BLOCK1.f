      BLOCKDATA BLOCK1
C     *
C     INITIALIZATION FOR PARAMETERS IN MNDO, AM1, AND PM3.
C     SOME OF THE DATA STATEMENTS HAVE BEEN COPIED FROM MOPAC(6.0)
C     WRITTEN BY J.J.P.STEWART.
C     *
C     COMMON BLOCKS: PARDER, EXPHAT, MNDOBL, MNDOEL, AM1, AM1GAU, PM3,
C     PM3GAU.
C     *
C     SPECIAL CONVENTION.
C     ELEMENT 86 IS A CONNECTION ATOM FOR QM/MM TREATMENTS.
C     CORE CHARGE 1, ONE ELECTRON IN A 2S ORBITAL.
C     EXPERIMENTAL HEAT OF FORMATION AND MASS FROM METHYL (CH3).
C     *
      USE LIMIT, ONLY: LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     CHARACTER REFMN*80,REFM3*80,REFAM*80,REFPM3*80
C *** GENERAL COMMON BLOCKS
      COMMON
     ./PARDER/ CORE(LMZ),EHEAT(LMZ),EISOL(LMZ)
C    ./REFS  / REFMN(LMZ),REFM3(LMZ),REFAM(LMZ),REFPM3(LMZ)
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
C *** CORE CHARGES.
C     LINES      ATOMIC NUMBERS
C     1 - 3          1 - 18
C     4 - 5         19 - 36
C     6 - 7         37 - 54
C     8 - B         55 - 86
      DATA CORE/
     1  1.D0,2.D0,
     2  1.D0,2.D0,3.D0,4.D0,5.D0,6.D0,7.D0,8.D0,
     3  1.D0,2.D0,3.D0,4.D0,5.D0,6.D0,7.D0,8.D0,
     4  1.D0,2.D0,3.D0,4.D0,5.D0,6.D0,7.D0,8.D0,9.D0,10.D0,11.D0,2.D0,
     5            3.D0,4.D0,5.D0,6.D0,7.D0,8.D0,
     6  1.D0,2.D0,3.D0,4.D0,5.D0,6.D0,7.D0,8.D0,9.D0,10.D0,11.D0,2.D0,
     7            3.D0,4.D0,5.D0,6.D0,7.D0,8.D0,
     8  1.D0,2.D0,3.D0,3.D0,3.D0,3.D0,3.D0,3.D0,3.D0,
     9            3.D0,3.D0,3.D0,3.D0,3.D0,3.D0,3.D0,
     A            3.D0,4.D0,5.D0,6.D0,7.D0,8.D0,9.D0,10.D0,11.D0,2.D0,
     B            3.D0,4.D0,5.D0,6.D0,7.D0,1.D0/
C     *
C     ENTHALPIES OF FORMATION OF GASEOUS ATOMS ARE TAKEN FROM  ANNUAL
C     REPORTS,1974,71B,P 117.  THERE ARE SOME SIGNIFICANT DIFFERENCES
C     BETWEEN THE VALUES REPORTED THERE AND THE VALUES PREVIOUSLY IN
C     THE BLOCK DATA OF THIS PROGRAM.  ONLY THE THIRD  ROW ELEMENTS
C     HAVE BEEN UPDATED.
C     *
C     ALL THE OTHER ELEMENTS ARE TAKEN FROM CRC HANDBOOK 1981-1982.
      DATA EXHEAT(1)  / 52.102D0/
      DATA EXHEAT(2)  /  0.000D0/
C
      DATA EXHEAT(3)  / 38.410D0/
      DATA EXHEAT(4)  / 76.960D0/
      DATA EXHEAT(5)  /135.700D0/
      DATA EXHEAT(6)  /170.890D0/
      DATA EXHEAT(7)  /113.000D0/
      DATA EXHEAT(8)  / 59.559D0/
      DATA EXHEAT(9)  / 18.890D0/
      DATA EXHEAT(10) /  0.000D0/
C
      DATA EXHEAT(11) / 25.650D0/
      DATA EXHEAT(12) / 35.000D0/
      DATA EXHEAT(13) / 79.490D0/
      DATA EXHEAT(14) /108.390D0/
      DATA EXHEAT(15) / 75.570D0/
      DATA EXHEAT(16) / 66.400D0/
      DATA EXHEAT(17) / 28.990D0/
      DATA EXHEAT(18) /  0.000D0/
C
      DATA EXHEAT(19) / 21.420D0/
      DATA EXHEAT(20) / 42.600D0/
      DATA EXHEAT(21) / 90.300D0/
      DATA EXHEAT(22) /112.300D0/
      DATA EXHEAT(23) /122.900D0/
      DATA EXHEAT(24) / 95.000D0/
      DATA EXHEAT(25) / 67.700D0/
      DATA EXHEAT(26) / 99.300D0/
      DATA EXHEAT(27) /102.400D0/
      DATA EXHEAT(28) /102.800D0/
      DATA EXHEAT(29) / 80.700D0/
      DATA EXHEAT(30) / 31.170D0/
      DATA EXHEAT(31) / 65.400D0/
      DATA EXHEAT(32) / 89.500D0/
      DATA EXHEAT(33) / 72.300D0/
      DATA EXHEAT(34) / 54.300D0/
      DATA EXHEAT(35) / 26.740D0/
      DATA EXHEAT(36) /  0.000D0/
C
      DATA EXHEAT(37) / 19.600D0/
      DATA EXHEAT(38) / 39.100D0/
      DATA EXHEAT(39) /101.500D0/
      DATA EXHEAT(40) /145.500D0/
      DATA EXHEAT(41) /172.400D0/
      DATA EXHEAT(42) /157.300D0/
      DATA EXHEAT(43) /  0.000D0/
      DATA EXHEAT(44) /155.500D0/
      DATA EXHEAT(45) /133.000D0/
      DATA EXHEAT(46) / 90.000D0/
      DATA EXHEAT(47) / 68.100D0/
      DATA EXHEAT(48) / 26.720D0/
      DATA EXHEAT(49) / 58.000D0/
      DATA EXHEAT(50) / 72.200D0/
      DATA EXHEAT(51) / 63.200D0/
      DATA EXHEAT(52) / 47.000D0/
      DATA EXHEAT(53) / 25.517D0/
      DATA EXHEAT(54) /  0.000D0/
C
      DATA EXHEAT(55) / 18.700D0/
      DATA EXHEAT(56) / 42.500D0/
      DATA EXHEAT(57) /  0.000D0/
      DATA EXHEAT(58) /101.300D0/
      DATA EXHEAT(59) /  0.000D0/
      DATA EXHEAT(60) /  0.000D0/
      DATA EXHEAT(61) /  0.000D0/
      DATA EXHEAT(62) / 49.400D0/
      DATA EXHEAT(63) /  0.000D0/
      DATA EXHEAT(64) /  0.000D0/
      DATA EXHEAT(65) /  0.000D0/
      DATA EXHEAT(66) /  0.000D0/
      DATA EXHEAT(67) /  0.000D0/
      DATA EXHEAT(68) / 75.800D0/
      DATA EXHEAT(69) /  0.000D0/
      DATA EXHEAT(70) / 36.350D0/
      DATA EXHEAT(71) /  0.000D0/
      DATA EXHEAT(72) /148.000D0/
      DATA EXHEAT(73) /186.900D0/
      DATA EXHEAT(74) /203.100D0/
      DATA EXHEAT(75) /185.000D0/
      DATA EXHEAT(76) /188.000D0/
      DATA EXHEAT(77) /160.000D0/
      DATA EXHEAT(78) /135.200D0/
      DATA EXHEAT(79) / 88.000D0/
      DATA EXHEAT(80) / 14.690D0/
      DATA EXHEAT(81) / 43.550D0/
      DATA EXHEAT(82) / 46.620D0/
      DATA EXHEAT(83) / 50.100D0/
      DATA EXHEAT(84) /  0.000D0/
      DATA EXHEAT(85) /  0.000D0/
      DATA EXHEAT(86) / 34.800D0/
C     *
C *** START OF MNDO PARAMETERS.
C     *
C     LIST OF ELEMENTS WITH MNDO PARAMETERS.
      DATA IMM/1,1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,0,
     1         0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,1,0,
     2         0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,
     3         25*0,                 1,0,1,0,0,0,1/
C     DATA FOR ELEMENT  1        HYDROGEN
C     DATA REFMN  ( 1)/'  H: (MNDO):  M.J.S. DEWAR, W. THIEL, J. AM. CHE
C    1M. SOC., 99, 4899 (1977)        '/
      DATA USSM   ( 1)/     -11.9062760D0/
      DATA UPPM   ( 1)/       0.0000000D0/
      DATA BETASM ( 1)/      -6.9890640D0/
      DATA BETAPM ( 1)/      -6.9890640D0/
      DATA ZSM    ( 1)/       1.3319670D0/
      DATA ZPM    ( 1)/       1.3319670D0/
      DATA ALPM   ( 1)/       2.5441341D0/
      DATA EISOLM ( 1)/     -11.9062760D0/
      DATA GSSM   ( 1)/      12.8480000D0/
      DATA GSPM   ( 1)/       0.0000000D0/
      DATA GPPM   ( 1)/       0.0000000D0/
      DATA GP2M   ( 1)/       0.0000000D0/
      DATA HSPM   ( 1)/       0.0000000D0/
      DATA DDM    ( 1)/       0.0000000D0/
      DATA QQM    ( 1)/       0.0000000D0/
      DATA AMM    ( 1)/       0.4721793D0/
      DATA ADM    ( 1)/       0.4721793D0/
      DATA AQM    ( 1)/       0.4721793D0/
C     DATA FOR ELEMENT  2        HELIUM
C     DATA REFMN  ( 2)/'  He: (MNDO):  M. KOLB, W. THIEL, J. COMP. CHEM.
C    1, 14, 37, (1993)                '/
      DATA USSM   ( 2)/     -54.4030000D0/
      DATA UPPM   ( 2)/       0.0000000D0/
      DATA BETASM ( 2)/     -28.6938694D0/
      DATA BETAPM ( 2)/     -28.6938694D0/
      DATA ZSM    ( 2)/       2.0179064D0/
      DATA ZPM    ( 2)/       2.0179064D0/
      DATA ALPM   ( 2)/       3.2620769D0/
      DATA EISOLM ( 2)/     -78.9830000D0/
      DATA GSSM   ( 2)/      29.8230000D0/
      DATA GSPM   ( 2)/       0.0000000D0/
      DATA GPPM   ( 2)/       0.0000000D0/
      DATA GP2M   ( 2)/       0.0000000D0/
      DATA HSPM   ( 2)/       0.0000000D0/
      DATA DDM    ( 2)/       0.0000000D0/
      DATA QQM    ( 2)/       0.0000000D0/
      DATA AMM    ( 2)/       1.0960309D0/
      DATA ADM    ( 2)/       1.0960309D0/
      DATA AQM    ( 2)/       1.0960309D0/
C     DATA FOR ELEMENT  3        LITHIUM
C     DATA REFMN  ( 3)/' Li: (MNDO):  TAKEN FROM MNDOC BY W.THIEL,
C    1QCPE NO.438, V. 2, P.63, (1982).'/
      DATA USSM   ( 3)/      -5.1280000D0/
      DATA UPPM   ( 3)/      -2.7212000D0/
      DATA BETASM ( 3)/      -1.3500400D0/
      DATA BETAPM ( 3)/      -1.3500400D0/
      DATA ZSM    ( 3)/       0.7023800D0/
      DATA ZPM    ( 3)/       0.7023800D0/
      DATA ALPM   ( 3)/       1.2501400D0/
      DATA EISOLM ( 3)/      -5.1280000D0/
      DATA GSSM   ( 3)/       7.3000000D0/
      DATA GSPM   ( 3)/       5.4200000D0/
      DATA GPPM   ( 3)/       5.0000000D0/
      DATA GP2M   ( 3)/       4.5200000D0/
      DATA HSPM   ( 3)/       0.8300000D0/
      DATA DDM    ( 3)/       2.0549783D0/
      DATA QQM    ( 3)/       1.7437069D0/
      DATA AMM    ( 3)/       0.2682837D0/
      DATA ADM    ( 3)/       0.2269793D0/
      DATA AQM    ( 3)/       0.2614581D0/
C     DATA FOR ELEMENT  4        BERYLLIUM
C     DATA REFMN  ( 4)/' Be: (MNDO):  M.J.S. DEWAR, H.S. RZEPA, J. AM. C
C    1HEM. SOC., 100, 777, (1978)     '/
      DATA USSM   ( 4)/     -16.6023780D0/
      DATA UPPM   ( 4)/     -10.7037710D0/
      DATA BETASM ( 4)/      -4.0170960D0/
      DATA BETAPM ( 4)/      -4.0170960D0/
      DATA ZSM    ( 4)/       1.0042100D0/
      DATA ZPM    ( 4)/       1.0042100D0/
      DATA ALPM   ( 4)/       1.6694340D0/
      DATA EISOLM ( 4)/     -24.2047560D0/
      DATA GSSM   ( 4)/       9.0000000D0/
      DATA GSPM   ( 4)/       7.4300000D0/
      DATA GPPM   ( 4)/       6.9700000D0/
      DATA GP2M   ( 4)/       6.2200000D0/
      DATA HSPM   ( 4)/       1.2800000D0/
      DATA DDM    ( 4)/       1.4373245D0/
      DATA QQM    ( 4)/       1.2196103D0/
      DATA AMM    ( 4)/       0.3307607D0/
      DATA ADM    ( 4)/       0.3356142D0/
      DATA AQM    ( 4)/       0.3846373D0/
C     DATA FOR ELEMENT  5        BORON
C     DATA REFMN  ( 5)/'  B: (MNDO):  M.J.S. DEWAR, M.L. MCKEE, J. AM. C
C    1HEM. SOC., 99, 5231 (1977)      '/
      DATA USSM   ( 5)/     -34.5471300D0/
      DATA UPPM   ( 5)/     -23.1216900D0/
      DATA BETASM ( 5)/      -8.2520540D0/
      DATA BETAPM ( 5)/      -8.2520540D0/
      DATA ZSM    ( 5)/       1.5068010D0/
      DATA ZPM    ( 5)/       1.5068010D0/
      DATA ALPM   ( 5)/       2.1349930D0/
      DATA EISOLM ( 5)/     -64.3159500D0/
      DATA GSSM   ( 5)/      10.5900000D0/
      DATA GSPM   ( 5)/       9.5600000D0/
      DATA GPPM   ( 5)/       8.8600000D0/
      DATA GP2M   ( 5)/       7.8600000D0/
      DATA HSPM   ( 5)/       1.8100000D0/
      DATA DDM    ( 5)/       0.9579073D0/
      DATA QQM    ( 5)/       0.8128113D0/
      DATA AMM    ( 5)/       0.3891951D0/
      DATA ADM    ( 5)/       0.4904730D0/
      DATA AQM    ( 5)/       0.5556979D0/
C     DATA FOR ELEMENT  6        CARBON
C     DATA REFMN  ( 6)/'  C: (MNDO):  M.J.S. DEWAR, W. THIEL, J. AM. CHE
C    1M. SOC., 99, 4899 (1977)        '/
      DATA USSM   ( 6)/     -52.2797450D0/
      DATA UPPM   ( 6)/     -39.2055580D0/
      DATA BETASM ( 6)/     -18.9850440D0/
      DATA BETAPM ( 6)/      -7.9341220D0/
      DATA ZSM    ( 6)/       1.7875370D0/
      DATA ZPM    ( 6)/       1.7875370D0/
      DATA ALPM   ( 6)/       2.5463800D0/
      DATA EISOLM ( 6)/    -120.5006060D0/
      DATA GSSM   ( 6)/      12.2300000D0/
      DATA GSPM   ( 6)/      11.4700000D0/
      DATA GPPM   ( 6)/      11.0800000D0/
      DATA GP2M   ( 6)/       9.8400000D0/
      DATA HSPM   ( 6)/       2.4300000D0/
      DATA DDM    ( 6)/       0.8074662D0/
      DATA QQM    ( 6)/       0.6851578D0/
      DATA AMM    ( 6)/       0.4494671D0/
      DATA ADM    ( 6)/       0.6149474D0/
      DATA AQM    ( 6)/       0.6685897D0/
C     DATA FOR ELEMENT  7        NITROGEN
C     DATA REFMN  ( 7)/'  N: (MNDO):  M.J.S. DEWAR, W. THIEL, J. AM. CHE
C    1M. SOC., 99, 4899 (1977)        '/
C     CORRECTION OF EISOLM AND AQM VALUES ACCORDING TO
C     J.J.P. STEWART, QCPE BULL. 9, 80 (1989).
      DATA USSM   ( 7)/     -71.9321220D0/
      DATA UPPM   ( 7)/     -57.1723190D0/
      DATA BETASM ( 7)/     -20.4957580D0/
      DATA BETAPM ( 7)/     -20.4957580D0/
      DATA ZSM    ( 7)/       2.2556140D0/
      DATA ZPM    ( 7)/       2.2556140D0/
      DATA ALPM   ( 7)/       2.8613420D0/
      DATA EISOLM ( 7)/    -202.5662010D0/
      DATA GSSM   ( 7)/      13.5900000D0/
      DATA GSPM   ( 7)/      12.6600000D0/
      DATA GPPM   ( 7)/      12.9800000D0/
      DATA GP2M   ( 7)/      11.5900000D0/
      DATA HSPM   ( 7)/       3.1400000D0/
      DATA DDM    ( 7)/       0.6399037D0/
      DATA QQM    ( 7)/       0.5429763D0/
      DATA AMM    ( 7)/       0.4994487D0/
      DATA ADM    ( 7)/       0.7843643D0/
      DATA AQM    ( 7)/       0.8126445D0/
C     DATA FOR ELEMENT  8        OXYGEN
C     DATA REFMN  ( 8)/'  O: (MNDO):  M.J.S. DEWAR, W. THIEL, J. AM. CHE
C    1M. SOC., 99, 4899 (1977)        '/
      DATA USSM   ( 8)/     -99.6443090D0/
      DATA UPPM   ( 8)/     -77.7974720D0/
      DATA BETASM ( 8)/     -32.6880820D0/
      DATA BETAPM ( 8)/     -32.6880820D0/
      DATA ZSM    ( 8)/       2.6999050D0/
      DATA ZPM    ( 8)/       2.6999050D0/
      DATA ALPM   ( 8)/       3.1606040D0/
      DATA EISOLM ( 8)/    -317.8685060D0/
      DATA GSSM   ( 8)/      15.4200000D0/
      DATA GSPM   ( 8)/      14.4800000D0/
      DATA GPPM   ( 8)/      14.5200000D0/
      DATA GP2M   ( 8)/      12.9800000D0/
      DATA HSPM   ( 8)/       3.9400000D0/
      DATA DDM    ( 8)/       0.5346024D0/
      DATA QQM    ( 8)/       0.4536252D0/
      DATA AMM    ( 8)/       0.5667034D0/
      DATA ADM    ( 8)/       0.9592562D0/
      DATA AQM    ( 8)/       0.9495934D0/
C     DATA FOR ELEMENT  9        FLUORINE
C     DATA REFMN  ( 9)/'  F: (MNDO):  M.J.S. DEWAR, H.S. RZEPA, J. AM. C
C    1HEM. SOC., 100, 777 (1978)      '/
      DATA USSM   ( 9)/    -131.0715480D0/
      DATA UPPM   ( 9)/    -105.7821370D0/
      DATA BETASM ( 9)/     -48.2904660D0/
      DATA BETAPM ( 9)/     -36.5085400D0/
      DATA ZSM    ( 9)/       2.8484870D0/
      DATA ZPM    ( 9)/       2.8484870D0/
      DATA ALPM   ( 9)/       3.4196606D0/
      DATA EISOLM ( 9)/    -476.6837810D0/
      DATA GSSM   ( 9)/      16.9200000D0/
      DATA GSPM   ( 9)/      17.2500000D0/
      DATA GPPM   ( 9)/      16.7100000D0/
      DATA GP2M   ( 9)/      14.9100000D0/
      DATA HSPM   ( 9)/       4.8300000D0/
      DATA DDM    ( 9)/       0.5067166D0/
      DATA QQM    ( 9)/       0.4299633D0/
      DATA AMM    ( 9)/       0.6218302D0/
      DATA ADM    ( 9)/       1.0850301D0/
      DATA AQM    ( 9)/       1.0343643D0/
C     DATA FOR THE SODIUM-LIKE SPARKLE
C     DATA REFMN  (11)/' Na: (MNDO):  SODIUM-LIKE SPARKLE.   USE WITH CA
C    1RE.                             '/
C     DATA EHEAT  (11)/       0.0000000D0/
C     DATA ALPAM1 (11)/       1.6680000D0/
C     DATA EISOLA (11)/       0.0000000D0/
C     DATA AMAM1  (11)/       0.5000000D0/
C     DATA ALPM   (11)/       1.6600000D0/
C     DATA EISOLM (11)/       0.0000000D0/
C     DATA AMM    (11)/       0.5000000D0/
C     DATA FOR ELEMENT 12        MAGNESIUM
C     DATA REFMN  (12)/' Mg: (MNDO):  A.A.VOITYUK, ZHURNAL STRUKTURNOI K
C    1HIMII 28, 128 (1987)            '/
      DATA USSM   (12)/     -15.0400000D0/
      DATA UPPM   (12)/     - 9.2640000D0/
      DATA BETASM (12)/      -2.5860000D0/
      DATA BETAPM (12)/      -2.8420000D0/
      DATA ZSM    (12)/       1.0490000D0/
      DATA ZPM    (12)/       0.8890000D0/
C     DATA ZDM    (12)/       1.0000000D0/
      DATA ALPM   (12)/       1.8130000D0/
      DATA EISOLM (12)/     -22.6900000D0/
      DATA GSSM   (12)/       7.3900000D0/
      DATA GSPM   (12)/       6.5700000D0/
      DATA GPPM   (12)/       6.6800000D0/
      DATA GP2M   (12)/       5.9000000D0/
      DATA HSPM   (12)/       0.8200000D0/
      DATA DDM    (12)/       2.0360000D0/
      DATA QQM    (12)/       1.8820000D0/
      DATA AMM    (12)/       0.2715915D0/
      DATA ADM    (12)/       0.2269632D0/
      DATA AQM    (12)/       0.2925688D0/
C     DATA FOR ELEMENT 13        ALUMINUM
C     DATA REFMN  (13)/' Al: (MNDO):  L.P. DAVIS, ET.AL.  J. COMP. CHEM.
C    1, 2, 433 (1981) SEE MANUAL.     '/
C     THE MONOCENTRIC INTEGRALS HSP AND GSP FOR ALUMINIUM ARE ONLY
C     ESTIMATES. A VALUE OF G1 FOR AL IS NEEDED TO RESOLVE OLEARIS
C     INTEGRALS.
      DATA USSM   (13)/     -23.8070970D0/
      DATA UPPM   (13)/     -17.5198780D0/
      DATA BETASM (13)/      -2.6702840D0/
      DATA BETAPM (13)/      -2.6702840D0/
      DATA ZSM    (13)/       1.4441610D0/
      DATA ZPM    (13)/       1.4441610D0/
C     DATA ZDM    (13)/       1.0000000D0/
      DATA ALPM   (13)/       1.8688394D0/
      DATA EISOLM (13)/     -44.4840720D0/
      DATA GSSM   (13)/       8.0900000D0/
      DATA GSPM   (13)/       6.6300000D0/
      DATA GPPM   (13)/       5.9800000D0/
      DATA GP2M   (13)/       5.4000000D0/
      DATA HSPM   (13)/       0.7000000D0/
      DATA DDM    (13)/       1.3992387D0/
      DATA QQM    (13)/       1.1586797D0/
      DATA AMM    (13)/       0.2973172D0/
      DATA ADM    (13)/       0.2635574D0/
      DATA AQM    (13)/       0.3673560D0/
C     DATA FOR ELEMENT 14          SILICON
C     DATA REFMN  (14)/' Si: (MNDO): M.J.S.DEWAR, ET. AL. ORGANOMETALLIC
C    1S  5, 375 (1986)                '/
      DATA USSM   (14)/     -37.0375330D0/
      DATA UPPM   (14)/     -27.7696780D0/
      DATA BETASM (14)/      -9.0868040D0/
      DATA BETAPM (14)/      -1.0758270D0/
      DATA ZSM    (14)/       1.3159860D0/
      DATA ZPM    (14)/       1.7099430D0/
C     DATA ZDM    (14)/       1.0000000D0/
      DATA ALPM   (14)/       2.2053160D0/
      DATA EISOLM (14)/     -82.8394220D0/
      DATA GSSM   (14)/       9.8200000D0/
      DATA GSPM   (14)/       8.3600000D0/
      DATA GPPM   (14)/       7.3100000D0/
      DATA GP2M   (14)/       6.5400000D0/
      DATA HSPM   (14)/       1.3200000D0/
      DATA DDM    (14)/       1.2580349D0/
      DATA QQM    (14)/       0.9785824D0/
      DATA AMM    (14)/       0.3608967D0/
      DATA ADM    (14)/       0.3664244D0/
      DATA AQM    (14)/       0.4506740D0/
C     DATA FOR ELEMENT 15        PHOSPHORUS
C     DATA REFMN  (15)/'  P: (MNDO): M.J.S.DEWAR, M.L.MCKEE, H.S.RZEPA,
C    1J. AM. CHEM. SOC., 100 3607 1978'/
      DATA USSM   (15)/     -56.1433600D0/
      DATA UPPM   (15)/     -42.8510800D0/
      DATA BETASM (15)/      -6.7916000D0/
      DATA BETAPM (15)/      -6.7916000D0/
      DATA ZSM    (15)/       2.1087200D0/
      DATA ZPM    (15)/       1.7858100D0/
C     DATA ZDM    (15)/       1.0000000D0/
      DATA ALPM   (15)/       2.4152800D0/
      DATA EISOLM (15)/    -152.9599600D0/
      DATA GSSM   (15)/      11.5600000D0/
      DATA GSPM   (15)/      10.0800000D0/
      DATA GPPM   (15)/       8.6400000D0/
      DATA GP2M   (15)/       7.6800000D0/
      DATA HSPM   (15)/       1.9200000D0/
      DATA DDM    (15)/       1.0129699D0/
      DATA QQM    (15)/       0.9370090D0/
      DATA AMM    (15)/       0.4248438D0/
      DATA ADM    (15)/       0.4882420D0/
      DATA AQM    (15)/       0.4979406D0/
C     DATA FOR ELEMENT 16        SULFUR
C     DATA REFMN  (16)/'  S: (MNDO): M.J.S.DEWAR, C.H. REYNOLDS, J. COM
C    1P. CHEM. 7, 140-143 (1986)      '/
      DATA USSM   (16)/     -72.2422810D0/
      DATA UPPM   (16)/     -56.9732070D0/
      DATA BETASM (16)/     -10.7616700D0/
      DATA BETAPM (16)/     -10.1084330D0/
      DATA ZSM    (16)/       2.3129620D0/
      DATA ZPM    (16)/       2.0091460D0/
C     DATA ZDM    (16)/       1.0000000D0/
      DATA ALPM   (16)/       2.4780260D0/
      DATA EISOLM (16)/    -226.0123900D0/
      DATA GSSM   (16)/      12.8800000D0/
      DATA GSPM   (16)/      11.2600000D0/
      DATA GPPM   (16)/       9.9000000D0/
      DATA GP2M   (16)/       8.8300000D0/
      DATA HSPM   (16)/       2.2600000D0/
      DATA DDM    (16)/       0.9189935D0/
      DATA QQM    (16)/       0.8328514D0/
      DATA AMM    (16)/       0.4733554D0/
      DATA ADM    (16)/       0.5544502D0/
      DATA AQM    (16)/       0.5585244D0/
C     DATA FOR ELEMENT 17        CHLORINE
C     DATA REFMN  (17)/' Cl: (MNDO): M.J.S.DEWAR, H.S.RZEPA, J. COMP. CH
C    1EM., 4, 158 (1983)              '/
      DATA USSM   (17)/    -100.2271660D0/
      DATA UPPM   (17)/     -77.3786670D0/
      DATA BETASM (17)/     -14.2623200D0/
      DATA BETAPM (17)/     -14.2623200D0/
      DATA ZSM    (17)/       3.7846450D0/
      DATA ZPM    (17)/       2.0362630D0/
C     DATA ZDM    (17)/       1.0000000D0/
      DATA ALPM   (17)/       2.5422010D0/
      DATA EISOLM (17)/    -353.1176670D0/
      DATA GSSM   (17)/      15.0300000D0/
      DATA GSPM   (17)/      13.1600000D0/
      DATA GPPM   (17)/      11.3000000D0/
      DATA GP2M   (17)/       9.9700000D0/
      DATA HSPM   (17)/       2.4200000D0/
      DATA DDM    (17)/       0.4986870D0/
      DATA QQM    (17)/       0.8217603D0/
      DATA AMM    (17)/       0.5523705D0/
      DATA ADM    (17)/       0.8061220D0/
      DATA AQM    (17)/       0.6053435D0/
C     DATA FOR THE POTASSIUM-LIKE SPARKLE
C     DATA REFAM  (19)/' K:  (AM1):  POTASSIUM-LIKE SPARKLE.   USE WITH
C    1 CARE.                          '/
C     DATA EHEAT  (19)/       0.0000000D0/
C     DATA ALPAM1 (19)/       1.4050000D0/
C     DATA EISOLA (19)/       0.0000000D0/
C     DATA AMAM1  (19)/       0.5000000D0/
C     DATA ALPM   (19)/       1.3960000D0/
C     DATA EISOLM (19)/       0.0000000D0/
C     DATA AMM    (19)/       0.5000000D0/
C     DATA FOR ELEMENT 30        ZINC
C     DATA REFMN  (30)/' Zn: (MNDO):  M.J.S. DEWAR, K.M. MERZ, ORGANOMET
C    1ALLICS, 5, 1494-1496 (1986)     '/
      DATA USSM  ( 30)/     -20.8397160D0/
      DATA UPPM  ( 30)/     -19.6252240D0/
      DATA BETASM( 30)/      -1.0000000D0/
      DATA BETAPM( 30)/      -2.0000000D0/
      DATA ZSM   ( 30)/       2.0473590D0/
      DATA ZPM   ( 30)/       1.4609460D0/
C     DATA ZDM   ( 30)/       1.0000000D0/
      DATA ALPM  ( 30)/       1.5064570D0/
      DATA EISOLM( 30)/     -29.8794320D0/
      DATA GSSM  ( 30)/      11.8000000D0/
      DATA GSPM  ( 30)/      11.1820180D0/
      DATA GPPM  ( 30)/      13.3000000D0/
      DATA GP2M  ( 30)/      12.9305200D0/
      DATA HSPM  ( 30)/       0.4846060D0/
      DATA DDM   ( 30)/       1.3037826D0/
      DATA QQM   ( 30)/       1.4520183D0/
      DATA AMM   ( 30)/       0.4336641D0/
      DATA ADM   ( 30)/       0.2375912D0/
      DATA AQM   ( 30)/       0.2738858D0/
C     DATA FOR ELEMENT 32        GERMANIUM
C     DATA REFMN  (32)/' Ge: (MNDO): M.J.S.DEWAR, G.L.GRADY, E.F.HEALY,O
C    1RGANOMETALLICS 6, 186-189 (1987)'/
      DATA USSM  ( 32)/     -33.9493670D0/
      DATA UPPM  ( 32)/     -27.4251050D0/
      DATA BETASM( 32)/      -4.5164790D0/
      DATA BETAPM( 32)/      -1.7555170D0/
      DATA ZSM   ( 32)/       1.2931800D0/
      DATA ZPM   ( 32)/       2.0205640D0/
      DATA ALPM  ( 32)/       1.9784980D0/
      DATA EISOLM( 32)/     -76.2489440D0/
      DATA GSSM  ( 32)/       9.8000000D0/
      DATA GSPM  ( 32)/       8.3000000D0/
      DATA GPPM  ( 32)/       7.3000000D0/
      DATA GP2M  ( 32)/       6.5000000D0/
      DATA HSPM  ( 32)/       1.3000000D0/
      DATA DDM   ( 32)/       1.2556091D0/
      DATA QQM   ( 32)/       1.0498655D0/
      DATA AMM   ( 32)/       0.3601617D0/
      DATA ADM   ( 32)/       0.3643722D0/
      DATA AQM   ( 32)/       0.4347337D0/
C     DATA FOR ELEMENT 35        BROMINE
C     DATA REFMN  (35)/' Br: (MNDO): M.J.S.DEWAR, E.F. HEALY, J. COMP. C
C    1HEM., 4, 542 (1983)             '/
      DATA USSM   (35)/     -99.9864405D0/
      DATA UPPM   (35)/     -75.6713075D0/
      DATA BETASM (35)/      -8.9171070D0/
      DATA BETAPM (35)/      -9.9437400D0/
      DATA ZSM    (35)/       3.8543019D0/
      DATA ZPM    (35)/       2.1992091D0/
C     DATA ZDM    (35)/       1.0000000D0/
      DATA ALPM   (35)/       2.4457051D0/
      DATA EISOLM (35)/    -346.6812500D0/
      DATA GSSM   (35)/      15.03643948D0/
      DATA GSPM   (35)/      13.03468242D0/
      DATA GPPM   (35)/      11.27632539D0/
      DATA GP2M   (35)/       9.85442552D0/
      DATA HSPM   (35)/       2.45586832D0/
      DATA DDM    (35)/       0.6051074D0/
      DATA QQM    (35)/       0.9645873D0/
      DATA AMM    (35)/       0.5526068D0/
      DATA ADM    (35)/       0.7258330D0/
      DATA AQM    (35)/       0.5574589D0/
C     DATA FOR ELEMENT 50        TIN
C     DATA REFMN  (50)/' Sn: (MNDO): M.J.S.DEWAR,G.L.GRADY,J.J.P.STEWART
C    1, J.AM.CHEM.SOC.,106,6771 (1984)'/
      DATA USSM   (50)/     -40.8518020D0/
      DATA UPPM   (50)/     -28.5602490D0/
      DATA BETASM (50)/      -3.2351470D0/
      DATA BETAPM (50)/      -4.2904160D0/
      DATA ZSM    (50)/       2.0803800D0/
      DATA ZPM    (50)/       1.9371060D0/
      DATA ALPM   (50)/       1.8008140D0/
      DATA EISOLM (50)/     -92.3241020D0/
      DATA GSSM   (50)/       9.8000000D0/
      DATA GSPM   (50)/       8.3000000D0/
      DATA GPPM   (50)/       7.3000000D0/
      DATA GP2M   (50)/       6.5000000D0/
      DATA HSPM   (50)/       1.3000000D0/
      DATA DDM    (50)/       1.5697766D0/
      DATA QQM    (50)/       1.3262292D0/
      DATA AMM    (50)/       0.3601617D0/
      DATA ADM    (50)/       0.3219998D0/
      DATA AQM    (50)/       0.3713827D0/
C     DATA FOR ELEMENT 53        IODINE
C     DATA REFMN  (53)/'  I: (MNDO): M.J.S.DEWAR, E.F. HEALY, J.J.P. STE
C    1WART, J.COMP.CHEM., 5,358 (1984)'/
      DATA USSM   (53)/    -100.0030538D0/
      DATA UPPM   (53)/     -74.6114692D0/
      DATA BETASM (53)/      -7.4144510D0/
      DATA BETAPM (53)/      -6.1967810D0/
      DATA ZSM    (53)/       2.2729610D0/
      DATA ZPM    (53)/       2.1694980D0/
C     DATA ZDM    (53)/       1.0000000D0/
      DATA ALPM   (53)/       2.2073200D0/
      DATA EISOLM (53)/    -340.5983600D0/
      DATA GSSM   (53)/      15.04044855D0/
      DATA GSPM   (53)/      13.05655798D0/
      DATA GPPM   (53)/      11.14778369D0/
      DATA GP2M   (53)/       9.91409071D0/
      DATA HSPM   (53)/       2.45638202D0/
      DATA DDM    (53)/       1.4253233D0/
      DATA QQM    (53)/       1.1841707D0/
      DATA AMM    (53)/       0.5527541D0/
      DATA ADM    (53)/       0.4593451D0/
      DATA AQM    (53)/       0.4585376D0/
C     DATA FOR ELEMENT 80        MERCURY
C     DATA REFMN  (80)/' Hg: (MNDO): M.J.S.DEWAR,  ET. AL. ORGANOMETALLI
C    1CS 4, 1964 (1985) SEE MANUAL    '/
      DATA USSM   (80)/     -19.8095740D0/
      DATA UPPM   (80)/     -13.1025300D0/
      DATA BETASM (80)/      -0.4045250D0/
      DATA BETAPM (80)/      -6.2066830D0/
      DATA ZSM    (80)/       2.2181840D0/
      DATA ZPM    (80)/       2.0650380D0/
      DATA ALPM   (80)/       1.3356410D0/
      DATA EISOLM (80)/     -28.8191480D0/
      DATA GSSM   (80)/      10.8000000D0/
      DATA GSPM   (80)/       9.3000000D0/
      DATA GPPM   (80)/      14.3000000D0/
      DATA GP2M   (80)/      13.5000000D0/
      DATA HSPM   (80)/       1.3000000D0/
      DATA DDM    (80)/       1.7378048D0/
      DATA QQM    (80)/       1.4608064D0/
      DATA AMM    (80)/       0.3969129D0/
      DATA ADM    (80)/       0.3047694D0/
      DATA AQM    (80)/       0.3483102D0/
C     DATA FOR ELEMENT 82        LEAD
C     DATA REFMN  (82)/' Pb: (MNDO): M.J.S.DEWAR, ET.AL ORGANOMETALLICS
C    14, 1973-1980 (1985)             '/
      DATA USSM   (82)/     -47.3196920D0/
      DATA UPPM   (82)/     -28.8475600D0/
      DATA BETASM (82)/      -8.0423870D0/
      DATA BETAPM (82)/      -3.0000000D0/
      DATA ZSM    (82)/       2.4982860D0/
      DATA ZPM    (82)/       2.0820710D0/
      DATA ALPM   (82)/       1.7283330D0/
      DATA EISOLM (82)/    -105.8345040D0/
      DATA GSSM   (82)/       9.8000000D0/
      DATA GSPM   (82)/       8.3000000D0/
      DATA GPPM   (82)/       7.3000000D0/
      DATA GP2M   (82)/       6.5000000D0/
      DATA HSPM   (82)/       1.3000000D0/
      DATA DDM    (82)/       1.5526624D0/
      DATA QQM    (82)/       1.4488558D0/
      DATA AMM    (82)/       0.3601617D0/
      DATA ADM    (82)/       0.3239309D0/
      DATA AQM    (82)/       0.3502057D0/
C     DATA FOR ELEMENT 86        CONNECTION ATOM FOR QM/MM TREATMENTS
C     DATA REFMN  (86)/' **: (MNDO): I. ANTES, W. THIEL, J. PHYS. CHEM.
C    1A 103, 9290-9295 (1999)         '/
      DATA USSM   (86)/     -13.37151483D0/
      DATA UPPM   (86)/       0.00000000D0/
      DATA BETASM (86)/     -12.67678925D0/
      DATA BETAPM (86)/     -12.67678925D0/
      DATA ZSM    (86)/       0.93258788D0/
      DATA ZPM    (86)/       0.93258788D0/
      DATA ALPM   (86)/       1.64501728D0/
      DATA EISOLM (86)/     -13.37151483D0/
      DATA GSSM   (86)/      11.89194026D0/
      DATA GSPM   (86)/       0.00000000D0/
      DATA GPPM   (86)/       0.00000000D0/
      DATA GP2M   (86)/       0.00000000D0/
      DATA HSPM   (86)/       0.00000000D0/
      DATA DDM    (86)/       0.00000000D0/
      DATA QQM    (86)/       0.00000000D0/
      DATA AMM    (86)/       0.43704301D0/
      DATA ADM    (86)/       0.43704301D0/
      DATA AQM    (86)/       0.43704301D0/
C     ALTERNATIVE PARAMETER SET (UNPUBLISHED, NOT USED)
C     DATA USSM   (86)/     -12.75916989D0/
C     DATA BETASM (86)/     -11.99647103D0/
C     DATA BETAPM (86)/     -11.99647103D0/
C     DATA ZSM    (86)/       0.93559232D0/
C     DATA ZPM    (86)/       0.93559232D0/
C     DATA ALPM   (86)/       1.65897919D0/
C     DATA EISOLM (86)/     -12.75916989D0/
C     DATA GSSM   (86)/      10.85347036D0/
C     *
C *** START OF AM1 PARAMETERS.
C     *
C     NUMBER OF GAUSSIAN USED PER ATOM.
C     AM1 ACTUALLY USES ZERO GAUSSIANS FOR Zn, Ge, Sn, AND Hg.
C     FOR TECHNICAL REASONS IM1 IS SET TO 1 IN THESE CASES,
C     AND THE CORRESPONDING PARAMETERS ARE DEFINED TO BE ZERO.
      DATA IM1/3,0,0,0,3,4,3,2,2,0,2,3,1,3,3,3,2,0,
     1         0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,2,0,
     2         0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,2,0,
     3         25*0,                 1,0,0,0,0,0,4/
C     DATA FOR ELEMENT  1       AM1:   HYDROGEN
C     DATA REFAM  ( 1)/'  H: (AM1): M.J.S. DEWAR ET AL, J. AM. CHEM. SOC
C    1. 107 3902-3909 (1985)          '/
      DATA USSAM1( 1)/     -11.3964270D0/
      DATA USSAMD( 1)/     -11.2237910D0/
      DATA UPPAM1( 1)/       0.0000000D0/
      DATA UPPAMD( 1)/       0.0000000D0/
      DATA BETASA( 1)/      -6.1737870D0/
      DATA BSAMD ( 1)/      -6.3762650D0/
      DATA BETAPA( 1)/      -6.1737870D0/
      DATA BPAMD ( 1)/      -6.3762650D0/
      DATA ZSAM1 ( 1)/       1.1880780D0/
      DATA ZPAM1 ( 1)/       1.1880780D0/
      DATA ALPAM1( 1)/       2.8823240D0/
      DATA ALPAMD( 1)/       3.5777560D0/
      DATA EISOLA( 1)/     -11.3964270D0/
      DATA GSSAM1( 1)/      12.8480000D0/
      DATA GSPAM1( 1)/       0.0000000D0/
      DATA GPPAM1( 1)/       0.0000000D0/
      DATA GP2AM1( 1)/       0.0000000D0/
      DATA HSPAM1( 1)/       0.0000000D0/
      DATA DDAM1 ( 1)/       0.0000000D0/
      DATA QQAM1 ( 1)/       0.0000000D0/
      DATA AMAM1 ( 1)/       0.4721793D0/
      DATA ADAM1 ( 1)/       0.4721793D0/
      DATA AQAM1 ( 1)/       0.4721793D0/
      DATA GUESA1( 1,1)/       0.1227960D0/
      DATA GUESA2( 1,1)/       5.0000000D0/
      DATA GUESA3( 1,1)/       1.2000000D0/
      DATA GUESA1( 1,2)/       0.0050900D0/
      DATA GUESA2( 1,2)/       5.0000000D0/
      DATA GUESA3( 1,2)/       1.8000000D0/
      DATA GUESA1( 1,3)/      -0.0183360D0/
      DATA GUESA2( 1,3)/       2.0000000D0/
      DATA GUESA3( 1,3)/       2.1000000D0/
      DATA AM1C6 ( 1)/       0.1600000D0/
      DATA AM1R0 ( 1)/       0.1110000D0/
C     DATA FOR ELEMENT  5       AM1:   BORON  *
C     DATA REFAM  ( 5)/'  B: (AM1):  M.J.S. DEWAR, C. JIE, E. G. ZOEBISC
C    1H ORGANOMETALLICS 7, 513 (1988) '/
      DATA USSAM1(  5)/     -34.4928700D0/
      DATA UPPAM1(  5)/     -22.6315250D0/
      DATA BETASA(  5)/      -9.5991140D0/
      DATA BETAPA(  5)/      -6.2737570D0/
      DATA ZSAM1 (  5)/       1.6117090D0/
      DATA ZPAM1 (  5)/       1.5553850D0/
      DATA ALPAM1(  5)/       2.4469090D0/
      DATA EISOLA(  5)/     -63.7172650D0/
      DATA GSSAM1(  5)/      10.5900000D0/
      DATA GSPAM1(  5)/       9.5600000D0/
      DATA GPPAM1(  5)/       8.8600000D0/
      DATA GP2AM1(  5)/       7.8600000D0/
      DATA HSPAM1(  5)/       1.8100000D0/
      DATA DDAM1 (  5)/       0.9107622D0/
      DATA QQAM1 (  5)/       0.7874223D0/
      DATA AMAM1 (  5)/       0.3891951D0/
      DATA ADAM1 (  5)/       0.5045152D0/
      DATA AQAM1 (  5)/       0.5678856D0/
C     DATA FOR ELEMENT  6       AM1:   CARBON
C     DATA REFAM  ( 6)/'  C: (AM1): M.J.S. DEWAR ET AL, J. AM. CHEM. SOC
C    1. 107 3902-3909 (1985)          '/
      DATA USSAM1( 6)/     -52.0286580D0/
      DATA USSAMD( 6)/     -52.1837980D0/
      DATA UPPAM1( 6)/     -39.6142390D0/
      DATA UPPAMD( 6)/     -39.3684130D0/
      DATA BETASA( 6)/     -15.7157830D0/
      DATA BSAMD ( 6)/     -15.6823410D0/
      DATA BETAPA( 6)/      -7.7192830D0/
      DATA BPAMD ( 6)/      -7.8047620D0/
      DATA ZSAM1 ( 6)/       1.8086650D0/
      DATA ZPAM1 ( 6)/       1.6851160D0/
      DATA ALPAM1( 6)/       2.6482740D0/
      DATA ALPAMD( 6)/       2.6255060D0/
      DATA EISOLA( 6)/    -120.8157940D0/
      DATA GSSAM1( 6)/      12.2300000D0/
      DATA GSPAM1( 6)/      11.4700000D0/
      DATA GPPAM1( 6)/      11.0800000D0/
      DATA GP2AM1( 6)/       9.8400000D0/
      DATA HSPAM1( 6)/       2.4300000D0/
      DATA DDAM1 ( 6)/       0.8236736D0/
      DATA QQAM1 ( 6)/       0.7268015D0/
      DATA AMAM1 ( 6)/       0.4494671D0/
      DATA ADAM1 ( 6)/       0.6082946D0/
      DATA AQAM1 ( 6)/       0.6423492D0/
      DATA GUESA1( 6,1)/       0.0113550D0/
      DATA GUESA2( 6,1)/       5.0000000D0/
      DATA GUESA3( 6,1)/       1.6000000D0/
      DATA GUESA1( 6,2)/       0.0459240D0/
      DATA GUESA2( 6,2)/       5.0000000D0/
      DATA GUESA3( 6,2)/       1.8500000D0/
      DATA GUESA1( 6,3)/      -0.0200610D0/
      DATA GUESA2( 6,3)/       5.0000000D0/
      DATA GUESA3( 6,3)/       2.0500000D0/
      DATA GUESA1( 6,4)/      -0.0012600D0/
      DATA GUESA2( 6,4)/       5.0000000D0/
      DATA GUESA3( 6,4)/       2.6500000D0/
      DATA AM1C6 ( 6)/       1.6500000D0/
      DATA AM1R0 ( 6)/       0.1610000D0/
C     DATA FOR ELEMENT  7       AM1:   NITROGEN
C     DATA REFAM  ( 7)/'  N: (AM1): M.J.S. DEWAR ET AL, J. AM. CHEM. SOC
C    1. 107 3902-3909 (1985)          '/
      DATA USSAM1( 7)/     -71.8600000D0/
      DATA USSAMD( 7)/     -71.9978450D0/
      DATA UPPAM1( 7)/     -57.1675810D0/
      DATA UPPAMD( 7)/     -57.4017180D0/
      DATA BETASA( 7)/     -20.2991100D0/
      DATA BSAMD ( 7)/     -20.0924080D0/
      DATA BETAPA( 7)/     -18.2386660D0/
      DATA BPAMD ( 7)/     -18.4706790D0/
      DATA ZSAM1 ( 7)/       2.3154100D0/
      DATA ZPAM1 ( 7)/       2.1579400D0/
      DATA ALPAM1( 7)/       2.9472860D0/
      DATA ALPAMD( 7)/       2.9687370D0/
      DATA EISOLA( 7)/    -202.4077430D0/
      DATA GSSAM1( 7)/      13.5900000D0/
      DATA GSPAM1( 7)/      12.6600000D0/
      DATA GPPAM1( 7)/      12.9800000D0/
      DATA GP2AM1( 7)/      11.5900000D0/
      DATA HSPAM1( 7)/       3.1400000D0/
      DATA DDAM1 ( 7)/       0.6433247D0/
      DATA QQAM1 ( 7)/       0.5675528D0/
      DATA AMAM1 ( 7)/       0.4994487D0/
      DATA ADAM1 ( 7)/       0.7820840D0/
      DATA AQAM1 ( 7)/       0.7883498D0/
      DATA GUESA1( 7,1)/       0.0252510D0/
      DATA GUESA2( 7,1)/       5.0000000D0/
      DATA GUESA3( 7,1)/       1.5000000D0/
      DATA GUESA1( 7,2)/       0.0289530D0/
      DATA GUESA2( 7,2)/       5.0000000D0/
      DATA GUESA3( 7,2)/       2.1000000D0/
      DATA GUESA1( 7,3)/      -0.0058060D0/
      DATA GUESA2( 7,3)/       2.0000000D0/
      DATA GUESA3( 7,3)/       2.4000000D0/
      DATA AM1C6 ( 7)/       1.1100000D0/
      DATA AM1R0 ( 7)/       0.1550000D0/
C     DATA FOR ELEMENT  8       AM1:   OXYGEN
C     DATA REFAM  ( 8)/'  O: (AM1): M.J.S. DEWAR ET AL, J. AM. CHEM. SOC
C    1. 107 3902-3909 (1985)          '/
      DATA USSAM1( 8)/     -97.8300000D0/
      DATA USSAMD( 8)/     -97.6105880D0/
      DATA UPPAM1( 8)/     -78.2623800D0/
      DATA UPPAMD( 8)/     -78.5897000D0/
      DATA BETASA( 8)/     -29.2727730D0/
      DATA BSAMD ( 8)/     -29.5024810D0/
      DATA BETAPA( 8)/     -29.2727730D0/
      DATA BPAMD ( 8)/     -29.4953800D0/
      DATA ZSAM1 ( 8)/       3.1080320D0/
      DATA ZPAM1 ( 8)/       2.5240390D0/
      DATA ALPAM1( 8)/       4.4553710D0/
      DATA ALPAMD( 8)/       4.6336990D0/
      DATA EISOLA( 8)/    -316.0995200D0/
      DATA GSSAM1( 8)/      15.4200000D0/
      DATA GSPAM1( 8)/      14.4800000D0/
      DATA GPPAM1( 8)/      14.5200000D0/
      DATA GP2AM1( 8)/      12.9800000D0/
      DATA HSPAM1( 8)/       3.9400000D0/
      DATA DDAM1 ( 8)/       0.4988896D0/
      DATA QQAM1 ( 8)/       0.4852322D0/
      DATA AMAM1 ( 8)/       0.5667034D0/
      DATA ADAM1 ( 8)/       0.9961066D0/
      DATA AQAM1 ( 8)/       0.9065223D0/
      DATA GUESA1( 8,1)/       0.2809620D0/
      DATA GUESA2( 8,1)/       5.0000000D0/
      DATA GUESA3( 8,1)/       0.8479180D0/
      DATA GUESA1( 8,2)/       0.0814300D0/
      DATA GUESA2( 8,2)/       7.0000000D0/
      DATA GUESA3( 8,2)/       1.4450710D0/
      DATA AM1C6 ( 8)/       0.7000000D0/
      DATA AM1R0 ( 8)/       0.1490000D0/
C     DATA FOR ELEMENT  9       AM1:   FLUORINE  *
C     DATA REFAM  ( 9)/'  F: (AM1): M.J.S.DEWAR AND E.G.ZOEBISCH, J.MOL.
C    1STRUCT. 180, 1 (1988)           '/
      DATA USSAM1( 9)/    -136.1055790D0/
      DATA UPPAM1( 9)/    -104.8898850D0/
      DATA BETASA( 9)/     -69.5902770D0/
      DATA BETAPA( 9)/     -27.9223600D0/
      DATA ZSAM1 ( 9)/       3.7700820D0/
      DATA ZPAM1 ( 9)/       2.4946700D0/
      DATA ALPAM1( 9)/       5.5178000D0/
      DATA EISOLA( 9)/    -482.2905830D0/
      DATA GSSAM1( 9)/      16.9200000D0/
      DATA GSPAM1( 9)/      17.2500000D0/
      DATA GPPAM1( 9)/      16.7100000D0/
      DATA GP2AM1( 9)/      14.9100000D0/
      DATA HSPAM1( 9)/       4.8300000D0/
      DATA DDAM1 ( 9)/       0.4145203D0/
      DATA QQAM1 ( 9)/       0.4909446D0/
      DATA AMAM1 ( 9)/       0.6218302D0/
      DATA ADAM1 ( 9)/       1.2088792D0/
      DATA AQAM1 ( 9)/       0.9449355D0/
      DATA GUESA1( 9,1)/       0.2420790D0/
      DATA GUESA2( 9,1)/       4.8000000D0/
      DATA GUESA3( 9,1)/       0.9300000D0/
      DATA GUESA1( 9,2)/       0.0036070D0/
      DATA GUESA2( 9,2)/       4.6000000D0/
      DATA GUESA3( 9,2)/       1.6600000D0/
C     DATA FOR ELEMENT 11        SODIUM
C     DATA REFMN  (11)/' Na: (AM1):   E.N. BROTHERS, K.M. MERZ, J.PHYS.C
C    1HEM. B 106, 2779 (2002).        '/
      DATA USSAM1 (11)/      -5.25553620D0/
      DATA UPPAM1 (11)/      -2.08127810D0/
      DATA BETASA (11)/      -1.45369440D0/
      DATA BETAPA (11)/      -0.22980640D0/
      DATA ZSAM1  (11)/       0.67977920D0/
      DATA ZPAM1  (11)/       1.21704680D0/
C     DATA ZDAM1  (11)/       1.00000000D0/
      DATA ALPAM1 (11)/       2.24871640D0/
      DATA EISOLA (11)/      -5.25553620D0/
      DATA GSSAM1 (11)/       7.34591780D0/
      DATA GSPAM1 (11)/       8.40425500D0/
      DATA GPPAM1 (11)/       4.11305160D0/
      DATA GP2AM1 (11)/       4.03709570D0/
      DATA HSPAM1 (11)/       0.24876990D0/
      DATA DDAM1  (11)/       1.58997550D0/
      DATA QQAM1  (11)/       1.37490190D0/
      DATA AMAM1  (11)/       0.26997130D0/
      DATA ADAM1  (11)/       0.16180470D0/
      DATA AQAM1  (11)/       0.18657380D0/
      DATA GUESA1(11,1)/      0.53226680D0/
      DATA GUESA2(11,1)/      0.48003040D0/
      DATA GUESA3(11,1)/      1.16810550D0/
      DATA GUESA1(11,2)/      0.92235980D0/
      DATA GUESA2(11,2)/      1.90767760D0/
      DATA GUESA3(11,2)/      1.15376700D0/
C     DATA FOR ELEMENT 12        MAGNESIUM
C     DATA REFMN  (12)/' Mg: (AM1):   M.C. HUTTER, J.R. REIMERS, N.S. HU
C    1SH, J.PHYS.CHEM. 102, 8080(1998)'/
      DATA USSAM1 (12)/     -14.96959313D0/
      DATA UPPAM1 (12)/     -11.56229248D0/
      DATA BETASA (12)/      -1.25974355D0/
      DATA BETAPA (12)/      -0.77836604D0/
      DATA ZSAM1  (12)/       1.22339270D0/
      DATA ZPAM1  (12)/       1.02030798D0/
C     DATA ZDAM1  (12)/       1.00000000D0/
      DATA ALPAM1 (12)/       1.67049799D0/
      DATA EISOLA (12)/     -22.43786349D0/
      DATA GSSAM1 (12)/       7.50132277D0/
CJPC  DATA GSPAM1 (12)/       6.22766502D0/
CJPC  DATA GPPAM1 (12)/       6.34591536D0/
      DATA GSPAM1 (12)/       6.34591536D0/
      DATA GPPAM1 (12)/       4.77534467D0/
      DATA GP2AM1 (12)/       4.34017279D0/
      DATA HSPAM1 (12)/       0.48930466D0/
      DATA DDAM1  (12)/       1.75012115D0/
      DATA QQAM1  (12)/       1.64001467D0/
      DATA AMAM1  (12)/       0.27568257D0/
      DATA ADAM1  (12)/       0.19957236D0/
CJPC  DATA AQAM1  (12)/       0.46064361D0/
      DATA AQAM1  (12)/       0.26440152D0/
      DATA GUESA1(12,1)/      2.55017735D0/
      DATA GUESA2(12,1)/      4.29397225D0/
      DATA GUESA3(12,1)/      0.79989601D0/
      DATA GUESA1(12,2)/     -0.00565806D0/
      DATA GUESA2(12,2)/      2.96053910D0/
      DATA GUESA3(12,2)/      1.47499983D0/
      DATA GUESA1(12,3)/     -0.00610286D0/
      DATA GUESA2(12,3)/      2.61416919D0/
      DATA GUESA3(12,3)/      2.42604040D0/
C     DATA FOR ELEMENT 13          :   ALUMINUM  *
C     DATA REFAM  (13)/' Al: (AM1):  M. J. S. Dewar, A. J. Holder, Organ
C    1ometallics, 9, 508-511 (1990).  '/
      DATA USSAM1(13)/     -24.3535850D0/
      DATA UPPAM1(13)/     -18.3636450D0/
      DATA BETASA(13)/      -3.8668220D0/
      DATA BETAPA(13)/      -2.3171460D0/
      DATA ZSAM1 (13)/       1.5165930D0/
      DATA ZPAM1 (13)/       1.3063470D0/
C     DATA ZD    (13)/       1.0000000D0/
      DATA ALPAM1(13)/       1.9765860D0/
      DATA EISOLA(13)/     -46.4208150D0/
      DATA GSSAM1(13)/       8.0900000D0/
      DATA GSPAM1(13)/       6.6300000D0/
      DATA GPPAM1(13)/       5.9800000D0/
      DATA GP2AM1(13)/       5.4000000D0/
      DATA HSPAM1(13)/       0.7000000D0/
      DATA DDAM1 (13)/       1.4040443D0/
      DATA QQAM1 (13)/       1.2809154D0/
      DATA AMAM1 (13)/       0.2973172D0/
      DATA ADAM1 (13)/       0.2630229D0/
      DATA AQAM1 (13)/       0.3427832D0/
      DATA GUESA1(13,1)/       0.0900000D0/
      DATA GUESA2(13,1)/      12.3924430D0/
      DATA GUESA3(13,1)/       2.0503940D0/
C     DATA FOR ELEMENT 14       AM1:   SILICON  *
C     DATA REFAM  (14)/' Si: (AM1): M.J.S.DEWAR, C. JIE, ORGANOMETALLICS
C    1, 6, 1486-1490 (1987).          '/
      DATA USSAM1(14)/     -33.9536220D0/
      DATA UPPAM1(14)/     -28.9347490D0/
      DATA BETASA(14)/      -3.784852D0/
      DATA BETAPA(14)/      -1.968123D0/
      DATA ZSAM1 (14)/       1.830697D0/
      DATA ZPAM1 (14)/       1.2849530D0/
C     DATA ZD    (14)/       1.0000000D0/
      DATA ALPAM1(14)/       2.257816D0/
      DATA EISOLA(14)/     -79.0017420D0/
      DATA GSSAM1(14)/       9.8200000D0/
      DATA GSPAM1(14)/       8.3600000D0/
      DATA GPPAM1(14)/       7.3100000D0/
      DATA GP2AM1(14)/       6.5400000D0/
      DATA HSPAM1(14)/       1.3200000D0/
      DATA DDAM1 (14)/       1.1631107D0/
      DATA QQAM1 (14)/       1.3022422D0/
      DATA AMAM1 (14)/       0.3608967D0/
      DATA ADAM1 (14)/       0.3829813D0/
      DATA AQAM1 (14)/       0.3712106D0/
      DATA GUESA1(14,1)/       0.25D0/
      DATA GUESA2(14,1)/       9.000D0/
      DATA GUESA3(14,1)/       0.911453D0/
      DATA GUESA1(14,2)/       0.061513D0/
      DATA GUESA2(14,2)/       5.00D0/
      DATA GUESA3(14,2)/       1.995569D0/
      DATA GUESA1(14,3)/       0.0207890D0/
      DATA GUESA2(14,3)/       5.00D0/
      DATA GUESA3(14,3)/       2.990610D0/
C     DATA FOR ELEMENT 15      PHOSPHORUS
C     DATA REFAM  (15)/'  P: (AM1): M.J.S.DEWAR AND C.JIE, J.MOL.STRUCT.
C    1 187, 1 (1989)                  '/
      DATA USSAM1(15)/     -42.0298630D0/
      DATA UPPAM1(15)/     -34.0307090D0/
      DATA BETASA(15)/     - 6.3537640D0/
      DATA BETAPA(15)/      -6.5907090D0/
      DATA ZSAM1 (15)/       1.9812800D0/
      DATA ZPAM1 (15)/       1.8751500D0/
C     DATA ZD    (15)/       1.0000000D0/
      DATA ALPAM1(15)/       2.4553220D0/
      DATA EISOLA(15)/    -124.4368355D0/
      DATA GSSAM1(15)/      11.5600050D0/
      DATA GSPAM1(15)/       5.2374490D0/
      DATA GPPAM1(15)/       7.8775890D0/
      DATA GP2AM1(15)/       7.3076480D0/
      DATA HSPAM1(15)/       0.7792380D0/
      DATA DDAM1 (15)/       1.0452022D0/
      DATA QQAM1 (15)/       0.8923660D0/
      DATA AMAM1 (15)/       0.4248440D0/
      DATA ADAM1 (15)/       0.3275319D0/
      DATA AQAM1 (15)/       0.4386854D0/
      DATA GUESA1(15,1)/      -0.0318270D0/
      DATA GUESA2(15,1)/       6.0000000D0/
      DATA GUESA3(15,1)/       1.4743230D0/
      DATA GUESA1(15,2)/       0.0184700D0/
      DATA GUESA2(15,2)/       7.0000000D0/
      DATA GUESA3(15,2)/       1.7793540D0/
      DATA GUESA1(15,3)/       0.0332900D0/
      DATA GUESA2(15,3)/       9.0000000D0/
      DATA GUESA3(15,3)/       3.0065760D0/
C     DATA FOR ELEMENT 16          :   SULFUR  *
C     DATA REFAM  (16)/'  S: (AM1): M.J.S.DEWAR, Y-C YUAN, THEOCHEM, IN
C    1 PRESS                          '/
      DATA USSAM1(16)/     -56.6940560D0/
      DATA UPPAM1(16)/     -48.7170490D0/
      DATA BETASA(16)/      -3.9205660D0/
      DATA BETAPA(16)/      -7.9052780D0/
      DATA ZSAM1 (16)/       2.3665150D0/
      DATA ZPAM1 (16)/       1.6672630D0/
C     DATA ZD    (16)/       1.0000000D0/
      DATA ALPAM1(16)/       2.4616480D0/
      DATA EISOLA(16)/    -191.7321930D0/
      DATA GSSAM1(16)/      11.7863290D0/
      DATA GSPAM1(16)/       8.6631270D0/
      DATA GPPAM1(16)/      10.0393080D0/
      DATA GP2AM1(16)/       7.7816880D0/
      DATA HSPAM1(16)/       2.5321370D0/
      DATA DDAM1 (16)/       0.9004265D0/
      DATA QQAM1 (16)/       1.0036329D0/
      DATA AMAM1 (16)/       0.4331617D0/
      DATA ADAM1 (16)/       0.5907115D0/
      DATA AQAM1 (16)/       0.6454943D0/
      DATA GUESA1(16,1)/      -0.5091950D0/
      DATA GUESA2(16,1)/       4.5936910D0/
      DATA GUESA3(16,1)/       0.7706650D0/
      DATA GUESA1(16,2)/      -0.0118630D0/
      DATA GUESA2(16,2)/       5.8657310D0/
      DATA GUESA3(16,2)/       1.5033130D0/
      DATA GUESA1(16,3)/       0.0123340D0/
      DATA GUESA2(16,3)/      13.5573360D0/
      DATA GUESA3(16,3)/       2.0091730D0/
C     DATA FOR ELEMENT 17       AM1:   CHLORINE  *
C     DATA REFAM  (17)/' Cl: (AM1): M.J.S.DEWAR AND E.G.ZOEBISCH, J.MOL.
C    1STRUCT. 180, 1 (1988)           '/
      DATA USSAM1(17)/    -111.6139480D0/
      DATA UPPAM1(17)/     -76.6401070D0/
      DATA BETASA(17)/     -24.5946700D0/
      DATA BETAPA(17)/     -14.6372160D0/
      DATA ZSAM1 (17)/       3.6313760D0/
      DATA ZPAM1 (17)/       2.0767990D0/
C     DATA ZD    (17)/       1.0000000D0/
      DATA ALPAM1(17)/       2.9193680D0/
      DATA EISOLA(17)/    -372.1984310D0/
      DATA GSSAM1(17)/      15.0300000D0/
      DATA GSPAM1(17)/      13.1600000D0/
      DATA GPPAM1(17)/      11.3000000D0/
      DATA GP2AM1(17)/       9.9700000D0/
      DATA HSPAM1(17)/       2.4200000D0/
      DATA DDAM1 (17)/       0.5406286D0/
      DATA QQAM1 (17)/       0.8057208D0/
      DATA AMAM1 (17)/       0.5523705D0/
      DATA ADAM1 (17)/       0.7693200D0/
      DATA AQAM1 (17)/       0.6133369D0/
      DATA GUESA1(17,1)/       0.0942430D0/
      DATA GUESA2(17,1)/       4.0000000D0/
      DATA GUESA3(17,1)/       1.3000000D0/
      DATA GUESA1(17,2)/       0.0271680D0/
      DATA GUESA2(17,2)/       4.0000000D0/
      DATA GUESA3(17,2)/       2.1000000D0/
C     DATA FOR ELEMENT 30        ZINC
C     DATA REFAM  (30)/' Zn: (   ):  M.J.S. DEWAR, K.M. MERZ, ORGANOMET
C    1ALLICS, 7, 522-524 (1988)       '/
      DATA USSAM1(30)/     -21.0400080D0/
      DATA UPPAM1(30)/     -17.6555740D0/
      DATA BETASA(30)/      -1.9974290D0/
      DATA BETAPA(30)/      -4.7581190D0/
      DATA ZSAM1 (30)/       1.9542990D0/
      DATA ZPAM1 (30)/       1.3723650D0/
C     DATA ZD    (30)/       1.0000000D0/
      DATA ALPAM1(30)/       1.4845630D0/
      DATA EISOLA(30)/     -30.2800160D0/
      DATA GSSAM1(30)/      11.8000000D0/
      DATA GSPAM1(30)/      11.1820180D0/
      DATA GPPAM1(30)/      13.3000000D0/
      DATA GP2AM1(30)/      12.9305200D0/
      DATA HSPAM1(30)/       0.4846060D0/
      DATA DDAM1 (30)/       1.3581113D0/
      DATA QQAM1 (30)/       1.5457406D0/
      DATA AMAM1 (30)/       0.4336641D0/
      DATA ADAM1 (30)/       0.2317423D0/
      DATA AQAM1 (30)/       0.2621165D0/
      DATA GUESA1(30,1)/     0.0000000D0/
      DATA GUESA2(30,1)/     0.0000000D0/
      DATA GUESA3(30,1)/     0.0000000D0/
C     DATA FOR ELEMENT 32        GERMANIUM
C     DATA REFAM  (32)/' Ge: (   ): M.J.S.Dewar and C.Jie, Organometalli
C    1cs, 8, 1544, (1989)             '/
      DATA USSAM1(32)/     -34.1838890D0/
      DATA UPPAM1(32)/     -28.6408110D0/
      DATA BETASA(32)/      -4.3566070D0/
      DATA BETAPA(32)/      -0.9910910D0/
      DATA ZSAM1 (32)/       1.2196310D0/
      DATA ZPAM1 (32)/       1.9827940D0/
      DATA ALPAM1(32)/       2.1364050D0/
      DATA EISOLA(32)/     -78.7084810D0/
      DATA GSSAM1(32)/      10.1686050D0/
      DATA GSPAM1(32)/       8.1444730D0/
      DATA GPPAM1(32)/       6.6719020D0/
      DATA GP2AM1(32)/       6.2697060D0/
      DATA HSPAM1(32)/       0.9370930D0/
      DATA DDAM1 (32)/       1.2472095D0/
      DATA QQAM1 (32)/       1.0698642D0/
      DATA AMAM1 (32)/       0.3737084D0/
      DATA ADAM1 (32)/       0.3180309D0/
      DATA AQAM1 (32)/       0.3485612D0/
      DATA GUESA1(32,1)/     0.0000000D0/
      DATA GUESA2(32,1)/     0.0000000D0/
      DATA GUESA3(32,1)/     0.0000000D0/
C     DATA FOR ELEMENT 35       AM1:   BROMINE  *
C     DATA REFAM  (35)/' Br: (AM1): M.J.S.DEWAR AND E.G.ZOEBISCH, J.MOL.
C    1STRUCT. 180, 1 (1988)           '/
      DATA USSAM1(35)/    -104.6560630D0/
      DATA UPPAM1(35)/     -74.9300520D0/
      DATA BETASA(35)/     -19.3998800D0/
      DATA BETAPA(35)/      -8.9571950D0/
      DATA ZSAM1 (35)/       3.0641330D0/
      DATA ZPAM1 (35)/       2.0383330D0/
C     DATA ZD    (35)/       1.0000000D0/
      DATA ALPAM1(35)/       2.5765460D0/
      DATA EISOLA(35)/    -352.3142087D0/
      DATA GSSAM1(35)/      15.0364395D0/
      DATA GSPAM1(35)/      13.0346824D0/
      DATA GPPAM1(35)/      11.2763254D0/
      DATA GP2AM1(35)/       9.8544255D0/
      DATA HSPAM1(35)/       2.4558683D0/
      DATA DDAM1 (35)/       0.8458104D0/
      DATA QQAM1 (35)/       1.0407133D0/
      DATA AMAM1 (35)/       0.5526071D0/
      DATA ADAM1 (35)/       0.6024598D0/
      DATA AQAM1 (35)/       0.5307555D0/
      DATA GUESA1(35,1)/       0.0666850D0/
      DATA GUESA2(35,1)/       4.0000000D0/
      DATA GUESA3(35,1)/       1.5000000D0/
      DATA GUESA1(35,2)/       0.0255680D0/
      DATA GUESA2(35,2)/       4.0000000D0/
      DATA GUESA3(35,2)/       2.3000000D0/
C     DATA FOR ELEMENT 50        TIN
C     DATA REFAD  (50)/' Sn: (AM1): M.J.S. DEWAR ET AL, ORGANOMETALLICS,
C    1 10, 431, (1991)                '/
      DATA USSAM1(50)/     -35.4967410D0/
      DATA UPPAM1(50)/     -28.0976360D0/
      DATA ZSAM1 (50)/       2.5993760D0/
      DATA ZPAM1 (50)/       1.6959620D0/
      DATA BETASA(50)/      -3.2350000D0/
      DATA BETAPA(50)/      -2.5778900D0/
      DATA ALPAM1(50)/       1.8369360D0/
      DATA EISOLA(50)/     -80.6887540D0/
      DATA GSSAM1(50)/       9.8000000D0/
      DATA GSPAM1(50)/       8.3000000D0/
      DATA GPPAM1(50)/       7.3000000D0/
      DATA GP2AM1(50)/       6.5000000D0/
      DATA HSPAM1(50)/       1.3000000D0/
C     DATA HPPAM1(50)/       0.4000000D0/
      DATA DDAM1 (50)/       1.1528230D0/
      DATA QQAM1 (50)/       1.5148020D0/
      DATA AMAM1 (50)/       0.3601620D0/
      DATA ADAM1 (50)/       0.3823780D0/
      DATA AQAM1 (50)/       0.3400760D0/
      DATA GUESA1(50,1)/     0.0000000D0/
      DATA GUESA2(50,1)/     0.0000000D0/
      DATA GUESA3(50,1)/     0.0000000D0/
C     DATA FOR ELEMENT 53       AM1:   IODINE  *
C     DATA REFAM  (53)/'  I: (AM1): M.J.S.DEWAR AND E.G.ZOEBISCH, J.MOL.
C    1STRUCT. 180, 1 (1988)           '/
      DATA USSAM1(53)/    -103.5896630D0/
      DATA UPPAM1(53)/     -74.4299970D0/
      DATA BETASA(53)/      -8.4433270D0/
      DATA BETAPA(53)/      -6.3234050D0/
      DATA ZSAM1 (53)/       2.1028580D0/
      DATA ZPAM1 (53)/       2.1611530D0/
C     DATA ZD    (53)/       1.0000000D0/
      DATA ALPAM1(53)/       2.2994240D0/
      DATA EISOLA(53)/    -346.8642857D0/
      DATA GSSAM1(53)/      15.0404486D0/
      DATA GSPAM1(53)/      13.0565580D0/
      DATA GPPAM1(53)/      11.1477837D0/
      DATA GP2AM1(53)/       9.9140907D0/
      DATA HSPAM1(53)/       2.4563820D0/
      DATA DDAM1 (53)/       1.4878778D0/
      DATA QQAM1 (53)/       1.1887388D0/
      DATA AMAM1 (53)/       0.5527544D0/
      DATA ADAM1 (53)/       0.4497523D0/
      DATA AQAM1 (53)/       0.4631775D0/
      DATA GUESA1(53,1)/       0.0043610D0/
      DATA GUESA2(53,1)/       2.3000000D0/
      DATA GUESA3(53,1)/       1.8000000D0/
      DATA GUESA1(53,2)/       0.0157060D0/
      DATA GUESA2(53,2)/       3.0000000D0/
      DATA GUESA3(53,2)/       2.2400000D0/
C     DATA FOR ELEMENT 80        MERCURY
C     DATA REFAM  (80)/' Hg: (   ): M.J.S.Dewar and C.Jie, Organometalli
C    1cs 8, 1547, (1989)              '/
      DATA USSAM1(80)/     -19.9415780D0/
      DATA UPPAM1(80)/     -11.1108700D0/
      DATA BETASA(80)/      -0.9086570D0/
      DATA BETAPA(80)/      -4.9093840D0/
      DATA ZSAM1 (80)/       2.0364130D0/
      DATA ZPAM1 (80)/       1.9557660D0/
      DATA ALPAM1(80)/       1.4847340D0/
      DATA EISOLA(80)/     -29.0831560D0/
      DATA GSSAM1(80)/      10.8000000D0/
      DATA GSPAM1(80)/       9.3000000D0/
      DATA GPPAM1(80)/      14.3000000D0/
      DATA GP2AM1(80)/      13.5000000D0/
      DATA HSPAM1(80)/       1.3000000D0/
      DATA DDAM1 (80)/       1.8750829D0/
      DATA QQAM1 (80)/       1.5424241D0/
      DATA AMAM1 (80)/       0.3969129D0/
      DATA ADAM1 (80)/       0.2926605D0/
      DATA AQAM1 (80)/       0.3360599D0/
      DATA GUESA1(80,1)/     0.0000000D0/
      DATA GUESA2(80,1)/     0.0000000D0/
      DATA GUESA3(80,1)/     0.0000000D0/
C     DATA FOR ELEMENT 86        CONNECTION ATOM FOR QM/MM TREATMENTS
C     DATA REFAM  (86)/' **: (AM1): I. ANTES, W. THIEL, J. PHYS. CHEM. A
C    1 103, 9290-9295 (1999)          '/
      DATA USSAM1(86)/     -11.83958268D0/
      DATA UPPAM1(86)/       0.00000000D0/
      DATA BETASA(86)/     -11.20857778D0/
      DATA BETAPA(86)/     -11.20857778D0/
      DATA ZSAM1 (86)/       1.15644960D0/
      DATA ZPAM1 (86)/       1.15644960D0/
      DATA ALPAM1(86)/       1.46212100D0/
      DATA EISOLA(86)/     -11.83958268D0/
      DATA GSSAM1(86)/      11.39327472D0/
      DATA GSPAM1(86)/       0.00000000D0/
      DATA GPPAM1(86)/       0.00000000D0/
      DATA GP2AM1(86)/       0.00000000D0/
      DATA HSPAM1(86)/       0.00000000D0/
      DATA DDAM1 (86)/       0.00000000D0/
      DATA QQAM1 (86)/       0.00000000D0/
      DATA AMAM1 (86)/       0.41871645D0/
      DATA ADAM1 (86)/       0.41871645D0/
      DATA AQAM1 (86)/       0.41871645D0/
      DATA GUESA1(86,1)/    -0.30401893D0/
      DATA GUESA2(86,1)/     5.07141717D0/
      DATA GUESA3(86,1)/     1.62230340D0/
      DATA GUESA1(86,2)/     0.05480783D0/
      DATA GUESA2(86,2)/     5.21442318D0/
      DATA GUESA3(86,2)/     1.97508768D0/
      DATA GUESA1(86,3)/    -0.25449049D0/
      DATA GUESA2(86,3)/     5.10879333D0/
      DATA GUESA3(86,3)/     2.26466753D0/
      DATA GUESA1(86,4)/    -0.11453638D0/
      DATA GUESA2(86,4)/     4.70889286D0/
      DATA GUESA3(86,4)/     2.92552463D0/
C     *
C     START OF MNDO-PM3 PARAMETER SET.
C     *
C     NUMBER OF GAUSSIANS USED PER ATOM.
C     PM3 ACTUALLY USES ZERO GAUSSIANS FOR CD.
C     FOR TECHNICAL REASONS IM1 IS SET TO 1 IN THIS CASE,
C     AND THE CORRESPONDING PARAMETERS ARE DEFINED TO BE ZERO.
      DATA IM3/2,0,2,2,0,2,2,2,2,0,2,2,2,2,2,2,2,0,
     1         0,0,0,0,0,0,0,0,0,0,0,2,2,2,2,2,2,0,
     2         0,0,0,0,0,0,0,0,0,0,0,1,2,2,2,2,2,0,
     3         25*0,                 2,2,2,2,0,0,2/
C     DATA FOR ELEMENT  1        HYDROGEN
C     DATA REFPM3 ( 1)/'  H: (PM3): J. J. P. STEWART, J. COMP. CHEM.
C    1 10, 209 (1989)                 '/
      DATA USSPM3(  1)/     -13.0733210D0/
      DATA USSPMD(  1)/     -13.0540760D0/
      DATA UPPPM3(  1)/       0.0000000D0/
      DATA UPPPMD(  1)/       0.0000000D0/
      DATA BETASP(  1)/      -5.6265120D0/
      DATA BSPMD (  1)/      -5.6289010D0/
      DATA BETAPP(  1)/      -5.6265120D0/
      DATA BPPMD (  1)/      -5.6289010D0/
      DATA ZSPM3 (  1)/       0.9678070D0/
      DATA ZPPM3 (  1)/       0.9678070D0/
      DATA ALPPM3(  1)/       3.3563860D0/
      DATA ALPPMD(  1)/       3.4175320D0/
      DATA EISOLP(  1)/     -13.0733210D0/
      DATA GSSPM3(  1)/      14.7942080D0/
      DATA GSPPM3(  1)/       0.0000000D0/
      DATA GPPPM3(  1)/       0.0000000D0/
      DATA GP2PM3(  1)/       0.0000000D0/
      DATA HSPPM3(  1)/       0.0000000D0/
      DATA DDPM3 (  1)/       0.0000000D0/
      DATA QQPM3 (  1)/       0.0000000D0/
      DATA AMPM3 (  1)/       0.5437048D0/
      DATA ADPM3 (  1)/       0.5437048D0/
      DATA AQPM3 (  1)/       0.5437048D0/
      DATA GUESP1(  1,1)/       1.1287500D0/
      DATA GUESP2(  1,1)/       5.0962820D0/
      DATA GUESP3(  1,1)/       1.5374650D0/
      DATA GUESP1(  1,2)/      -1.0603290D0/
      DATA GUESP2(  1,2)/       6.0037880D0/
      DATA GUESP3(  1,2)/       1.5701890D0/
      DATA PM3C6 ( 1)/       0.1600000D0/
      DATA PM3R0 ( 1)/       0.1110000D0/
C     DATA FOR ELEMENT  3        LITHIUM
C     DATA REFMN  ( 3)/' Li: (PM3): E. ANDERS, R. KOCH, P. FREUNSCHT, J.
C    1 COMP. CHEM. 14, 1301 (1993)    '/
      DATA USSPM3 ( 3)/      -5.3000000D0/
      DATA UPPPM3 ( 3)/      -3.4000000D0/
      DATA BETASP ( 3)/      -0.5500000D0/
      DATA BETAPP ( 3)/      -1.5000000D0/
      DATA ZSPM3  ( 3)/       0.6500000D0/
      DATA ZPPM3  ( 3)/       0.7500000D0/
      DATA ALPPM3 ( 3)/       1.2550000D0/
      DATA EISOLP ( 3)/      -5.3000000D0/
      DATA GSSPM3 ( 3)/       4.5000000D0/
      DATA GSPPM3 ( 3)/       3.0000000D0/
      DATA GPPPM3 ( 3)/       5.2500000D0/
      DATA GP2PM3 ( 3)/       4.5000000D0/
      DATA HSPPM3 ( 3)/       0.1500000D0/
      DATA DDPM3  ( 3)/       2.0357652D0/
      DATA QQPM3  ( 3)/       1.6329932D0/
      DATA AMPM3  ( 3)/       0.1653804D0/
      DATA ADPM3  ( 3)/       0.1156738D0/
      DATA AQPM3  ( 3)/       0.3166080D0/
      DATA GUESP1(  3,1)/      -0.4500000D0/
      DATA GUESP2(  3,1)/       5.0000000D0/
      DATA GUESP3(  3,1)/       1.0000000D0/
      DATA GUESP1(  3,2)/       0.8000000D0/
      DATA GUESP2(  3,2)/       6.5000000D0/
      DATA GUESP3(  3,2)/       1.0000000D0/
C     DATA FOR ELEMENT  4        BERYLLIUM
C     DATA REFPM3( 4)/ ' Be: (PM3): J. J. P. STEWART, J. COMP. CHEM.
C    1 12, 320 (1991)                 '/
      DATA USSPM3(  4)/     -17.2647520D0/
      DATA UPPPM3(  4)/     -11.3042430D0/
      DATA BETASP(  4)/      -3.9620530D0/
      DATA BETAPP(  4)/      -2.7806840D0/
      DATA ZSPM3 (  4)/       0.8774390D0/
      DATA ZPPM3 (  4)/       1.5087550D0/
      DATA ALPPM3(  4)/       1.5935360D0/
      DATA EISOLP(  4)/     -25.5166530D0/
      DATA GSSPM3(  4)/       9.0128510D0/
      DATA GSPPM3(  4)/       6.5761990D0/
      DATA GPPPM3(  4)/       6.0571820D0/
      DATA GP2PM3(  4)/       9.0052190D0/
      DATA HSPPM3(  4)/       0.5446790D0/
      DATA DDPM3 (  4)/       1.0090531D0/
      DATA QQPM3 (  4)/       0.8117586D0/
      DATA AMPM3 (  4)/       0.3312330D0/
      DATA ADPM3 (  4)/       0.2908996D0/
      DATA AQPM3 (  4)/       0.3530008D0/
      DATA GUESP1(  4,1)/       1.6315720D0/
      DATA GUESP2(  4,1)/       2.6729620D0/
      DATA GUESP3(  4,1)/       1.7916860D0/
      DATA GUESP1(  4,2)/      -2.1109590D0/
      DATA GUESP2(  4,2)/       1.9685940D0/
      DATA GUESP3(  4,2)/       1.7558710D0/
C     DATA FOR ELEMENT  6      CARBON
C     DATA REFPM3 ( 6)/'  C: (PM3): J. J. P. STEWART, J. COMP. CHEM.
C    1 10, 209 (1989)                 '/
      DATA USSPM3(  6)/     -47.2703200D0/
      DATA USSPMD(  6)/     -47.2754310D0/
      DATA UPPPM3(  6)/     -36.2669180D0/
      DATA UPPPMD(  6)/     -36.2689160D0/
      DATA BETASP(  6)/     -11.9100150D0/
      DATA BSPMD (  6)/     -11.9414660D0/
      DATA BETAPP(  6)/      -9.8027550D0/
      DATA BPPMD (  6)/      -9.8197600D0/
      DATA ZSPM3 (  6)/       1.5650850D0/
      DATA ZPPM3 (  6)/       1.8423450D0/
      DATA ALPPM3(  6)/       2.7078070D0/
      DATA ALPPMD(  6)/       2.7211520D0/
      DATA EISOLP(  6)/    -111.2299170D0/
      DATA GSSPM3(  6)/      11.2007080D0/
      DATA GSPPM3(  6)/      10.2650270D0/
      DATA GPPPM3(  6)/      10.7962920D0/
      DATA GP2PM3(  6)/       9.0425660D0/
      DATA HSPPM3(  6)/       2.2909800D0/
      DATA DDPM3 (  6)/       0.8332396D0/
      DATA QQPM3 (  6)/       0.6647750D0/
      DATA AMPM3 (  6)/       0.4116394D0/
      DATA ADPM3 (  6)/       0.5885862D0/
      DATA AQPM3 (  6)/       0.7647667D0/
      DATA GUESP1(  6,1)/       0.0501070D0/
      DATA GUESP2(  6,1)/       6.0031650D0/
      DATA GUESP3(  6,1)/       1.6422140D0/
      DATA GUESP1(  6,2)/       0.0507330D0/
      DATA GUESP2(  6,2)/       6.0029790D0/
      DATA GUESP3(  6,2)/       0.8924880D0/
      DATA PM3C6 ( 6)/       1.6500000D0/
      DATA PM3R0 ( 6)/       0.1610000D0/
C     DATA FOR ELEMENT  7      NITROGEN
C     DATA REFPM3 ( 7)/'  N: (PM3): J. J. P. STEWART, J. COMP. CHEM.
C    1 10, 209 (1989)                 '/
      DATA USSPM3(  7)/     -49.3356720D0/
      DATA USSPMD(  7)/     -49.3484600D0/
      DATA UPPPM3(  7)/     -47.5097360D0/
      DATA UPPPMD(  7)/     -47.5437680D0/
      DATA BETASP(  7)/     -14.0625210D0/
      DATA BSPMD (  7)/     -14.0684110D0/
      DATA BETAPP(  7)/     -20.0438480D0/
      DATA BPPMD (  7)/     -20.0392920D0/
      DATA ZSPM3 (  7)/       2.0280940D0/
      DATA ZPPM3 (  7)/       2.3137280D0/
      DATA ALPPM3(  7)/       2.8305450D0/
      DATA ALPPMD(  7)/       3.0604040D0/
      DATA EISOLP(  7)/    -157.6137755D0/
      DATA GSSPM3(  7)/      11.9047870D0/
      DATA GSPPM3(  7)/       7.3485650D0/
      DATA GPPPM3(  7)/      11.7546720D0/
      DATA GP2PM3(  7)/      10.8072770D0/
      DATA HSPPM3(  7)/       1.1367130D0/
      DATA DDPM3 (  7)/       0.6577006D0/
      DATA QQPM3 (  7)/       0.5293383D0/
      DATA AMPM3 (  7)/       0.4375151D0/
      DATA ADPM3 (  7)/       0.5030995D0/
      DATA AQPM3 (  7)/       0.7364933D0/
      DATA GUESP1(  7,1)/       1.5016740D0/
      DATA GUESP2(  7,1)/       5.9011480D0/
      DATA GUESP3(  7,1)/       1.7107400D0/
      DATA GUESP1(  7,2)/      -1.5057720D0/
      DATA GUESP2(  7,2)/       6.0046580D0/
      DATA GUESP3(  7,2)/       1.7161490D0/
      DATA PM3C6 ( 7)/       1.1100000D0/
      DATA PM3R0 ( 7)/       0.1550000D0/
C     DATA FOR ELEMENT  8      OXYGEN
C     DATA REFPM3 ( 8)/'  O: (PM3): J. J. P. STEWART, J. COMP. CHEM.
C    1 10, 209 (1989)                 '/
      DATA USSPM3(  8)/     -86.9930020D0/
      DATA USSPMD(  8)/     -86.9603020D0/
      DATA UPPPM3(  8)/     -71.8795800D0/
      DATA UPPPMD(  8)/     -71.9268450D0/
      DATA BETASP(  8)/     -45.2026510D0/
      DATA BSPMD (  8)/     -45.2343020D0/
      DATA BETAPP(  8)/     -24.7525150D0/
      DATA BPPMD (  8)/     -24.7880370D0/
      DATA ZSPM3 (  8)/       3.7965440D0/
      DATA ZPPM3 (  8)/       2.3894020D0/
      DATA ALPPM3(  8)/       3.2171020D0/
      DATA ALPPMD(  8)/       3.3878060D0/
      DATA EISOLP(  8)/    -289.3422065D0/
      DATA GSSPM3(  8)/      15.7557600D0/
      DATA GSPPM3(  8)/      10.6211600D0/
      DATA GPPPM3(  8)/      13.6540160D0/
      DATA GP2PM3(  8)/      12.4060950D0/
      DATA HSPPM3(  8)/       0.5938830D0/
      DATA DDPM3 (  8)/       0.4086173D0/
      DATA QQPM3 (  8)/       0.5125738D0/
      DATA AMPM3 (  8)/       0.5790430D0/
      DATA ADPM3 (  8)/       0.5299630D0/
      DATA AQPM3 (  8)/       0.8179630D0/
      DATA GUESP1(  8,1)/      -1.1311280D0/
      DATA GUESP2(  8,1)/       6.0024770D0/
      DATA GUESP3(  8,1)/       1.6073110D0/
      DATA GUESP1(  8,2)/       1.1378910D0/
      DATA GUESP2(  8,2)/       5.9505120D0/
      DATA GUESP3(  8,2)/       1.5983950D0/
      DATA PM3C6 ( 8)/       0.7000000D0/
      DATA PM3R0 ( 8)/       0.1490000D0/
C     DATA FOR ELEMENT  9      FLUORINE
C     DATA REFPM3 ( 9)/'  F: (PM3): J. J. P. STEWART, J. COMP. CHEM.
C    1 10, 209 (1989)                 '/
      DATA USSPM3(  9)/    -110.4353030D0/
      DATA UPPPM3(  9)/    -105.6850470D0/
      DATA BETASP(  9)/     -48.4059390D0/
      DATA BETAPP(  9)/     -27.7446600D0/
      DATA ZSPM3 (  9)/       4.7085550D0/
      DATA ZPPM3 (  9)/       2.4911780D0/
      DATA ALPPM3(  9)/       3.3589210D0/
      DATA EISOLP(  9)/    -437.5171690D0/
      DATA GSSPM3(  9)/      10.4966670D0/
      DATA GSPPM3(  9)/      16.0736890D0/
      DATA GPPPM3(  9)/      14.8172560D0/
      DATA GP2PM3(  9)/      14.4183930D0/
      DATA HSPPM3(  9)/       0.7277630D0/
      DATA DDPM3 (  9)/       0.3125302D0/
      DATA QQPM3 (  9)/       0.4916328D0/
      DATA AMPM3 (  9)/       0.3857650D0/
      DATA ADPM3 (  9)/       0.6768503D0/
      DATA AQPM3 (  9)/       0.6120047D0/
      DATA GUESP1(  9,1)/      -0.0121660D0/
      DATA GUESP2(  9,1)/       6.0235740D0/
      DATA GUESP3(  9,1)/       1.8568590D0/
      DATA GUESP1(  9,2)/      -0.0028520D0/
      DATA GUESP2(  9,2)/       6.0037170D0/
      DATA GUESP3(  9,2)/       2.6361580D0/
C     DATA FOR ELEMENT 11        SODIUM
C     DATA REFMN  (11)/' Na: (PM3):   E.N. BROTHERS, K.M. MERZ, J.PHYS.C
C    1HEM. B 106, 2779 (2002).        '/
      DATA USSPM3 (11)/      -5.29459290D0/
      DATA UPPPM3 (11)/      -2.45965640D0/
      DATA BETASP (11)/      -0.15108700D0/
      DATA BETAPP (11)/      -0.21840960D0/
      DATA ZSPM3  (11)/       1.13750110D0/
      DATA ZPPM3  (11)/       1.18774330D0/
C     DATA ZDPM3  (11)/       1.00000000D0/
      DATA ALPPM3 (11)/       2.36771690D0/
      DATA EISOLP (11)/      -5.29459290D0/
      DATA GSSPM3 (11)/       3.95586920D0/
      DATA GSPPM3 (11)/       7.19291090D0/
      DATA GPPPM3 (11)/       5.33639630D0/
      DATA GP2PM3 (11)/       5.05880740D0/
      DATA HSPPM3 (11)/       0.56878890D0/
      DATA DDPM3  (11)/       1.73523770D0/
      DATA QQPM3  (11)/       1.40882300D0/
      DATA AMPM3  (11)/       0.14538290D0/
      DATA ADPM3  (11)/       0.21210630D0/
      DATA AQPM3  (11)/       0.25758290D0/
      DATA GUESP1(11,1)/      0.64336550D0/
      DATA GUESP2(11,1)/      1.54650540D0/
      DATA GUESP3(11,1)/      0.99766990D0/
      DATA GUESP1(11,2)/      1.08717880D0/
      DATA GUESP2(11,2)/      1.45290000D0/
      DATA GUESP3(11,2)/      1.45060990D0/
C     DATA FOR ELEMENT 12        MAGNESIUM
C     DATA REFPM3(12)/ ' Mg: (PM3): J. J. P. STEWART, J. COMP. CHEM.
C    1 12, 320 (1991)                 '/
      DATA USSPM3( 12)/     -14.6236880D0/
      DATA UPPPM3( 12)/     -14.1734600D0/
      DATA BETASP( 12)/      -2.0716910D0/
      DATA BETAPP( 12)/      -0.5695810D0/
      DATA ZSPM3 ( 12)/       0.6985520D0/
      DATA ZPPM3 ( 12)/       1.4834530D0/
      DATA ALPPM3( 12)/       1.3291470D0/
      DATA EISOLP( 12)/     -22.5530760D0/
      DATA GSSPM3( 12)/       6.6943000D0/
      DATA GSPPM3( 12)/       6.7939950D0/
      DATA GPPPM3( 12)/       6.9104460D0/
      DATA GP2PM3( 12)/       7.0908230D0/
      DATA HSPPM3( 12)/       0.5433000D0/
      DATA DDPM3 ( 12)/       1.1403950D0/
      DATA QQPM3 ( 12)/       1.1279899D0/
      DATA AMPM3 ( 12)/       0.2460235D0/
      DATA ADPM3 ( 12)/       0.2695751D0/
      DATA AQPM3 ( 12)/       0.2767522D0/
      DATA GUESP1( 12,1)/       2.1170500D0/
      DATA GUESP2( 12,1)/       6.0094770D0/
      DATA GUESP3( 12,1)/       2.0844060D0/
      DATA GUESP1( 12,2)/      -2.5477670D0/
      DATA GUESP2( 12,2)/       4.3953700D0/
      DATA GUESP3( 12,2)/       2.0636740D0/
C     DATA FOR ELEMENT 13      ALUMINUM
C     DATA REFPM3 (13)/' Al: (PM3): J. J. P. STEWART, J. COMP. CHEM.
C    1 10, 209 (1989)                 '/
      DATA USSPM3( 13)/     -24.8454040D0/
      DATA UPPPM3( 13)/     -22.2641590D0/
      DATA BETASP( 13)/      -0.5943010D0/
      DATA BETAPP( 13)/      -0.9565500D0/
      DATA ZSPM3 ( 13)/       1.7028880D0/
      DATA ZPPM3 ( 13)/       1.0736290D0/
      DATA ALPPM3( 13)/       1.5217030D0/
      DATA EISOLP( 13)/     -46.8647630D0/
      DATA GSSPM3( 13)/       5.7767370D0/
      DATA GSPPM3( 13)/      11.6598560D0/
      DATA GPPPM3( 13)/       6.3477900D0/
      DATA GP2PM3( 13)/       6.1210770D0/
      DATA HSPPM3( 13)/       4.0062450D0/
      DATA DDPM3 ( 13)/       1.2102799D0/
      DATA QQPM3 ( 13)/       1.5585645D0/
      DATA AMPM3 ( 13)/       0.2123020D0/
      DATA ADPM3 ( 13)/       0.6418584D0/
      DATA AQPM3 ( 13)/       0.2262838D0/
      DATA GUESP1( 13,1)/      -0.4730900D0/
      DATA GUESP2( 13,1)/       1.9158250D0/
      DATA GUESP3( 13,1)/       1.4517280D0/
      DATA GUESP1( 13,2)/      -0.1540510D0/
      DATA GUESP2( 13,2)/       6.0050860D0/
      DATA GUESP3( 13,2)/       2.5199970D0/
C     DATA FOR ELEMENT 14      SILICON
C     DATA REFPM3 (14)/' Si: (PM3): J. J. P. STEWART, J. COMP. CHEM.
C    1 10, 209 (1989)                 '/
C     CORRECTION SEE QCPE BULLETIN 10, 34 (1990)
      DATA USSPM3( 14)/     -26.7634830D0/
      DATA UPPPM3( 14)/     -22.8136350D0/
      DATA BETASP( 14)/      -2.8621450D0/
      DATA BETAPP( 14)/      -3.9331480D0/
      DATA ZSPM3 ( 14)/       1.6350750D0/
      DATA ZPPM3 ( 14)/       1.3130880D0/
      DATA ALPPM3( 14)/       2.1358090D0/
      DATA EISOLP( 14)/     -67.7882140D0/
      DATA GSSPM3( 14)/       5.0471960D0/
      DATA GSPPM3( 14)/       5.9490570D0/
      DATA GPPPM3( 14)/       6.7593670D0/
      DATA GP2PM3( 14)/       5.1612970D0/
      DATA HSPPM3( 14)/       0.9198320D0/
      DATA DDPM3 ( 14)/       1.3144550D0/
      DATA QQPM3 ( 14)/       1.2743396D0/
      DATA AMPM3 ( 14)/       0.1854905D0/
      DATA ADPM3 ( 14)/       0.3060715D0/
      DATA AQPM3 ( 14)/       0.4877432D0/
      DATA GUESP1( 14,1)/      -0.3906000D0/
      DATA GUESP2( 14,1)/       6.0000540D0/
      DATA GUESP3( 14,1)/       0.6322620D0/
      DATA GUESP1( 14,2)/       0.0572590D0/
      DATA GUESP2( 14,2)/       6.0071830D0/
      DATA GUESP3( 14,2)/       2.0199870D0/
COLD  DATA GUESP1( 14,3)/       0.0207890D0/
COLD  DATA GUESP2( 14,3)/       5.0000000D0/
COLD  DATA GUESP3( 14,3)/       2.9906100D0/
C     DATA FOR ELEMENT 15      PHOSPHORUS
C     DATA REFPM3 (15)/'  P: (PM3): J. J. P. STEWART, J. COMP. CHEM.
C    1 10, 209 (1989)                 '/
      DATA USSPM3( 15)/     -40.4130960D0/
      DATA UPPPM3( 15)/     -29.5930520D0/
      DATA BETASP( 15)/     -12.6158790D0/
      DATA BETAPP( 15)/      -4.1600400D0/
      DATA ZSPM3 ( 15)/       2.0175630D0/
      DATA ZPPM3 ( 15)/       1.5047320D0/
      DATA ALPPM3( 15)/       1.9405340D0/
      DATA EISOLP( 15)/    -117.9591740D0/
      DATA GSSPM3( 15)/       7.8016150D0/
      DATA GSPPM3( 15)/       5.1869490D0/
      DATA GPPPM3( 15)/       6.6184780D0/
      DATA GP2PM3( 15)/       6.0620020D0/
      DATA HSPPM3( 15)/       1.5428090D0/
      DATA DDPM3 ( 15)/       1.0644947D0/
      DATA QQPM3 ( 15)/       1.1120386D0/
      DATA AMPM3 ( 15)/       0.2867187D0/
      DATA ADPM3 ( 15)/       0.4309446D0/
      DATA AQPM3 ( 15)/       0.3732517D0/
      DATA GUESP1( 15,1)/      -0.6114210D0/
      DATA GUESP2( 15,1)/       1.9972720D0/
      DATA GUESP3( 15,1)/       0.7946240D0/
      DATA GUESP1( 15,2)/      -0.0939350D0/
      DATA GUESP2( 15,2)/       1.9983600D0/
      DATA GUESP3( 15,2)/       1.9106770D0/
C     DATA FOR ELEMENT 16      SULFUR
C     DATA REFPM3 (16)/'  S: (PM3): J. J. P. STEWART, J. COMP. CHEM.
C    1 10, 209 (1989)                 '/
      DATA USSPM3( 16)/     -49.8953710D0/
      DATA UPPPM3( 16)/     -44.3925830D0/
      DATA BETASP( 16)/      -8.8274650D0/
      DATA BETAPP( 16)/      -8.0914150D0/
      DATA ZSPM3 ( 16)/       1.8911850D0/
      DATA ZPPM3 ( 16)/       1.6589720D0/
      DATA ALPPM3( 16)/       2.2697060D0/
      DATA EISOLP( 16)/    -183.4537395D0/
      DATA GSSPM3( 16)/       8.9646670D0/
      DATA GSPPM3( 16)/       6.7859360D0/
      DATA GPPPM3( 16)/       9.9681640D0/
      DATA GP2PM3( 16)/       7.9702470D0/
      DATA HSPPM3( 16)/       4.0418360D0/
      DATA DDPM3 ( 16)/       1.1214313D0/
      DATA QQPM3 ( 16)/       1.0086488D0/
      DATA AMPM3 ( 16)/       0.3294622D0/
      DATA ADPM3 ( 16)/       0.6679118D0/
      DATA AQPM3 ( 16)/       0.6137472D0/
      DATA GUESP1( 16,1)/      -0.3991910D0/
      DATA GUESP2( 16,1)/       6.0006690D0/
      DATA GUESP3( 16,1)/       0.9621230D0/
      DATA GUESP1( 16,2)/      -0.0548990D0/
      DATA GUESP2( 16,2)/       6.0018450D0/
      DATA GUESP3( 16,2)/       1.5799440D0/
C     DATA FOR ELEMENT 17      CHLORINE
C     DATA REFPM3 (17)/' Cl: (PM3): J. J. P. STEWART, J. COMP. CHEM.
C    1 10, 209 (1989)                 '/
      DATA USSPM3( 17)/    -100.6267470D0/
      DATA UPPPM3( 17)/     -53.6143960D0/
      DATA BETASP( 17)/     -27.5285600D0/
      DATA BETAPP( 17)/     -11.5939220D0/
      DATA ZSPM3 ( 17)/       2.2462100D0/
      DATA ZPPM3 ( 17)/       2.1510100D0/
      DATA ALPPM3( 17)/       2.5172960D0/
      DATA EISOLP( 17)/    -315.1949480D0/
      DATA GSSPM3( 17)/      16.0136010D0/
      DATA GSPPM3( 17)/       8.0481150D0/
      DATA GPPPM3( 17)/       7.5222150D0/
      DATA GP2PM3( 17)/       7.5041540D0/
      DATA HSPPM3( 17)/       3.4811530D0/
      DATA DDPM3 ( 17)/       0.9175856D0/
      DATA QQPM3 ( 17)/       0.7779230D0/
      DATA AMPM3 ( 17)/       0.5885190D0/
      DATA ADPM3 ( 17)/       0.6814522D0/
      DATA AQPM3 ( 17)/       0.3643694D0/
      DATA GUESP1( 17,1)/      -0.1715910D0/
      DATA GUESP2( 17,1)/       6.0008020D0/
      DATA GUESP3( 17,1)/       1.0875020D0/
      DATA GUESP1( 17,2)/      -0.0134580D0/
      DATA GUESP2( 17,2)/       1.9666180D0/
      DATA GUESP3( 17,2)/       2.2928910D0/
C     DATA FOR ELEMENT 30        ZINC
C     DATA REFPM3(30)/ ' Zn: (PM3): J. J. P. STEWART, J. COMP. CHEM.
C    1 12, 320 (1991)                 '/
      DATA USSPM3( 30)/     -18.5321980D0/
      DATA UPPPM3( 30)/     -11.0474090D0/
      DATA BETASP( 30)/      -0.7155780D0/
      DATA BETAPP( 30)/      -6.3518640D0/
      DATA ZSPM3 ( 30)/       1.8199890D0/
      DATA ZPPM3 ( 30)/       1.5069220D0/
      DATA ALPPM3( 30)/       1.3501260D0/
      DATA EISOLP( 30)/     -27.3872000D0/
      DATA GSSPM3( 30)/       9.6771960D0/
      DATA GSPPM3( 30)/       7.7362040D0/
      DATA GPPPM3( 30)/       4.9801740D0/
      DATA GP2PM3( 30)/       4.6696560D0/
      DATA HSPPM3( 30)/       0.6004130D0/
      DATA DDPM3 ( 30)/       1.5005758D0/
      DATA QQPM3 ( 30)/       1.4077174D0/
      DATA AMPM3 ( 30)/       0.3556485D0/
      DATA ADPM3 ( 30)/       0.2375689D0/
      DATA AQPM3 ( 30)/       0.2661069D0/
      DATA GUESP1( 30,1)/      -0.1112340D0/
      DATA GUESP2( 30,1)/       6.0014780D0/
      DATA GUESP3( 30,1)/       1.5160320D0/
      DATA GUESP1( 30,2)/      -0.1323700D0/
      DATA GUESP2( 30,2)/       1.9958390D0/
      DATA GUESP3( 30,2)/       2.5196420D0/
C     DATA FOR ELEMENT 31        GALLIUM
C     DATA REFPM3(31)/ ' Ga: (PM3): J. J. P. STEWART, J. COMP. CHEM.
C    1 12, 320 (1991)                 '/
      DATA USSPM3( 31)/     -29.8555930D0/
      DATA UPPPM3( 31)/     -21.8753710D0/
      DATA BETASP( 31)/      -4.9456180D0/
      DATA BETAPP( 31)/      -0.4070530D0/
      DATA ZSPM3 ( 31)/       1.8470400D0/
      DATA ZPPM3 ( 31)/       0.8394110D0/
      DATA ALPPM3( 31)/       1.6051150D0/
      DATA EISOLP( 31)/     -57.3280250D0/
      DATA GSSPM3( 31)/       8.4585540D0/
      DATA GSPPM3( 31)/       8.9256190D0/
      DATA GPPPM3( 31)/       5.0868550D0/
      DATA GP2PM3( 31)/       4.9830450D0/
      DATA HSPPM3( 31)/       2.0512600D0/
      DATA DDPM3 ( 31)/       0.9776692D0/
      DATA QQPM3 ( 31)/       2.5271534D0/
      DATA AMPM3 ( 31)/       0.3108620D0/
      DATA ADPM3 ( 31)/       0.5129360D0/
      DATA AQPM3 ( 31)/       0.1546208D0/
      DATA GUESP1( 31,1)/      -0.5601790D0/
      DATA GUESP2( 31,1)/       5.6232730D0/
      DATA GUESP3( 31,1)/       1.5317800D0/
      DATA GUESP1( 31,2)/      -0.2727310D0/
      DATA GUESP2( 31,2)/       1.9918430D0/
      DATA GUESP3( 31,2)/       2.1838640D0/
C     DATA FOR ELEMENT 32        GERMANIUM
C     DATA REFPM3(32)/ ' Ge: (PM3): J. J. P. STEWART, J. COMP. CHEM.
C    1 12, 320 (1991)                 '/
      DATA USSPM3( 32)/     -35.4671955D0/
      DATA UPPPM3( 32)/     -31.5863583D0/
      DATA BETASP( 32)/      -5.3250024D0/
      DATA BETAPP( 32)/      -2.2501567D0/
      DATA ZSPM3 ( 32)/       2.2373526D0/
      DATA ZPPM3 ( 32)/       1.5924319D0/
      DATA ALPPM3( 32)/       1.9723370D0/
      DATA EISOLP( 32)/     -84.0156006D0/
      DATA GSSPM3( 32)/       5.3769635D0/
      DATA GSPPM3( 32)/      10.2095293D0/
      DATA GPPPM3( 32)/       7.6718647D0/
      DATA GP2PM3( 32)/       6.9242663D0/
      DATA HSPPM3( 32)/       1.3370204D0/
      DATA DDPM3 ( 32)/       1.1920304D0/
      DATA QQPM3 ( 32)/       1.3321263D0/
      DATA AMPM3 ( 32)/       0.1976098D0/
      DATA ADPM3 ( 32)/       0.3798182D0/
      DATA AQPM3 ( 32)/       0.3620669D0/
      DATA GUESP1( 32,1)/       0.9631726D0/
      DATA GUESP2( 32,1)/       6.0120134D0/
      DATA GUESP3( 32,1)/       2.1633655D0/
      DATA GUESP1( 32,2)/      -0.9593891D0/
      DATA GUESP2( 32,2)/       5.7491802D0/
      DATA GUESP3( 32,2)/       2.1693724D0/
C     DATA FOR ELEMENT 33        ARSENIC
C     DATA REFPM3(33)/ ' As: (PM3): J. J. P. STEWART, J. COMP. CHEM.
C    1 12, 320 (1991)                 '/
      DATA USSPM3( 33)/     -38.5074240D0/
      DATA UPPPM3( 33)/     -35.1524150D0/
      DATA BETASP( 33)/      -8.2321650D0/
      DATA BETAPP( 33)/      -5.0173860D0/
      DATA ZSPM3 ( 33)/       2.6361770D0/
      DATA ZPPM3 ( 33)/       1.7038890D0/
      DATA ALPPM3( 33)/       1.7944770D0/
      DATA EISOLP( 33)/    -122.6326140D0/
      DATA GSSPM3( 33)/       8.7890010D0/
      DATA GSPPM3( 33)/       5.3979830D0/
      DATA GPPPM3( 33)/       8.2872500D0/
      DATA GP2PM3( 33)/       8.2103460D0/
      DATA HSPPM3( 33)/       1.9510340D0/
      DATA DDPM3 ( 33)/       0.9679655D0/
      DATA QQPM3 ( 33)/       1.2449874D0/
      DATA AMPM3 ( 33)/       0.3230063D0/
      DATA ADPM3 ( 33)/       0.5042239D0/
      DATA AQPM3 ( 33)/       0.2574219D0/
      DATA GUESP1( 33,1)/      -0.4600950D0/
      DATA GUESP2( 33,1)/       1.9831150D0/
      DATA GUESP3( 33,1)/       1.0867930D0/
      DATA GUESP1( 33,2)/      -0.0889960D0/
      DATA GUESP2( 33,2)/       1.9929440D0/
      DATA GUESP3( 33,2)/       2.1400580D0/
C     DATA FOR ELEMENT 34        SELENIUM
C     DATA REFPM3(34)/ ' Se: (PM3): J. J. P. STEWART, J. COMP. CHEM.
C    1 12, 320 (1991)                 '/
      DATA USSPM3( 34)/     -55.3781350D0/
      DATA UPPPM3( 34)/     -49.8230760D0/
      DATA BETASP( 34)/      -6.1578220D0/
      DATA BETAPP( 34)/      -5.4930390D0/
      DATA ZSPM3 ( 34)/       2.8280510D0/
      DATA ZPPM3 ( 34)/       1.7325360D0/
      DATA ALPPM3( 34)/       3.0439570D0/
      DATA EISOLP( 34)/    -192.7748115D0/
      DATA GSSPM3( 34)/       7.4325910D0/
      DATA GSPPM3( 34)/      10.0604610D0/
      DATA GPPPM3( 34)/       9.5683260D0/
      DATA GP2PM3( 34)/       7.7242890D0/
      DATA HSPPM3( 34)/       4.0165580D0/
      DATA DDPM3 ( 34)/       0.8719813D0/
      DATA QQPM3 ( 34)/       1.2244019D0/
      DATA AMPM3 ( 34)/       0.2731566D0/
      DATA ADPM3 ( 34)/       0.7509697D0/
      DATA AQPM3 ( 34)/       0.5283737D0/
      DATA GUESP1( 34,1)/       0.0478730D0/
      DATA GUESP2( 34,1)/       6.0074000D0/
      DATA GUESP3( 34,1)/       2.0817170D0/
      DATA GUESP1( 34,2)/       0.1147200D0/
      DATA GUESP2( 34,2)/       6.0086720D0/
      DATA GUESP3( 34,2)/       1.5164230D0/
C     DATA FOR ELEMENT 35      BROMINE
C     DATA REFPM3 (35)/' Br: (PM3): J. J. P. STEWART, J. COMP. CHEM.
C    1 10, 209 (1989)                 '/
      DATA USSPM3( 35)/    -116.6193110D0/
      DATA UPPPM3( 35)/     -74.2271290D0/
      DATA BETASP( 35)/     -31.1713420D0/
      DATA BETAPP( 35)/      -6.8140130D0/
      DATA ZSPM3 ( 35)/       5.3484570D0/
      DATA ZPPM3 ( 35)/       2.1275900D0/
      DATA ALPPM3( 35)/       2.5118420D0/
      DATA EISOLP( 35)/    -352.5398970D0/
      DATA GSSPM3( 35)/      15.9434250D0/
      DATA GSPPM3( 35)/      16.0616800D0/
      DATA GPPPM3( 35)/       8.2827630D0/
      DATA GP2PM3( 35)/       7.8168490D0/
      DATA HSPPM3( 35)/       0.5788690D0/
      DATA DDPM3 ( 35)/       0.2759025D0/
      DATA QQPM3 ( 35)/       0.9970532D0/
      DATA AMPM3 ( 35)/       0.5859399D0/
      DATA ADPM3 ( 35)/       0.6755383D0/
      DATA AQPM3 ( 35)/       0.3823719D0/
      DATA GUESP1( 35,1)/       0.9604580D0/
      DATA GUESP2( 35,1)/       5.9765080D0/
      DATA GUESP3( 35,1)/       2.3216540D0/
      DATA GUESP1( 35,2)/      -0.9549160D0/
      DATA GUESP2( 35,2)/       5.9447030D0/
      DATA GUESP3( 35,2)/       2.3281420D0/
C     DATA FOR ELEMENT 48        CADMIUM
C     DATA REFPM3(48)/ ' Cd: (PM3): J. J. P. STEWART, J. COMP. CHEM.
C    1 12, 320 (1991)                 '/
      DATA USSPM3( 48)/     -15.8285840D0/
      DATA UPPPM3( 48)/       8.7497950D0/
      DATA BETASP( 48)/      -8.5819440D0/
      DATA BETAPP( 48)/      -0.6010340D0/
      DATA ZSPM3 ( 48)/       1.6793510D0/
      DATA ZPPM3 ( 48)/       2.0664120D0/
      DATA ALPPM3( 48)/       1.5253820D0/
      DATA EISOLP( 48)/     -22.4502080D0/
      DATA GSSPM3( 48)/       9.2069600D0/
      DATA GSPPM3( 48)/       8.2315390D0/
      DATA GPPPM3( 48)/       4.9481040D0/
      DATA GP2PM3( 48)/       4.6696560D0/
      DATA HSPPM3( 48)/       1.6562340D0/
      DATA DDPM3 ( 48)/       1.5982681D0/
      DATA QQPM3 ( 48)/       1.2432402D0/
      DATA AMPM3 ( 48)/       0.3383668D0/
      DATA ADPM3 ( 48)/       0.3570290D0/
      DATA AQPM3 ( 48)/       0.2820582D0/
      DATA GUESP1( 48,1)/     0.0000000D0/
      DATA GUESP2( 48,1)/     0.0000000D0/
      DATA GUESP3( 48,1)/     0.0000000D0/
C     DATA FOR ELEMENT 49        INDIUM
C     DATA REFPM3(49)/ ' In: (PM3): J. J. P. STEWART, J. COMP. CHEM.
C    1 12, 320 (1991)                 '/
      DATA USSPM3( 49)/     -26.1762050D0/
      DATA UPPPM3( 49)/     -20.0058220D0/
      DATA BETASP( 49)/      -2.9933190D0/
      DATA BETAPP( 49)/      -1.8289080D0/
      DATA ZSPM3 ( 49)/       2.0161160D0/
      DATA ZPPM3 ( 49)/       1.4453500D0/
      DATA ALPPM3( 49)/       1.4183850D0/
      DATA EISOLP( 49)/     -51.9750470D0/
      DATA GSSPM3( 49)/       6.5549000D0/
      DATA GSPPM3( 49)/       8.2298730D0/
      DATA GPPPM3( 49)/       6.2992690D0/
      DATA GP2PM3( 49)/       4.9842110D0/
      DATA HSPPM3( 49)/       2.6314610D0/
      DATA DDPM3 ( 49)/       1.5766241D0/
      DATA QQPM3 ( 49)/       1.7774563D0/
      DATA AMPM3 ( 49)/       0.2409004D0/
      DATA ADPM3 ( 49)/       0.4532655D0/
      DATA AQPM3 ( 49)/       0.3689812D0/
      DATA GUESP1( 49,1)/      -0.3431380D0/
      DATA GUESP2( 49,1)/       1.9940340D0/
      DATA GUESP3( 49,1)/       1.6255160D0/
      DATA GUESP1( 49,2)/      -0.1095320D0/
      DATA GUESP2( 49,2)/       5.6832170D0/
      DATA GUESP3( 49,2)/       2.8670090D0/
C     DATA FOR ELEMENT 50        TIN
C     DATA REFPM3(50)/ ' Sn: (PM3): J. J. P. STEWART, J. COMP. CHEM.
C    1 12, 320 (1991)                 '/
      DATA USSPM3( 50)/     -34.5501920D0/
      DATA UPPPM3( 50)/     -25.8944190D0/
      DATA BETASP( 50)/      -2.7858020D0/
      DATA BETAPP( 50)/      -2.0059990D0/
      DATA ZSPM3 ( 50)/       2.3733280D0/
      DATA ZPPM3 ( 50)/       1.6382330D0/
      DATA ALPPM3( 50)/       1.6996500D0/
      DATA EISOLP( 50)/     -78.8877790D0/
      DATA GSSPM3( 50)/      10.1900330D0/
      DATA GSPPM3( 50)/       7.2353270D0/
      DATA GPPPM3( 50)/       5.6738100D0/
      DATA GP2PM3( 50)/       5.1822140D0/
      DATA HSPPM3( 50)/       1.0331570D0/
      DATA DDPM3 ( 50)/       1.3120038D0/
      DATA QQPM3 ( 50)/       1.5681814D0/
      DATA AMPM3 ( 50)/       0.3744959D0/
      DATA ADPM3 ( 50)/       0.3218163D0/
      DATA AQPM3 ( 50)/       0.2832529D0/
      DATA GUESP1( 50,1)/      -0.1503530D0/
      DATA GUESP2( 50,1)/       6.0056940D0/
      DATA GUESP3( 50,1)/       1.7046420D0/
      DATA GUESP1( 50,2)/      -0.0444170D0/
      DATA GUESP2( 50,2)/       2.2573810D0/
      DATA GUESP3( 50,2)/       2.4698690D0/
C     DATA FOR ELEMENT 51        ANTIMONY
C     DATA REFPM3(51)/ ' Sb: (PM3): J. J. P. STEWART, J. COMP. CHEM.
C    1 12, 320 (1991)                 '/
      DATA USSPM3( 51)/     -56.4321960D0/
      DATA UPPPM3( 51)/     -29.4349540D0/
      DATA BETASP( 51)/     -14.7942170D0/
      DATA BETAPP( 51)/      -2.8179480D0/
      DATA ZSPM3 ( 51)/       2.3430390D0/
      DATA ZPPM3 ( 51)/       1.8999920D0/
      DATA ALPPM3( 51)/       2.0343010D0/
      DATA EISOLP( 51)/    -148.9382890D0/
      DATA GSSPM3( 51)/       9.2382770D0/
      DATA GSPPM3( 51)/       5.2776800D0/
      DATA GPPPM3( 51)/       6.3500000D0/
      DATA GP2PM3( 51)/       6.2500000D0/
      DATA HSPPM3( 51)/       2.4244640D0/
      DATA DDPM3 ( 51)/       1.4091903D0/
      DATA QQPM3 ( 51)/       1.3521354D0/
      DATA AMPM3 ( 51)/       0.3395177D0/
      DATA ADPM3 ( 51)/       0.4589010D0/
      DATA AQPM3 ( 51)/       0.2423472D0/
      DATA GUESP1( 51,1)/       3.0020280D0/
      DATA GUESP2( 51,1)/       6.0053420D0/
      DATA GUESP3( 51,1)/       0.8530600D0/
      DATA GUESP1( 51,2)/      -0.0188920D0/
      DATA GUESP2( 51,2)/       6.0114780D0/
      DATA GUESP3( 51,2)/       2.7933110D0/
C     DATA FOR ELEMENT 52        TELLURIUM
C     DATA REFPM3(52)/ ' Te: (PM3): J. J. P. STEWART, J. COMP. CHEM.
C    1 12, 320 (1991)                 '/
      DATA USSPM3( 52)/     -44.9380360D0/
      DATA UPPPM3( 52)/     -46.3140990D0/
      DATA BETASP( 52)/      -2.6651460D0/
      DATA BETAPP( 52)/      -3.8954300D0/
      DATA ZSPM3 ( 52)/       4.1654920D0/
      DATA ZPPM3 ( 52)/       1.6475550D0/
      DATA ALPPM3( 52)/       2.4850190D0/
      DATA EISOLP( 52)/    -168.0945925D0/
      DATA GSSPM3( 52)/      10.2550730D0/
      DATA GSPPM3( 52)/       8.1691450D0/
      DATA GPPPM3( 52)/       7.7775920D0/
      DATA GP2PM3( 52)/       7.7551210D0/
      DATA HSPPM3( 52)/       3.7724620D0/
      DATA DDPM3 ( 52)/       0.3484177D0/
      DATA QQPM3 ( 52)/       1.5593085D0/
      DATA AMPM3 ( 52)/       0.3768862D0/
      DATA ADPM3 ( 52)/       1.1960743D0/
      DATA AQPM3 ( 52)/       0.2184786D0/
      DATA GUESP1( 52,1)/       0.0333910D0/
      DATA GUESP2( 52,1)/       5.9563790D0/
      DATA GUESP3( 52,1)/       2.2775750D0/
      DATA GUESP1( 52,2)/      -1.9218670D0/
      DATA GUESP2( 52,2)/       4.9732190D0/
      DATA GUESP3( 52,2)/       0.5242430D0/
C     DATA FOR ELEMENT 53      IODINE
C     DATA REFPM3 (53)/'  I: (PM3): J. J. P. STEWART, J. COMP. CHEM.
C    1 10, 209 (1989)                 '/
      DATA USSPM3( 53)/     -96.4540370D0/
      DATA UPPPM3( 53)/     -61.0915820D0/
      DATA BETASP( 53)/     -14.4942340D0/
      DATA BETAPP( 53)/      -5.8947030D0/
      DATA ZSPM3 ( 53)/       7.0010130D0/
      DATA ZPPM3 ( 53)/       2.4543540D0/
      DATA ALPPM3( 53)/       1.9901850D0/
      DATA EISOLP( 53)/    -288.3160860D0/
      DATA GSSPM3( 53)/      13.6319430D0/
      DATA GSPPM3( 53)/      14.9904060D0/
      DATA GPPPM3( 53)/       7.2883300D0/
      DATA GP2PM3( 53)/       5.9664070D0/
      DATA HSPPM3( 53)/       2.6300350D0/
      DATA DDPM3 ( 53)/       0.1581469D0/
      DATA QQPM3 ( 53)/       1.0467302D0/
      DATA AMPM3 ( 53)/       0.5009902D0/
      DATA ADPM3 ( 53)/       1.6699104D0/
      DATA AQPM3 ( 53)/       0.5153082D0/
      DATA GUESP1( 53,1)/      -0.1314810D0/
      DATA GUESP2( 53,1)/       5.2064170D0/
      DATA GUESP3( 53,1)/       1.7488240D0/
      DATA GUESP1( 53,2)/      -0.0368970D0/
      DATA GUESP2( 53,2)/       6.0101170D0/
      DATA GUESP3( 53,2)/       2.7103730D0/
C     DATA FOR ELEMENT 80        MERCURY
C     DATA REFPM3(80)/ ' Hg: (PM3): J. J. P. STEWART, J. COMP. CHEM.
C    1 12, 320 (1991)                 '/
      DATA USSPM3( 80)/     -17.7622290D0/
      DATA UPPPM3( 80)/     -18.3307510D0/
      DATA BETASP( 80)/      -3.1013650D0/
      DATA BETAPP( 80)/      -3.4640310D0/
      DATA ZSPM3 ( 80)/       1.4768850D0/
      DATA ZPPM3 ( 80)/       2.4799510D0/
      DATA ALPPM3( 80)/       1.5293770D0/
      DATA EISOLP( 80)/     -28.8997380D0/
      DATA GSSPM3( 80)/       6.6247200D0/
      DATA GSPPM3( 80)/      10.6392970D0/
      DATA GPPPM3( 80)/      14.7092830D0/
      DATA GP2PM3( 80)/      16.0007400D0/
      DATA HSPPM3( 80)/       2.0363110D0/
      DATA DDPM3 ( 80)/       1.2317811D0/
      DATA QQPM3 ( 80)/       1.2164033D0/
      DATA AMPM3 ( 80)/       0.2434664D0/
      DATA ADPM3 ( 80)/       0.4515472D0/
      DATA AQPM3 ( 80)/       0.2618394D0/
      DATA GUESP1( 80,1)/       1.0827200D0/
      DATA GUESP2( 80,1)/       6.4965980D0/
      DATA GUESP3( 80,1)/       1.1951460D0/
      DATA GUESP1( 80,2)/      -0.0965530D0/
      DATA GUESP2( 80,2)/       3.9262810D0/
      DATA GUESP3( 80,2)/       2.6271600D0/
C     DATA FOR ELEMENT 81        THALLIUM
C     DATA REFPM3(81)/ ' Tl: (PM3): J. J. P. STEWART, J. COMP. CHEM.
C    1 12, 320 (1991)                 '/
      DATA USSPM3( 81)/     -30.0531700D0/
      DATA UPPPM3( 81)/     -26.9206370D0/
      DATA BETASP( 81)/      -1.0844950D0/
      DATA BETAPP( 81)/      -7.9467990D0/
      DATA ZSPM3 ( 81)/       6.8679210D0/
      DATA ZPPM3 ( 81)/       1.9694450D0/
      DATA ALPPM3( 81)/       1.3409510D0/
      DATA EISOLP( 81)/     -56.6492050D0/
      DATA GSSPM3( 81)/      10.4604120D0/
      DATA GSPPM3( 81)/      11.2238830D0/
      DATA GPPPM3( 81)/       4.9927850D0/
      DATA GP2PM3( 81)/       8.9627270D0/
      DATA HSPPM3( 81)/       2.5304060D0/
      DATA DDPM3 ( 81)/       0.0781362D0/
      DATA QQPM3 ( 81)/       1.5317110D0/
      DATA AMPM3 ( 81)/       0.3844326D0/
      DATA ADPM3 ( 81)/       2.5741815D0/
      DATA AQPM3 ( 81)/       0.2213264D0/
      DATA GUESP1( 81,1)/      -1.3613990D0/
      DATA GUESP2( 81,1)/       3.5572260D0/
      DATA GUESP3( 81,1)/       1.0928020D0/
      DATA GUESP1( 81,2)/      -0.0454010D0/
      DATA GUESP2( 81,2)/       2.3069950D0/
      DATA GUESP3( 81,2)/       2.9650290D0/
C     DATA FOR ELEMENT 82        LEAD
C     DATA REFPM3(82)/ ' Pb: (PM3): J. J. P. STEWART, J. COMP. CHEM.
C    1 12, 320 (1991)                 '/
      DATA USSPM3( 82)/     -30.3227560D0/
      DATA UPPPM3( 82)/     -24.4258340D0/
      DATA BETASP( 82)/      -6.1260240D0/
      DATA BETAPP( 82)/      -1.3954300D0/
      DATA ZSPM3 ( 82)/       3.1412890D0/
      DATA ZPPM3 ( 82)/       1.8924180D0/
      DATA ALPPM3( 82)/       1.6200450D0/
      DATA EISOLP( 82)/     -73.4660775D0/
      DATA GSSPM3( 82)/       7.0119920D0/
      DATA GSPPM3( 82)/       6.7937820D0/
      DATA GPPPM3( 82)/       5.1837800D0/
      DATA GP2PM3( 82)/       5.0456510D0/
      DATA HSPPM3( 82)/       1.5663020D0/
      DATA DDPM3 ( 82)/       0.9866290D0/
      DATA QQPM3 ( 82)/       1.5940562D0/
      DATA AMPM3 ( 82)/       0.2576991D0/
      DATA ADPM3 ( 82)/       0.4527678D0/
      DATA AQPM3 ( 82)/       0.2150175D0/
      DATA GUESP1( 82,1)/      -0.1225760D0/
      DATA GUESP2( 82,1)/       6.0030620D0/
      DATA GUESP3( 82,1)/       1.9015970D0/
      DATA GUESP1( 82,2)/      -0.0566480D0/
      DATA GUESP2( 82,2)/       4.7437050D0/
      DATA GUESP3( 82,2)/       2.8618790D0/
C     DATA FOR ELEMENT 83        BISMUTH
C     DATA REFPM3(83)/ ' Bi: (PM3): J. J. P. STEWART, J. COMP. CHEM.
C    1 12, 320 (1991)                 '/
      DATA USSPM3( 83)/     -33.4959380D0/
      DATA UPPPM3( 83)/     -35.5210260D0/
      DATA BETASP( 83)/      -5.6072830D0/
      DATA BETAPP( 83)/      -5.8001520D0/
      DATA ZSPM3 ( 83)/       4.9164510D0/
      DATA ZPPM3 ( 83)/       1.9349350D0/
      DATA ALPPM3( 83)/       1.8574310D0/
      DATA EISOLP( 83)/    -109.2774910D0/
      DATA GSSPM3( 83)/       4.9894800D0/
      DATA GSPPM3( 83)/       6.1033080D0/
      DATA GPPPM3( 83)/       8.6960070D0/
      DATA GP2PM3( 83)/       8.3354470D0/
      DATA HSPPM3( 83)/       0.5991220D0/
      DATA DDPM3 ( 83)/       0.2798609D0/
      DATA QQPM3 ( 83)/       1.5590294D0/
      DATA AMPM3 ( 83)/       0.1833693D0/
      DATA ADPM3 ( 83)/       0.6776013D0/
      DATA AQPM3 ( 83)/       0.2586520D0/
      DATA GUESP1( 83,1)/       2.5816930D0/
      DATA GUESP2( 83,1)/       5.0940220D0/
      DATA GUESP3( 83,1)/       0.4997870D0/
      DATA GUESP1( 83,2)/       0.0603200D0/
      DATA GUESP2( 83,2)/       6.0015380D0/
      DATA GUESP3( 83,2)/       2.4279700D0/
C     DATA FOR ELEMENT 86        CONNECTION ATOM FOR QM/MM TREATMENTS
C     DATA REFPM3 (86)/' **: (PM3): I. ANTES, W. THIEL, J. PHYS. CHEM. A
C    1 103, 9290-9295 (1999)          '/
      DATA USSPM3(86)/     -11.59537348D0/
      DATA UPPPM3(86)/       0.00000000D0/
      DATA BETASP(86)/      -9.94923778D0/
      DATA BETAPP(86)/      -9.94923778D0/
      DATA ZSPM3 (86)/       1.43824357D0/
      DATA ZPPM3 (86)/       1.43824357D0/
      DATA ALPPM3(86)/       1.96719427D0/
      DATA EISOLP(86)/     -11.59537348D0/
      DATA GSSPM3(86)/      11.60459808D0/
      DATA GSPPM3(86)/       0.00000000D0/
      DATA GPPPM3(86)/       0.00000000D0/
      DATA GP2PM3(86)/       0.00000000D0/
      DATA HSPPM3(86)/       0.00000000D0/
      DATA DDPM3 (86)/       0.00000000D0/
      DATA QQPM3 (86)/       0.00000000D0/
      DATA AMPM3 (86)/       0.42648284D0/
      DATA ADPM3 (86)/       0.42648284D0/
      DATA AQPM3 (86)/       0.42648284D0/
      DATA GUESP1(86,1)/    -4.25077548D0/
      DATA GUESP2(86,1)/     6.90952417D0/
      DATA GUESP3(86,1)/     1.51754243D0/
      DATA GUESP1(86,2)/     4.40160107D0/
      DATA GUESP2(86,2)/     5.37194220D0/
      DATA GUESP3(86,2)/     1.48927868D0/
C     *
C *** START OF PDDG/MNDO PARAMETERS.
C     *
C     LIST OF ELEMENTS WITH PDDG/MNDO PARAMETERS.
      DATA IMPDM/1,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,1,0,
     1           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,
     2           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,
     3           25*0,                 0,0,0,0,0,0,0/
C     DATA FOR ELEMENT  1        HYDROGEN
C     DATA REFMN  ( 1)/'  H: (PDDG/MNDO):  M.P. REPASKY ET AL, J. COMP.
C     CHEM. 23, 1601 (2002)           '/
      DATA USSPDM ( 1)/     -11.7241140D0/
      DATA UPPPDM ( 1)/       0.0000000D0/
      DATA BSPDM  ( 1)/      -7.4935040D0/
      DATA BPPDM  ( 1)/      -7.4935040D0/
      DATA ZSPDM  ( 1)/       1.3224310D0/
      DATA ZPPDM  ( 1)/       1.3224310D0/
      DATA ALPPDM ( 1)/       2.4918130D0/
      DATA EISPDM ( 1)/     -12.0159190D0/
      DATA GSSPDM ( 1)/      12.8480000D0/
      DATA GSPPDM ( 1)/       0.0000000D0/
      DATA GPPPDM ( 1)/       0.0000000D0/
      DATA GP2PDM ( 1)/       0.0000000D0/
      DATA HSPPDM ( 1)/       0.0000000D0/
      DATA DDPDM  ( 1)/       0.0000000D0/
      DATA QQPDM  ( 1)/       0.0000000D0/
      DATA AMPDM  ( 1)/       0.4721793D0/
      DATA ADPDM  ( 1)/       0.4721793D0/
      DATA AQPDM  ( 1)/       0.4721793D0/
      DATA PAPDM  ( 1,1)/    -0.1088610D0/
      DATA PAPDM  ( 1,2)/    -0.0247060D0/
      DATA DAPDM  ( 1,1)/     0.4607210D0/
      DATA DAPDM  ( 1,2)/     1.2987310D0/
C     DATA FOR ELEMENT  6        CARBON
C     DATA REFMN  ( 6)/'  C: (PDDG/MNDO):  M.P. REPASKY ET AL, J. COMP.
C     CHEM. 23, 1601 (2002)           '/
      DATA USSPDM ( 6)/     -53.8375820D0/
      DATA UPPPDM ( 6)/     -39.9364090D0/
      DATA BSPDM  ( 6)/     -18.8413340D0/
      DATA BPPDM  ( 6)/      -7.9222340D0/
      DATA ZSPDM  ( 6)/       1.8098170D0/
      DATA ZPPDM  ( 6)/       1.8250080D0/
      DATA ALPPDM ( 6)/       2.5555220D0/
      DATA EISPDM ( 6)/    -123.8643930D0/
      DATA GSSPDM ( 6)/      12.2300000D0/
      DATA GSPPDM ( 6)/      11.4700000D0/
      DATA GPPPDM ( 6)/      11.0800000D0/
      DATA GP2PDM ( 6)/       9.8400000D0/
      DATA HSPPDM ( 6)/       2.4300000D0/
      DATA DDPDM  ( 6)/       0.7941580D0/
      DATA QQPDM  ( 6)/       0.6710900D0/
      DATA AMPDM  ( 6)/       0.4494671D0/
      DATA ADPDM  ( 6)/       0.6205810D0/
      DATA AQPDM  ( 6)/       0.6781010D0/
      DATA PAPDM  ( 6,1)/    -0.0068890D0/
      DATA PAPDM  ( 6,2)/    -0.0277510D0/
      DATA DAPDM  ( 6,1)/     1.1924560D0/
      DATA DAPDM  ( 6,2)/     1.3295220D0/
C     DATA FOR ELEMENT  7        NITROGEN
C     DATA REFMN  ( 7)/'  N: (PDDG/MNDO):  M.P. REPASKY ET AL, J. COMP.
C     CHEM. 23, 1601 (2002)           '/
      DATA USSPDM ( 7)/     -71.8718940D0/
      DATA UPPPDM ( 7)/     -58.2166170D0/
      DATA BSPDM  ( 7)/     -20.3757740D0/
      DATA BPPDM  ( 7)/     -21.0853730D0/
      DATA ZSPDM  ( 7)/       2.2314240D0/
      DATA ZPPDM  ( 7)/       2.2534600D0/
      DATA ALPPDM ( 7)/       2.8436780D0/
      DATA EISPDM ( 7)/    -206.4667560D0/
      DATA GSSPDM ( 7)/      13.5900000D0/
      DATA GSPPDM ( 7)/      12.6600000D0/
      DATA GPPPDM ( 7)/      12.9800000D0/
      DATA GP2PDM ( 7)/      11.5900000D0/
      DATA HSPPDM ( 7)/       3.1400000D0/
      DATA DDPDM  ( 7)/       0.6436240D0/
      DATA QQPDM  ( 7)/       0.5434950D0/
      DATA AMPDM  ( 7)/       0.4994487D0/
      DATA ADPDM  ( 7)/       0.7818860D0/
      DATA AQPDM  ( 7)/       0.8121110D0/
      DATA PAPDM  ( 7,1)/     0.0350270D0/
      DATA PAPDM  ( 7,2)/    -0.0017210D0/
      DATA DAPDM  ( 7,1)/     1.0116300D0/
      DATA DAPDM  ( 7,2)/     2.2784230D0/
C     DATA FOR ELEMENT  8        OXYGEN
C     DATA REFMN  ( 8)/'  O: (PDDG/MNDO):  M.P. REPASKY ET AL, J. COMP.
C     CHEM. 23, 1601 (2002)           '/
      DATA USSPDM ( 8)/     -97.8849700D0/
      DATA UPPPDM ( 8)/     -77.3426740D0/
      DATA BSPDM  ( 8)/     -33.6063360D0/
      DATA BPPDM  ( 8)/     -27.9844420D0/
      DATA ZSPDM  ( 8)/       2.5691720D0/
      DATA ZPPDM  ( 8)/       2.6971520D0/
      DATA ALPPDM ( 8)/       3.2388420D0/
      DATA EISPDM ( 8)/    -310.8796550D0/
      DATA GSSPDM ( 8)/      15.4200000D0/
      DATA GSPPDM ( 8)/      14.4800000D0/
      DATA GPPPDM ( 8)/      14.5200000D0/
      DATA GP2PDM ( 8)/      12.9800000D0/
      DATA HSPPDM ( 8)/       3.9400000D0/
      DATA DDPDM  ( 8)/       0.5473440D0/
      DATA QQPDM  ( 8)/       0.4540880D0/
      DATA AMPDM  ( 8)/       0.5667034D0/
      DATA ADPDM  ( 8)/       0.9471000D0/
      DATA AQPDM  ( 8)/       0.9489240D0/
      DATA PAPDM  ( 8,1)/     0.0863440D0/
      DATA PAPDM  ( 8,2)/     0.0304030D0/
      DATA DAPDM  ( 8,1)/     0.7254080D0/
      DATA DAPDM  ( 8,2)/     0.7203740D0/
C     DATA FOR ELEMENT  9        FLUORINE
C     DATA REFMN  ( 9)/'  F: (PDDG/MNDO):  I. TUBERT-BROHMAN ET AL, J.
C     COMP. CHEM., IN PRESS           '/
      DATA USSPDM ( 9)/    -134.2203790D0/
      DATA UPPPDM ( 9)/    -107.1559610D0/
      DATA BSPDM  ( 9)/     -67.8276120D0/
      DATA BPPDM  ( 9)/     -40.9248180D0/
      DATA ZSPDM  ( 9)/       4.3285190D0/
      DATA ZPPDM  ( 9)/       2.9050420D0/
      DATA ALPPDM ( 9)/       3.3223820D0/
      DATA EISPDM ( 9)/    -488.7032430D0/
      DATA GSSPDM ( 9)/      16.9200000D0/
      DATA GSPPDM ( 9)/      17.2500000D0/
      DATA GPPPDM ( 9)/      16.7100000D0/
      DATA GP2PDM ( 9)/      14.9100000D0/
      DATA HSPPDM ( 9)/       4.8300000D0/
      DATA DDPDM  ( 9)/       0.3615560D0/
      DATA QQPDM  ( 9)/       0.4215930D0/
      DATA AMPDM  ( 9)/       0.621830220D0/
      DATA ADPDM  ( 9)/       1.303600806D0/
      DATA AQPDM  ( 9)/       1.048409249D0/
      DATA PAPDM  ( 9,1)/    -0.0115790D0/
      DATA PAPDM  ( 9,2)/    -0.0129430D0/
      DATA DAPDM  ( 9,1)/     0.8346060D0/
      DATA DAPDM  ( 9,2)/     1.8756030D0/
C     DATA FOR ELEMENT 17        CHLORINE
C     DATA REFMN  (17)/' Cl: (PDDG/MNDO):  I. TUBERT-BROHMAN ET AL, J.
C     COMP. CHEM., IN PRESS           '/
      DATA USSPDM (17)/    -111.1336530D0/
      DATA UPPPDM (17)/     -78.0624930D0/
      DATA BSPDM  (17)/     -15.6633170D0/
      DATA BPPDM  (17)/     -15.3993310D0/
      DATA ZSPDM  (17)/       4.2124040D0/
      DATA ZPPDM  (17)/       2.0376470D0/
      DATA ALPPDM (17)/       2.6028460D0/
      DATA EISPDM (17)/    -378.9097270D0/
      DATA GSSPDM (17)/      15.0300000D0/
      DATA GSPPDM (17)/      13.1600000D0/
      DATA GPPPDM (17)/      11.3000000D0/
      DATA GP2PDM (17)/       9.9700000D0/
      DATA HSPPDM (17)/       2.4200000D0/
      DATA DDPDM  (17)/       0.4116090D0/
      DATA QQPDM  (17)/       0.8212020D0/
      DATA AMPDM  (17)/       0.552370221D0/
      DATA ADPDM  (17)/       0.902123237D0/
      DATA AQPDM  (17)/       0.605617221D0/
      DATA PAPDM  (17,1)/    -0.0171190D0/
      DATA PAPDM  (17,2)/     0.0054970D0/
      DATA DAPDM  (17,1)/     1.4663350D0/
      DATA DAPDM  (17,2)/     2.2368420D0/
C     DATA FOR ELEMENT 35        BROMINE
C     DATA REFMN  (35)/' Br: (PDDG/MNDO):  I. TUBERT-BROHMAN ET AL, J.
C     COMP. CHEM., IN PRESS           '/
      DATA USSPDM (35)/    -100.6370070D0/
      DATA UPPPDM (35)/     -76.0157350D0/
      DATA BSPDM  (35)/      -7.0541700D0/
      DATA BPPDM  (35)/     -10.2210300D0/
      DATA ZSPDM  (35)/       3.9999750D0/
      DATA ZPPDM  (35)/       2.2450400D0/
      DATA ALPPDM (35)/       2.4142650D0/
      DATA EISPDM (35)/    -349.5640960D0/
      DATA GSSPDM (35)/      15.03643948D0/
      DATA GSPPDM (35)/      13.03468242D0/
      DATA GPPPDM (35)/      11.27632539D0/
      DATA GP2PDM (35)/       9.85442552D0/
      DATA HSPPDM (35)/       2.45586832D0/
      DATA DDPDM  (35)/       0.5746230D0/
      DATA QQPDM  (35)/       0.9448920D0/
      DATA AMPDM  (35)/       0.552607090D0/
      DATA ADPDM  (35)/       0.747532768D0/
      DATA AQPDM  (35)/       0.564985796D0/
      DATA PAPDM  (35,1)/    -0.0171330D0/
      DATA PAPDM  (35,2)/    -0.0169640D0/
      DATA DAPDM  (35,1)/     2.2015390D0/
      DATA DAPDM  (35,2)/     2.2557640D0/
C     DATA FOR ELEMENT 53        IODINE
C     DATA REFMN  (53)/'  I: (PDDG/MNDO):  I. TUBERT-BROHMAN ET AL, J.
C     COMP. CHEM., IN PRESS           '/
      DATA USSPDM (53)/    -106.5884220D0/
      DATA UPPPDM (53)/     -75.2826050D0/
      DATA BSPDM  (53)/      -6.6983750D0/
      DATA BPPDM  (53)/      -5.6938140D0/
      DATA ZSPDM  (53)/       2.7184040D0/
      DATA ZPPDM  (53)/       2.4618130D0/
      DATA ALPPDM (53)/       2.2424460D0/
      DATA EISPDM (53)/    -356.0763980D0/
      DATA GSSPDM (53)/      15.04044855D0/
      DATA GSPPDM (53)/      13.05655798D0/
      DATA GPPPDM (53)/      11.14778369D0/
      DATA GP2PDM (53)/       9.91409071D0/
      DATA HSPPDM (53)/       2.45638202D0/
      DATA DDPDM  (53)/       1.2095290D0/
      DATA QQPDM  (53)/       1.0435590D0/
      DATA AMPDM  (53)/       0.552753708D0/
      DATA ADPDM  (53)/       0.498848160D0/
      DATA AQPDM  (53)/       0.504031750D0/
      DATA PAPDM  (53,1)/     0.0096160D0/
      DATA PAPDM  (53,2)/    -0.0075050D0/
      DATA DAPDM  (53,1)/     2.5723320D0/
      DATA DAPDM  (53,2)/     2.9364560D0/
C     *
C     START OF PDDG/PM3 PARAMETER SET.
C     *
C     NUMBER OF GAUSSIANS USED PER ATOM.
      DATA IMPDP/2,0,0,0,0,2,2,2,2,0,0,0,0,0,0,0,2,0,
     1           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,
     2           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,
     3           25*0,                 0,0,0,0,0,0,0/
C     DATA FOR ELEMENT  1        HYDROGEN
C     DATA REFPDP ( 1)/'  H: (PDDG/PM3):  M.P. REPASKY ET AL, J. COMP.
C     CHEM. 23, 1601 (2002)           '/
      DATA USSPDP(  1)/     -12.8932720D0/
      DATA UPPPDP(  1)/       0.0000000D0/
      DATA BSPDP (  1)/      -6.1526540D0/
      DATA BPPDP (  1)/      -6.1526540D0/
      DATA ZSPDP (  1)/       0.9727860D0/
      DATA ZPPDP (  1)/       0.9727860D0/
      DATA ALPPDP(  1)/       3.3816860D0/
      DATA EISPDP(  1)/     -13.1205660D0/
      DATA GSSPDP(  1)/      14.7942080D0/
      DATA GSPPDP(  1)/       0.0000000D0/
      DATA GPPPDP(  1)/       0.0000000D0/
      DATA GP2PDP(  1)/       0.0000000D0/
      DATA HSPPDP(  1)/       0.0000000D0/
      DATA DDPDP (  1)/       0.0000000D0/
      DATA QQPDP (  1)/       0.0000000D0/
      DATA AMPDP (  1)/       0.5437048D0/
      DATA ADPDP (  1)/       0.5437048D0/
      DATA AQPDP (  1)/       0.5437048D0/
      DATA PAPDP (  1,1)/     0.0571930D0/
      DATA PAPDP (  1,2)/    -0.0348230D0/
      DATA DAPDP (  1,1)/     0.6633950D0/
      DATA DAPDP (  1,2)/     1.0819010D0/
      DATA GUPDP1(  1,1)/       1.1222440D0/
      DATA GUPDP2(  1,1)/       4.7077900D0/
      DATA GUPDP3(  1,1)/       1.5470990D0/
      DATA GUPDP1(  1,2)/      -1.0697370D0/
      DATA GUPDP2(  1,2)/       5.8579950D0/
      DATA GUPDP3(  1,2)/       1.5678930D0/
C     DATA FOR ELEMENT  6      CARBON
C     DATA REFPDP ( 6)/'  C: (PDDG/PM3):  M.P. REPASKY ET AL, J. COMP.
C     CHEM. 23, 1601 (2002)           '/
      DATA USSPDP(  6)/     -48.2412410D0/
      DATA UPPPDP(  6)/     -36.4612560D0/
      DATA BSPDP (  6)/     -11.9528180D0/
      DATA BPPDP (  6)/      -9.9224110D0/
      DATA ZSPDP (  6)/       1.5678640D0/
      DATA ZPPDP (  6)/       1.8466590D0/
      DATA ALPPDP(  6)/       2.7257720D0/
      DATA EISPDP(  6)/    -113.4282420D0/
      DATA GSSPDP(  6)/      11.2007080D0/
      DATA GSPPDP(  6)/      10.2650270D0/
      DATA GPPPDP(  6)/      10.7962920D0/
      DATA GP2PDP(  6)/       9.0425660D0/
      DATA HSPPDP(  6)/       2.2909800D0/
      DATA DDPDP (  6)/       0.8314130D0/
      DATA QQPDP (  6)/       0.6632220D0/
      DATA AMPDP (  6)/       0.4116394D0/
      DATA ADPDP (  6)/       0.5892980D0/
      DATA AQPDP (  6)/       0.7659490D0/
      DATA PAPDP (  6,1)/    -0.0007430D0/
      DATA PAPDP (  6,2)/     0.0009850D0/
      DATA DAPDP (  6,1)/     0.8369150D0/
      DATA DAPDP (  6,2)/     1.5852360D0/
      DATA GUPDP1(  6,1)/       0.0489060D0/
      DATA GUPDP2(  6,1)/       5.7653400D0/
      DATA GUPDP3(  6,1)/       1.6822320D0/
      DATA GUPDP1(  6,2)/       0.0476970D0/
      DATA GUPDP2(  6,2)/       5.9737210D0/
      DATA GUPDP3(  6,2)/       0.8944060D0/
C     DATA FOR ELEMENT  7      NITROGEN
C     DATA REFPDP ( 7)/'  N: (PDDG/PM3):  M.P. REPASKY ET AL, J. COMP.
C     CHEM. 23, 1601 (2002)           '/
      DATA USSPDP(  7)/     -49.4545460D0/
      DATA UPPPDP(  7)/     -47.7574060D0/
      DATA BSPDP (  7)/     -14.1172300D0/
      DATA BPPDP (  7)/     -19.9385090D0/
      DATA ZSPDP (  7)/       2.0358070D0/
      DATA ZPPDP (  7)/       2.3243270D0/
      DATA ALPPDP(  7)/       2.8491240D0/
      DATA EISPDP(  7)/    -158.4162050D0/
      DATA GSSPDP(  7)/      11.9047870D0/
      DATA GSPPDP(  7)/       7.3485650D0/
      DATA GPPPDP(  7)/      11.7546720D0/
      DATA GP2PDP(  7)/      10.8072770D0/
      DATA HSPPDP(  7)/       1.1367130D0/
      DATA DDPDP (  7)/       0.6548550D0/
      DATA QQPDP (  7)/       0.5269240D0/
      DATA AMPDP (  7)/       0.4375151D0/
      DATA ADPDP (  7)/       0.5044210D0/
      DATA AQPDP (  7)/       0.7388760D0/
      DATA PAPDP (  7,1)/    -0.0031600D0/
      DATA PAPDP (  7,2)/     0.0125010D0/
      DATA DAPDP (  7,1)/     1.0041720D0/
      DATA DAPDP (  7,2)/     1.5163360D0/
      DATA GUPDP1(  7,1)/       1.5133200D0/
      DATA GUPDP2(  7,1)/       5.9043940D0/
      DATA GUPDP3(  7,1)/       1.7283760D0/
      DATA GUPDP1(  7,2)/      -1.5118920D0/
      DATA GUPDP2(  7,2)/       6.0300140D0/
      DATA GUPDP3(  7,2)/       1.7341080D0/
C     DATA FOR ELEMENT  8      OXYGEN
C     DATA REFPDP ( 8)/'  O: (PDDG/PM3):  M.P. REPASKY ET AL, J. COMP.
C     CHEM. 23, 1601 (2002)           '/
      DATA USSPDP(  8)/     -87.4125050D0/
      DATA UPPPDP(  8)/     -72.1830700D0/
      DATA BSPDP (  8)/     -44.8745530D0/
      DATA BPPDP (  8)/     -24.6019390D0/
      DATA ZSPDP (  8)/       3.8145650D0/
      DATA ZPPDP (  8)/       2.3180110D0/
      DATA ALPPDP(  8)/       3.2253090D0/
      DATA EISPDP(  8)/    -292.1887660D0/
      DATA GSSPDP(  8)/      15.7557600D0/
      DATA GSPPDP(  8)/      10.6211600D0/
      DATA GPPPDP(  8)/      13.6540160D0/
      DATA GP2PDP(  8)/      12.4060950D0/
      DATA HSPPDP(  8)/       0.5938830D0/
      DATA DDPDP (  8)/       0.4037410D0/
      DATA QQPDP (  8)/       0.5283600D0/
      DATA AMPDP (  8)/       0.5790430D0/
      DATA ADPDP (  8)/       0.5340360D0/
      DATA AQPDP (  8)/       0.8009090D0/
      DATA PAPDP (  8,1)/    -0.0010000D0/
      DATA PAPDP (  8,2)/    -0.0015220D0/
      DATA DAPDP (  8,1)/     1.3606850D0/
      DATA DAPDP (  8,2)/     1.3664070D0/
      DATA GUPDP1(  8,1)/      -1.1384550D0/
      DATA GUPDP2(  8,1)/       6.0000430D0/
      DATA GUPDP3(  8,1)/       1.6223620D0/
      DATA GUPDP1(  8,2)/       1.1460070D0/
      DATA GUPDP2(  8,2)/       5.9634940D0/
      DATA GUPDP3(  8,2)/       1.6147880D0/
C     DATA FOR ELEMENT  9        FLUORINE
C     DATA REFMN  ( 9)/'  F: (PDDG/PM3):  I. TUBERT-BROHMAN ET AL, J.
C     COMP. CHEM., IN PRESS           '/
      DATA USSPDP(  9)/    -111.4004320D0/
      DATA UPPPDP(  9)/    -106.3952640D0/
      DATA BSPDP (  9)/     -50.9373010D0/
      DATA BPPDP (  9)/     -31.6369760D0/
      DATA ZSPDP (  9)/       5.5380330D0/
      DATA ZPPDP (  9)/       2.5380660D0/
      DATA ALPPDP(  9)/       3.2005710D0/
      DATA EISPDP(  9)/    -442.4571330D0/
      DATA GSSPDP(  9)/      10.4966670D0/
      DATA GSPPDP(  9)/      16.0736890D0/
      DATA GPPPDP(  9)/      14.8172560D0/
      DATA GP2PDP(  9)/      14.4183930D0/
      DATA HSPPDP(  9)/       0.7277630D0/
      DATA DDPDP (  9)/       0.2466010D0/
      DATA QQPDP (  9)/       0.4825510D0/
      DATA AMPDP (  9)/       0.385764964D0/
      DATA ADPDP (  9)/       0.787856919D0/
      DATA AQPDP (  9)/       0.620499825D0/
      DATA PAPDP (  9,1)/    -0.0128660D0/
      DATA PAPDP (  9,2)/     0.0073150D0/
      DATA DAPDP (  9,1)/     1.3056810D0/
      DATA DAPDP (  9,2)/     1.8425720D0/
      DATA GUPDP1(  9,1)/      -0.0080790D0/
      DATA GUPDP2(  9,1)/       5.9389690D0/
      DATA GUPDP3(  9,1)/       1.8639490D0/
      DATA GUPDP1(  9,2)/      -0.0026590D0/
      DATA GUPDP2(  9,2)/       5.9251050D0/
      DATA GUPDP3(  9,2)/       2.3888640D0/
C     DATA FOR ELEMENT 17        CHLORINE
C     DATA REFMN  (17)/' Cl: (PDDG/PM3):  I. TUBERT-BROHMAN ET AL, J.
C     COMP. CHEM., IN PRESS           '/
      DATA USSPDP( 17)/     -95.0944340D0/
      DATA UPPPDP( 17)/     -53.9216510D0/
      DATA BSPDP ( 17)/     -26.9131290D0/
      DATA BPPDP ( 17)/     -14.9911780D0/
      DATA ZSPDP ( 17)/       2.5482680D0/
      DATA ZPPDP ( 17)/       2.2846240D0/
      DATA ALPPDP( 17)/       2.4976170D0/
      DATA EISPDP( 17)/    -305.7152010D0/
      DATA GSSPDP( 17)/      16.0136010D0/
      DATA GSPPDP( 17)/       8.0481150D0/
      DATA GPPPDP( 17)/       7.5222150D0/
      DATA GP2PDP( 17)/       7.5041540D0/
      DATA HSPPDP( 17)/       3.4811530D0/
      DATA DDPDP ( 17)/       0.8275610D0/
      DATA QQPDP ( 17)/       0.7324270D0/
      DATA AMPDP ( 17)/       0.588519168D0/
      DATA ADPDP ( 17)/       0.718221568D0/
      DATA AQPDP ( 17)/       0.217476025D0/
      DATA PAPDP ( 17,1)/    -0.0165520D0/
      DATA PAPDP ( 17,2)/    -0.0166460D0/
      DATA DAPDP ( 17,1)/     1.7276900D0/
      DATA DAPDP ( 17,2)/     1.7846550D0/
      DATA GUPDP1( 17,1)/      -0.1122220D0/
      DATA GUPDP2( 17,1)/       5.9637190D0/
      DATA GUPDP3( 17,1)/       1.0277190D0/
      DATA GUPDP1( 17,2)/      -0.0130610D0/
      DATA GUPDP2( 17,2)/       1.9995560D0/
      DATA GUPDP3( 17,2)/       2.2863770D0/
C     DATA FOR ELEMENT 35        BROMINE
C     DATA REFMN  (35)/' Br: (PDDG/PM3):  I. TUBERT-BROHMAN ET AL, J.
C     COMP. CHEM., IN PRESS           '/
      DATA USSPDP( 35)/    -115.8419630D0/
      DATA UPPPDP( 35)/     -74.2051460D0/
      DATA BSPDP ( 35)/     -21.5380440D0/
      DATA BPPDP ( 35)/      -8.5247640D0/
      DATA ZSPDP ( 35)/       4.3450790D0/
      DATA ZPPDP ( 35)/       2.1909610D0/
      DATA ALPPDP( 35)/       2.4246730D0/
      DATA EISPDP( 35)/    -351.0138870D0/
      DATA GSSPDP( 35)/      15.9434250D0/
      DATA GSPPDP( 35)/      16.0616800D0/
      DATA GPPPDP( 35)/       8.2827630D0/
      DATA GP2PDP( 35)/       7.8168490D0/
      DATA HSPPDP( 35)/       0.5788690D0/
      DATA DDPDP ( 35)/       0.4738600D0/
      DATA QQPDP ( 35)/       0.9682140D0/
      DATA AMPDP ( 35)/       0.585939789D0/
      DATA ADPDP ( 35)/       0.477815047D0/
      DATA AQPDP ( 35)/       0.390428870D0/
      DATA PAPDP ( 35,1)/    -0.0137720D0/
      DATA PAPDP ( 35,2)/     0.0088490D0/
      DATA DAPDP ( 35,1)/     1.8520300D0/
      DATA DAPDP ( 35,2)/     2.3389580D0/
      DATA GUPDP1( 35,1)/       0.9613620D0/
      DATA GUPDP2( 35,1)/       6.0136000D0/
      DATA GUPDP3( 35,1)/       2.3404450D0/
      DATA GUPDP1( 35,2)/      -0.9488340D0/
      DATA GUPDP2( 35,2)/       5.9763290D0/
      DATA GUPDP3( 35,2)/       2.3487450D0/
C     DATA FOR ELEMENT 53        IODINE
C     DATA REFMN  (53)/'  I: (PDDG/PM3):  I. TUBERT-BROHMAN ET AL, J.
C     COMP. CHEM., IN PRESS           '/
      DATA USSPDP( 53)/     -97.6641740D0/
      DATA UPPPDP( 53)/     -61.1671370D0/
      DATA BSPDP ( 53)/     -16.5926210D0/
      DATA BPPDP ( 53)/      -6.5998160D0/
      DATA ZSPDP ( 53)/       5.0628010D0/
      DATA ZPPDP ( 53)/       2.4177570D0/
      DATA ALPPDP( 53)/       1.9781700D0/
      DATA EISPDP( 53)/    -291.5378690D0/
      DATA GSSPDP( 53)/      13.6319430D0/
      DATA GSPPDP( 53)/      14.9904060D0/
      DATA GPPPDP( 53)/       7.2883300D0/
      DATA GP2PDP( 53)/       5.9664070D0/
      DATA HSPPDP( 53)/       2.6300350D0/
      DATA DDPDP ( 53)/       0.4072610D0/
      DATA QQPDP ( 53)/       1.0625740D0/
      DATA AMPDP ( 53)/       0.500989956D0/
      DATA ADPDP ( 53)/       0.939337579D0/
      DATA AQPDP ( 53)/       0.510317080D0/
      DATA PAPDP ( 53,1)/     0.0129010D0/
      DATA PAPDP ( 53,2)/    -0.0128250D0/
      DATA DAPDP ( 53,1)/     1.9942990D0/
      DATA DAPDP ( 53,2)/     2.2634170D0/
      DATA GUPDP1( 53,1)/      -0.1360030D0/
      DATA GUPDP2( 53,1)/       3.8529120D0/
      DATA GUPDP3( 53,1)/       1.6974550D0/
      DATA GUPDP1( 53,2)/      -0.0372870D0/
      DATA GUPDP2( 53,2)/       5.2292640D0/
      DATA GUPDP3( 53,2)/       2.7686690D0/
      END
