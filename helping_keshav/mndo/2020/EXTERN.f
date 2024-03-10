      SUBROUTINE EXTERN (IPAROK,IOPPRT)
C     *
C     INPUT OF EXTERNAL PARAMETERS FROM FILE NB14.
C     *
C     NOTATION. I=INPUT, O=OUTPUT.
C     IPAROK    CONVENTIONS FOR CONTENTS OF FILE NB14 (I).
C               = 0 DEFAULT (HANDLED AS IPAROK=1).
C               = 1 KEYWORD-ORIENTED INPUT ACCORDING TO MOPAC(6.0).
C               = 2 NUMERICAL INPUT, IPAROK=2, INTERNAL NUMBERING.
C               = 3 NUMERICAL INPUT, IPAROK=3, STANDARD NUMBERING.
C     IOPPRT    PRINTING OPTION (I).
C               =-5 NO PRINTING.
C     *
C     IN MOPAC(6.0) THE KEYWORD EXTERNAL=filename PROVIDES THE OPTION
C     TO READ ADDITIONAL PARAMETERS FROM AN EXTERNAL FILE filename.
C     WHEN USING MOPAC-TYPE INPUT, THE CONVENTIONS FROM MOPAC(6.0) ARE
C     FOLLOWED. THE EXTERNAL FILE CONTAINING THE PARAMETERS SHOULD BE
C     IDENTICAL FOR MOPAC(6.0) AND THIS PROGRAM, ALTHOUGH THE TECHNICAL
C     IMPLEMENTATION IS DIFFERENT (NUMBERING SCHEMES ETC, SEE CODE).
C     *
C     WHEN USING STANDARD INPUT, THE TYPE OF PARAMETER AND THE ELEMENT
C     ARE REPRESENTED NUMERICALLY. THE NUMBERING SCHEME FOR THE TYPES
C     OF PARAMETERS ARE GIVEN BELOW AS COMMENTS (TWO SCHEMES AVAILABLE,
C     IPAROK=2 AND 3).
C     *
      USE LIMIT, ONLY: LMZ, LMP
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (LMPAR=154)
      CHARACTER PARTYP(LMPAR)*5
      CHARACTER GETNAM*80,FILES*80
      CHARACTER TEXT*80,DUMMY*80
      CHARACTER KEYWRD*800
      CHARACTER*2 ELEMNT
      CHARACTER*2 ELCOPY
      COMMON
     ./ELEMTS/ ELEMNT(107)
     ./MOPAC / IMOPAC
     ./KEYWRD/ KEYWRD
     ./NBFILE/ NBF(20)
     ./PAR20 / X(LMP),NF(LMP),NS(LMP),NNPAR
C     DIMENSION NEW(LMZ)
      DIMENSION MOPPAR(LMPAR)
      SAVE PARTYP, MOPPAR
C *** INTERNAL NUMBERING SCHEME FOR THE PARAMETERS (VARIABLE IPARAM).
C     THIS NUMBERING SCHEME DIFFERS FROM MOPAC(6.0), BUT THE NAMES OF
C     THE PARAMETERS ARE THE SAME AS IN MOPAC(6.0) SO THAT THE INPUT
C     FILES ARE IDENTICAL FOR MOPAC(6.0) AND THIS PROGRAM WHEN USING
C     MOPAC-TYPE INPUT.
C      1 USS  ,   2 UPP  ,   3 ZS   ,   4 ZP   ,   5 BETAS,   6 BETAP,
C      7 ALPHA,   8 BETPI,   9 BETSH,  10 BETPH,  11 ALPS ,  12 ALPP ,
C     13 ALPPI,  14 ALPSH,  15 ALPPH,  16 FVAL1,  17 FVAL2,  18 GVAL1,
C     19 GVAL2,  20 FG   ,  21 UDD  ,  22 ZD   ,  23 BETAD,  24 ZSN  ,
C     25 ZPN  ,  26 ZDN  ,  27 POCOR,  28 GSCAL,  29 NUMAO,  30 XX30 ,
C     31 FN11 ,  32 FN21 ,  33 FN31 ,  34 FN12 ,  35 FN22 ,  36 FN32 ,
C     37 FN13 ,  38 FN23 ,  39 FN33 ,  40 FN14 ,  41 FN24 ,  42 FN34 ,
C     43 SOLV1,  44 SOLV2,  45 SOLV3,  46 SOLV4,  47 SOLV5,  48 SOLV6,
C     49 ZSCOR,  50 FSCOR,  51 BSCOR,  52 ASCOR,  53 DELTA,  54 OMEGA,
C     55 XSCAL,  56 XOFFL,  57 XOFFG,  58 ZSSCF,  59 ZPSCF,  60 ZDSCF,
C     61 BSSCF,  62 BPSCF,  63 BDSCF,  64 XUSS ,  65 XUPP ,  66 XUDD ,
C     67 ZSNMR,  68 ZPNMR,  69 ZDNMR,  70 BSNMR,  71 BPNMR,  72 BDNMR,
C     73 GSS  ,  74 GPP  ,  75 GSP  ,  76 GP2  ,  77 HSP  ,  78 HPP  ,
C     79 EHEAT,  80 F0DD ,  81 F2DD ,  82 F4DD ,  83 F0SD ,  84 G2SD ,
C     85 ALP01,  86 ALP06,  87 ALP07,  88 ALP08,  89 ALP09,  90 ALP14,
C     91 ALP15,  92 ALP16,  93 ALP17,  94 ALP35,  95 ALP53,  96 PDDG1,
C     97 PDDG2,  98 PDDG3,  99 PDDG4, 100 GNN  , 101 CPEZ , 102 CPEED,
C    103 CPEQ0, 104 MXCAT, 105 MXAN , 106 CPEFT, 107 ATPT , 108 ANPT ,
C    109 ATPS , 110 ANPS , 111 CPESL, 112 CPEEM, 113 CPEM0, 114 PQNEF,
C    115 PQDP , 116 PQQP , 117 PQOP , 118 PQB  , 119 CPERS, 120 CPERR,
C    121 CPERC, 122 CPEQM, 123 PQBO0, 124 PQBOZ, 125 CPEPS, 126 CPEDA,
C    127 CPEZB, 128 PQRL , 129 PQRH , 130 PQDS , 131 PQQS , 132 PQOS ,
C    133 SCRQ0, 134 SCRQ1, 135 SCRQ2, 136 SCNET, 137 SCNEB, 138 PBETS,
C    139 PBETP, 140 PBETD, 141 PALP , 142 SDDP , 143 SDQP , 144 SDDPQ,
C    145 SDQPQ, 146 SDRLD, 147 SDRLQ, 148 SDRHD, 149 SDRHQ, 150 SDNEF,
C    151 D3S6 , 151 D3S8 , 151 D3A1 , 151 D3A2
      DATA PARTYP/'USS  ','UPP  ','ZS   ','ZP   ','BETAS','BETAP',
     1            'ALP  ','BETPI','BETSH','BETPH','ALPS ','ALPP ',
     2            'ALPPI','ALPSH','ALPPH','FVAL1','FVAL2','GVAL1',
     3            'GVAL2','FG   ','UDD  ','ZD   ','BETAD','ZSN  ',
     4            'ZPN  ','ZDN  ','POCOR','GSCAL','NUMAO','XX30 ',
     5            'FN11 ','FN21 ','FN31 ','FN12 ','FN22 ','FN32 ',
     6            'FN13 ','FN23 ','FN33 ','FN14 ','FN24 ','FN34 ',
     7            'SOLV1','SOLV2','SOLV3','SOLV4','SOLV5','SOLV6',
     8            'ZSCOR','FSCOR','BSCOR','ASCOR','DELTA','OMEGA',
     9            'XSCAL','XOFFL','XOFFG','ZSSCF','ZPSCF','ZDSCF',
     A            'BSSCF','BPSCF','BDSCF','XUSS ','XUPP ','XUDD ',
     B            'ZSNMR','ZPNMR','ZDNMR','BSNMR','BPNMR','BDNMR',
     C            'GSS  ','GPP  ','GSP  ','GP2  ','HSP  ','HPP  ',
     D            'EHEAT','F0DD ','F2DD ','F4DD ','F0SD ','G2SD ',
     E            'ALP01','ALP06','ALP07','ALP08','ALP09','ALP14',
     F            'ALP15','ALP16','ALP17','ALP35','ALP53','PDDG1',
     G            'PDDG2','PDDG3','PDDG4','GNN  ','CPEZ ','CPEED',
     H            'CPEQ0','MXCAT','MXAN ','CPEFT','ATPT ','ANPT ',
     I            'ATPS ','ANPS ','CPESL','CPEEM','CPEM0','PQNEF',
     J            'PQDP ','PQQP ','PQOP ','PQB  ','CPERS','CPERR',
     K            'CPERC','CPEQM','PQBO0','PQBOZ','CPEPS','CPEDA',
     L            'CPEZB','PQRL ','PQRH ','PQDS ','PQQS ','PQOS ',
     M            'SCRQ0','SCRQ1','SCRQ2','SCNET','SCNEB','PBETS',
     N            'PBETP','PBETD','PALP ','SDDP ','SDQP ','SDDPQ',
     O            'SDQPQ','SDRLD','SDRLQ','SDRHD','SDRHQ','SDNEF',
     P            'D3S6' ,'D3S8' ,'D3A1' ,'D3A2' /

C *** MAPPING BETWEEN INTERNAL NUMBERING AND STANDARD NUMBERING.
C     PARAMETER I FOR INTERNAL NUMBERING CORRESPONDS TO PARAMETER
C     MOPPAR(I) IN STANDARD NUMBERING (PP97 CONVENTIONS).
C     THE NUMBERING SCHEMES ARE IDENTICAL UP TO I=72.
      DATA MOPPAR/   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,
     1              11,  12,  13,  14,  15,  16,  17,  18,  19,  20,
     2              21,  22,  23,  24,  25,  26,  27,  28,  29,  30,
     3              31,  32,  33,  34,  35,  36,  37,  38,  39,  40,
     4              41,  42,  43,  44,  45,  46,  47,  48,  49,  50,
     5              51,  52,  53,  54,  55,  56,  57,  58,  59,  60,
     6              61,  62,  63,  64,  65,  66,  67,  68,  69,  70,
     7              71,  72,  91,  92,  93,  94,  95,  96,  99, 121,
     8             122, 123, 124, 125, 201, 206, 207, 208, 209, 214,
     9             215, 216, 217, 235, 253,  73,  74,  75,  76, 998,
     A            1042,1043,1044,1054,1055,1056,1045,1046,1047,1048,
     B            1049,1050,1051,1060,1061,1062,1063,1064,1052,1053,
     C            1067,1068,1065,1066,1069,1070,1071,1072,1073,1074,
     D            1075,1076,1080,1081,1082,1083,1084,1090,1091,1092,
     E            1093,2000,2001,2002,2003,2004,2005,2006,2007,2008,
     F              79,  80,  81,  82/
C *** FILE NUMBERS.
      NB6    = NBF(6)
      NB14   = NBF(14)
C *** INITIALIZATION.
      KSTOP  = 0
      KSTOP1 = 0
      NNPAR  = 0
C     DO 10 I=1,LMZ
C  10 NEW(I) = 0
      IF(IOPPRT.GT.-5) WRITE(NB6,500)
      IF(IPAROK.GT.1) GO TO 100
C
C *** MOPAC-TYPE INPUT.
C     SOME OF THE CODE BELOW IS ADAPTED FROM MOPAC(6.0).
      IF(IMOPAC.NE.0) THEN
         I   = INDEX(KEYWRD,'EXTERNAL=')+9
         IF(I.GT.9) THEN
            J      = INDEX(KEYWRD(I:),' ')+I-1
            FILES  = GETNAM(KEYWRD(I:J))
            OPEN(NB14,STATUS='UNKNOWN',FILE=FILES)
         ENDIF
      ENDIF
      REWIND NB14
C *** INPUT LOOP.
   20 CONTINUE
      READ(NB14,510,ERR=300,END=200) TEXT
      IF(TEXT.EQ.' ') GO TO 200
      CALL UPPCAS(TEXT,80)
      IF(INDEX(TEXT,'END').NE.0) GO TO 200
C     TYPE OF PARAMETER: IPARAM=J.
      DO 30 J=1,LMPAR
      II = INDEX(TEXT,PARTYP(J))
      IF(II.GT.0) THEN
C        THE PARAMETER KEYWORD APPEARS ON THE INPUT LINE.
C        NOW MAKE SURE IT'S NOT JUST A SUBSTRING OF ANOTHER KEYWORD.
C        IF THE KEYWORD STARTS AT THE FIRST COLUMN OF THE INPUT LINE
C        OR THE PRECEEDING CHARACTER IS A BLANK, THE KEYWORD IS REAL.
         IF(II.EQ.1)                GO TO 40
         IF(TEXT(II-1:II-1).EQ.' ') GO TO 40
      ENDIF
   30 CONTINUE
      WRITE(NB6,520) TEXT
      KSTOP  = KSTOP+1
      GO TO 20
   40 IPARAM = J
C     ELEMENT: IELMNT=J.
      I      = II
      J      = INDEX(TEXT(I:),' ')+1
C     FIND THE NEXT CHARACTER THAT IS NOT BLANK.
      DO K=J,80
         IF(TEXT(K:K).NE.' ') GO TO 50
      ENDDO
      WRITE(NB6,530) TEXT
      KSTOP  = KSTOP+1
      GO TO 20
   50 DUMMY  = TEXT(K:)
      DO 60 J=1,LMZ
      ELCOPY = ELEMNT(J)
      CALL UPPCAS(ELCOPY,2)
      IF(INDEX(DUMMY(1:2),ELCOPY).NE.0) GO TO 70
   60 CONTINUE
      WRITE(NB6,530) TEXT
      KSTOP  = KSTOP+1
      GO TO 20
   70 IELMNT = J
C
C     IF WE ARE DEALING WITH A TWO-BODY PARAMETER, READ THE
C     SECOND ELEMENT AND LET IELMNT = 100*IELMNT1 + IELMNT2.
      IF(PARTYP(IPARAM)(1:4).EQ.'PQB '  .OR.
     1   PARTYP(IPARAM)(1:5).EQ.'PBETS' .OR.
     2   PARTYP(IPARAM)(1:5).EQ.'PBETP' .OR.
     3   PARTYP(IPARAM)(1:5).EQ.'PBETD' .OR.
     4   PARTYP(IPARAM)(1:4).EQ.'PALP') THEN
         J       = K+2
         DO K=J,80
            IF(TEXT(K:K).NE.' ') GOTO 80
         ENDDO
         WRITE(NB6,530) TEXT
         KSTOP  = KSTOP+1
         GO TO 20
   80    DUMMY  = TEXT(K:)
         DO 90 J=1,LMZ
         ELCOPY = ELEMNT(J)
         CALL UPPCAS(ELCOPY,2)
         IF(INDEX(DUMMY(1:2),ELCOPY).NE.0) GO TO 95
   90    CONTINUE
         WRITE(6,530) TEXT
         KSTOP  = KSTOP+1
         GO TO 20
   95    IELMNT = 100*IELMNT + J
      ENDIF
C     VALUE OF PARAMETER.
      PARAM = READA(DUMMY,INDEX(DUMMY,ELCOPY),0)
C     NEW ENTRY IN COMMON BLOCK.
      IF(NNPAR.LT.LMP) THEN
         NNPAR = NNPAR+1
         X (NNPAR) = PARAM
         NF(NNPAR) = MOPPAR(IPARAM)
         NS(NNPAR) = IELMNT
      ELSE
         KSTOP1 = KSTOP1+1
         GO TO 20
      ENDIF
C     NEW(IELMNT) = IELMNT
C     PRINTING SECTION.
      IF(IOPPRT.GT.-5) THEN
         IF(IELMNT.LT.100) THEN
            WRITE(NB6,540) PARTYP(IPARAM),ELEMNT(IELMNT),'  ',PARAM
         ELSE
            WRITE(NB6,540) PARTYP(IPARAM),ELEMNT(IELMNT/100),
     1                     ELEMNT(MOD(IELMNT,100)),PARAM
         ENDIF
      ENDIF
      GO TO 20
C     STANDARD INPUT WILL NOT WORK FOR TWO-BODY
C     PARAMETERS DUE TO LIMITATION OF FORMAT 600.
C     THUS PRINTING SECTION HAS NOT BEEN ADJUSTED.
C     *
C *** STANDARD INPUT.
C     OPTION IPAROK=2: INTERNAL NUMBERING FOR IPARAM (SEE ABOVE).
C     OPTION IPAROK=3: STANDARD NUMBERING SCHEME (PP97 CONVENTIONS)
C     WHICH IS DEFINED IN THE FOLLOWING LIST.
C       1 USS  ,   2 UPP  ,   3 ZS   ,   4 ZP   ,   5 BETAS,   6 BETAP,
C       7 ALPHA,   8 BETPI,   9 BETSH,  10 BETPH,  11 ALPS ,  12 ALPP ,
C      13 ALPPI,  14 ALPSH,  15 ALPPH,  16 FVAL1,  17 FVAL2,  18 GVAL1,
C      19 GVAL2,  20 FG   ,  21 UDD  ,  22 ZD   ,  23 BETAD,  24 ZSN  ,
C      25 ZPN  ,  26 ZDN  ,  27 POCOR,  28 GSCAL,  29 NUMAO,  30      ,
C      31 FN11 ,  32 FN21 ,  33 FN31 ,  34 FN12 ,  35 FN22 ,  36 FN32 ,
C      37 FN13 ,  38 FN23 ,  39 FN33 ,  40 FN14 ,  41 FN24 ,  42 FN34 ,
C      43 SOLV1,  44 SOLV2,  45 SOLV3,  46 SOLV4,  47 SOLV5,  48 SOLV6,
C      49 ZSCOR,  50 FSCOR,  51 BSCOR,  52 ASCOR,  53 DELTA,  54 OMEGA,
C      55 XSCAL,  56 XOFFL,  57 XOFFG,  58 ZSSCF,  59 ZPSCF,  60 ZDSCF,
C      61 BSSCF,  62 BPSCF,  63 BDSCF,  64 XUSS ,  65 XUPP ,  66 XUDD ,
C      67 ZSNMR,  68 ZPNMR,  69 ZDNMR,  70 BSNMR,  71 BPNMR,  72 BDNMR,
C      73 PDDG1,  74 PDDG2,  75 PDDG3,  76 PDDG4,  77 C6   ,  78 R0   ,
C      79 D3S6 ,  80 D3S8 ,  81 D3A1 ,  82 D3A2 ,  83      ,  84      ,
C      91 GSS  ,  92 GPP  ,  93 GSP  ,  94 GP2  ,  95 HSP  ,  96 HPP  ,
C      97      ,  98      ,  99 EHEAT, 100 CORE , 101      , 102      ,
C     121 F0DD , 122 F2DD , 123 F4DD , 124 F0SD , 125 G2SD , 126 F0PD ,
C     127 F2PD , 128 G1PD , 129 G3PD , 130      , 131      , 132      .
C     201 ALP01, 202 ALP02, 203 ALP03, 204 ALP04, 205 ALP05, 206 ALP06,
C     207 ALP07, 208 ALP08, 209 ALP09, 210 ALP10, etc up to  286 ALP86.
C     THE INPUT FOR IPARAM IS NUMERICAL, BOTH FOR IPAROK=2 AND 3.
C     THE INPUT FOR IELMNT IS ALSO NUMERICAL (ATOMIC NUMBER).
C     *
  100 CONTINUE
      REWIND NB14
C     INPUT LOOP.
  110 READ(NB14,600,ERR=300,END=200) IPARAM,IELMNT,PARAM
      IF(IPARAM.EQ.0) GO TO 200
C     CHECK FOR SIMPLE INPUT ERRORS.
      IF(IPARAM.LT.1 .OR. (IPARAM.GT.LMPAR .AND. IPAROK.EQ.2) .OR.
     1   IELMNT.LT.1 .OR. IELMNT.GT.LMZ) THEN
         KSTOP = KSTOP+1
         GO TO 110
      ENDIF
      IF(IPAROK.EQ.2) THEN
         IPARAM = MOPPAR(IPARAM)
      ENDIF
C     NEW ENTRY IN COMMON BLOCK.
      IF(NNPAR.LT.LMP) THEN
         NNPAR = NNPAR+1
         X (NNPAR) = PARAM
         NF(NNPAR) = IPARAM
         NS(NNPAR) = IELMNT
      ELSE
         KSTOP1 = KSTOP1+1
         GO TO 20
      ENDIF
C     NEW(IELMNT) = IELMNT
C     PRINTING SECTION.
      IF(IOPPRT.GT.-5) THEN
         IF(IPARAM.LE.72) THEN
            WRITE(NB6,540) PARTYP(IPARAM),ELEMNT(IELMNT),'  ',PARAM
         ELSE
            DO 120 I=73,LMPAR
            IF(MOPPAR(I).EQ.IPARAM) THEN
               WRITE(NB6,540) PARTYP(I),ELEMNT(IELMNT),'  ',PARAM
               GO TO 110
            ENDIF
  120       CONTINUE
            WRITE(NB6,550) IPARAM,ELEMNT(IELMNT),'  ',PARAM
         ENDIF
      ENDIF
      GO TO 110
C     *
C *** INPUT FINISHED.
  200 CONTINUE
      CLOSE(NB14)
      IF(KSTOP.GT.0 .OR. KSTOP1.GT.0) GO TO 300
      RETURN
C     *
C *** ERROR EXIT.
  300 CONTINUE
      WRITE(NB6,560)
      IF(KSTOP .GT.0) WRITE(NB6,570) KSTOP
      IF(KSTOP1.GT.0) WRITE(NB6,580) KSTOP1,LMP
      STOP 'EXTERN'
C
  500 FORMAT(///1X,'EXTERNAL PARAMETERS ARE READ FROM FILE NB14.',
     1       // 1X,' PARAMETER TYPE     ELEMENT(S)     PARAMETER'/)
  510 FORMAT(A)
  520 FORMAT(1X,A,'  FAULTY INPUT FOR PARAMETER TYPE')
  530 FORMAT(1X,A,'  FAULTY INPUT FOR ELEMENT')
  540 FORMAT(6X,A5,12X,A2,2X,A2,F16.8)
  550 FORMAT(6X,I3,14X,A2,2X,A2,F16.8)
  560 FORMAT(//1X,'ERROR WHEN READING EXTERNAL PARAMETER FILE'/)
  570 FORMAT(  1X,'NUMBER OF ILL-DEFINED PARAMETERS',I10/)
  580 FORMAT(//1X,'TOO MANY EXTERNAL PARAMETERS.',I5,
     1       / 1X,'NUMBER =',I5,5X,'MAXIMUM =',I5/)
  600 FORMAT(2I3,3X,F15.8)
      END
