BEGIN {
  freqline = 0 ;
  vibline = -1000 ;
  cartgradtable=-1000 ;
  intgradtable=-1000 ;
  ciline=-1000 ;
  shieldingsline=-1000 ;
  itoten = 0 ;
  iheat = 0 ;
  ioptcyc = 0 ;
  ihdlccyc = 0 ;
  ifinaltoten = 0 ;
  ifinalheat = 0 ;
  icitoten = 0 ;
  iciheat = 0 ;
  iroomheat = 0 ;
  izeroheat = 0 ;
  iintgrd = 0 ;
  icartgrd = 0 ;
  iqmmmgrd = 0 ;
  iintgrdopt = 0 ;
  icartgrdopt = 0 ;
  izpe = 0 ;
  iCIstate1 = 0 ;
  iCIstate2 = 0 ;
  iCIstate3 = 0 ;
  imeanpol = 0 ;
  ifirsthpol = 0 ;
  isecondhpol = 0 ;
  imol = 0 ;
}
( NR == 1 ) {
  ntoten = split(toten,atoten,":") ;
  nheat = split(heat,aheat,":") ;
  noptcyc = split(optcyc,aoptcyc,":") ;
  nhdlccyc = split(hdlccyc,ahdlccyc,":") ;
  nfinaltoten = split(finaltoten,afinaltoten,":") ;
  nfinalheat = split(finalheat,afinalheat,":") ;
  ncitoten = split(citoten,acitoten,":") ;
  nciheat = split(ciheat,aciheat,":") ;
  ndeltae  = split(deltae,adeltae,":") ;
  ndeltae2 = split(deltae2,adeltae2,":") ;
  ndeltae3 = split(deltae3,adeltae3,":") ;
  ndipole  = split(dipole,adipole,":") ;
  ndipole2 = split(dipole2,adipole2,":") ;
  ndipole3 = split(dipole3,adipole3,":") ;
  nfr  = split(fr,afr,":") ;
  nfr2 = split(fr2,afr2,":") ;
  nfr3 = split(fr3,afr3,":") ;
  nfp  = split(fp,afp,":") ;
  nfp2 = split(fp2,afp2,":") ;
  nfp3 = split(fp3,afp3,":") ;
  nfrp  = split(frp,afrp,":") ;
  nfrp2 = split(frp2,afrp2,":") ;
  nfrp3 = split(frp3,afrp3,":") ;
  nrot  = split(rot,arot,":") ;
  nrot2 = split(rot2,arot2,":") ;
  nrot3 = split(rot3,arot3,":") ;
  npop  = split(pop,apop,":") ;
  npop2 = split(pop2,apop2,":") ;
  npop3 = split(pop3,apop3,":") ;
  ncalref  = split(calref,acalref,":") ;
  ncalref2 = split(calref2,acalref2,":") ;
  ncalref3 = split(calref3,acalref3,":") ;
  nroomheat = split(roomheat,aroomheat,":") ;
  nzeroheat = split(zeroheat,azeroheat,":") ;
  nintgrd = split(intgrd,aintgrd,":") ;
  ncartgrd = split(cartgrd,acartgrd,":") ;
  nqmmmgrd = split(qmmmgrd,aqmmmgrd,":") ;
  nintgrdopt = split(intgrdopt,aintgrdopt,":") ;
  ncartgrdopt = split(cartgrdopt,acartgrdopt,":") ;
  nzpe = split(zpe,azpe,":")
  nCIstate1 = split(CIstate1,aCIstate1,":")
  nCIstate2 = split(CIstate2,aCIstate2,":")
  nCIstate3 = split(CIstate3,aCIstate3,":")
  nmeanpol = split(meanpol,ameanpol,":") ;
  nfirsthpol = split(firsthpol,afirsthpol,":") ;
  nsecondhpol = split(secondhpol,asecondhpol,":") ;
  nstate  = (ndeltae  > ndipole)  ? ndeltae  : ndipole;
  nstate2 = (ndeltae2 > ndipole2) ? ndeltae2 : ndipole2;
  nstate3 = (ndeltae3 > ndipole3) ? ndeltae3 : ndipole3;
  ntrans  = nfr;
  ntrans2 = nfr2;
  ntrans3 = nfr3;
  if (nfp != ntrans || nfrp != ntrans || nrot != ntrans)
    print "WARNING: Inconsistent number of transitions in opt file.";
  if (nfp2 != ntrans2 || nfrp2 != ntrans2 || nrot2 != ntrans2)
    print "WARNING: Molecule 2: Inconsistent number of transitions in opt file.";
  if (nfp3 != ntrans3 || nfrp3 != ntrans3 || nrot3 != ntrans3)
    print "WARNING: Molecule 3: Inconsistent number of transitions in opt file.";
}
/STANDARD INPUT USING KEYWORDS/{ imol = imol + 1 ; next ; }
/     TOTAL ENERGY/{ calctoten[++itoten] = $3 ; next ; }
/SCF HEAT OF FORMATION/{ calcheat[++iheat] = $5 ; next ; }
/HDLC OPTIMIZATION FINISHED AFTER/{ calchdlccyc[++ihdlccyc] = $5 ; next ; }
/OPTIMIZATION FINISHED AFTER/{ calcoptcyc[++ioptcyc] = $4 ; next ; }
/FINAL TOTAL ENERGY *[0-9\.+-][0-9\.+-]* *KCAL\/MOL/{ calcfinaltoten[++ifinaltoten] = $4 ; next ; }
/FINAL HEAT OF FORMATION *[0-9\.+-][0-9\.+-]* *KCAL\/MOL/{ calcfinalheat[++ifinalheat] = $5 ; next ; }
/CI TOTAL ENERGY *[0-9\.+-][0-9\.+-]* *KCAL\/MOL/{ calccitoten[++icitoten] = $4 ; next ; }
/CI HEAT OF FORMATION *[0-9\.+-][0-9\.+-]* *KCAL\/MOL/{ calcciheat[++iciheat] = $5 ; next ; }
/State dipole moments:/{ getline; getline;
  if (imol == 1)
  {
    for (i=1; i<=nstate; i++)
    {
      getline;
      calcdeltae[i] = $4;
      calcdipole[i] = $9;
    };
  }
  else if (imol == 2)
  {
    for (i=1; i<=nstate2; i++)
    {
      getline;
      calcdeltae2[i] = $4;
      calcdipole2[i] = $9;
    };
  }
  else if (imol == 3)
  {
    for (i=1; i<=nstate3; i++)
    {
      getline;
      calcdeltae3[i] = $4;
      calcdipole3[i] = $9;
    };
  }
  else
  {
    bad = bad sprintf("Unexpected molecule %d.", imol);
  }
  next;
}
/Properties of transitions   1 /{ getline; getline;
  if (imol == 1)
  {
    for (i=1; i<=ntrans; i++)
    {
      getline;
      calcfr[i]  = $7;
      calcfp[i]  = $8;
      calcfrp[i] = $9;
      calcrot[i] = $10;
    };
  }
  else if (imol == 2)
  {
    for (i=1; i<=ntrans2; i++)
    {
      getline;
      calcfr2[i]  = $7;
      calcfp2[i]  = $8;
      calcfrp2[i] = $9;
      calcrot2[i] = $10;
    };
  }
  else if (imol == 3)
  {
    for (i=1; i<=ntrans3; i++)
    {
      getline;
      calcfr3[i]  = $7;
      calcfp3[i]  = $8;
      calcfrp3[i] = $9;
      calcrot3[i] = $10;
    };
  }
  else
  {
    bad = bad sprintf("Unexpected molecule %d.", imol);
  }
  next;
}
/Population analysis of state/{ getline; getline; getline; getline;
  if (imol == 1)
  {
    for (i=1; i<=npop; i++)
    {
      getline;
      calcpop[i] = $2;
    };
  }
  else if (imol == 2)
  {
    for (i=1; i<=npop2; i++)
    {
      getline;
      calcpop2[i] = $2;
    };
  }
  else if (imol == 3)
  {
    for (i=1; i<=npop3; i++)
    {
      getline;
      calcpop3[i] = $2;
    };
  }
  else
  {
    bad = bad sprintf("Unexpected molecule %d.", imol);
  }
  next;
}
/COMPARISON OF THEORETICAL AND EXPERIMENTAL RESULTS./{ getline; getline; getline; getline;
  if (imol == 1)
  {
    for (i=1; i<=ncalref; i++)
    {
      getline;
      calccalref[i] = $6;
    };
  }
  else if (imol == 2)
  {
    for (i=1; i<=ncalref2; i++)
    {
      getline;
      calccalref2[i] = $6;
    };
  }
  else if (imol == 3)
  {
    for (i=1; i<=ncalref3; i++)
    {
      getline;
      calccalref3[i] = $6;
    };
  }
  else
  {
    bad = bad sprintf("Unexpected molecule %d.", imol);
  }
  next;
}
/HEAT OF FORMATION                 AT 298 K/{ calcroomheat[++iroomheat] = $7 ; next ; }
/HEAT OF FORMATION                 AT   0 K/{ calczeroheat[++izeroheat] = $7 ; next ; }
/FINAL OPTIMIZED GRADIENT NORM/{ calcintgrdopt[++iintgrdopt] = $5 ; next ; }
/   GRADIENT NORM  /{ calcintgrd[++iintgrd] = $3 ; next ; }
/FINAL CARTESIAN GRADIENT NORM/{ calccartgrdopt[++icartgrdopt] = $5 ; next ; }
/   CARTESIAN GRADIENT NORM  /{ calccartgrd[++icartgrd] = $4 ; next ; }
/   CARTESIAN QM\+MM GRADIENT NORM  /{ calcqmmmgrd[++iqmmmgrd] = $5 ; next ; }
/MEMORY REQUIRED DURING DERIVATIVES COMPUTATION IS/{ calcfreqmem = $7 ; next ; }
/MEMORY REQUIRED    :[ 0-9]* WORDS/{ calcfreqmem = $4 ; next ; }
/WORDS OF DISK SPACE WILL BE USED/{ calcfreqdisk = $1 ; next ; }
/DISK SPACE REQUIRED:[ 0-9]* WORDS/{ calcfreqdisk = $4 ; next ; }
/AVERAGE RESPON.E VECTORS PER VARIABLE/{ calcrespvec = $6 ; next ; }
/AVERAGE NUMBER OF USEFUL VECTORS PER VARIABLE/{ calcrespvec = $8 ; next ; }
/ *COORDINATES \(ANGSTROM\) *GRADIENTS \(KCAL\/\(MOL\*ANGSTROM\)\)/{ cartgradtable = NR + 3 ; next ; }
(cartgrda != "") && (NR == cartgradtable + cartgrda) { calccartgrdaX = $6 ; calccartgrdaY = $7 ; calccartgrdaZ = $8 ; next ; }
(cartgrdb != "") && (NR == cartgradtable + cartgrdb) { calccartgrdbX = $6 ; calccartgrdbY = $7 ; calccartgrdbZ = $8 ; next ; }
(cartgrdc != "") && (NR == cartgradtable + cartgrdc) { calccartgrdcX = $6 ; calccartgrdcY = $7 ; calccartgrdcZ = $8 ; next ; }
(cartgrdd != "") && (NR == cartgradtable + cartgrdd) { calccartgrddX = $6 ; calccartgrddY = $7 ; calccartgrddZ = $8 ; next ; }
/ GRADIENTS FOR GEOMETRY OPTIMIZATION/{ intgradtable = NR + 3 ; next ; }
(intgrda != "") && (NR == intgradtable + intgrda) { calcintgrdaX = $5 ; next ; }
(intgrdb != "") && (NR == intgradtable + intgrdb) { calcintgrdbX = $5 ; next ; }
(intgrdc != "") && (NR == intgradtable + intgrdc) { calcintgrdcX = $5 ; next ; }
(intgrdd != "") && (NR == intgradtable + intgrdd) { calcintgrddX = $5 ; next ; }
/ HEATS OF FORMATION \(IN KCAL\/MOL\) FOR THE CI STATES/{ ciline = NR + 3 ; next ; }
NR == ciline { calcCIstate1[++iCIstate1] = $1 ; calcCIstate2[++iCIstate2] = $2 ; calcCIstate3[++iCIstate3] = $3 ; next ; }
/MAXIMUM DEVIATION DUE TO SYMMETRIZATION/{ calcsymmdev = $6 ; next ; }
/EIGENVECTORS OF THE MASS-WEIGHTED CARTESIAN FORCE/{ freqline = NR + 5 ; next ; }
NR == freqline { calcfreq1 = $1 ; calcfreq2 = $2 ; calcfreq3 = $3 ; calcfreq4 = $4 ; 
calcfreq5 = $5 ; calcfreq6 = $6 ; calcfreq7 = $7 ; calcfreq8 = $8 ; 
calcfreq9 = $9 ; calcfreq10 = $10 ; next ; }
/ZERO-POINT VIBRATIONAL ENERGY/{ calczpe[++izpe] = $4 ; next ; }
/VIBRATIONAL INTENSITIES/{ vibline = NR + 4 - 7 ; next }
(modea != "") && (NR == (vibline + modea)) { calcwavea = $2 ; calcira = $3 ; next ; }
(modeb != "") && (NR == (vibline + modeb)) { calcwaveb = $2 ; calcirb = $3 ; next ; }
(modec != "") && (NR == (vibline + modec)) { calcwavec = $2 ; calcirc = $3 ; next ; }
(moded != "") && (NR == (vibline + moded)) { calcwaved = $2 ; calcird = $3 ; next ; }
/THREE-CENTER TERMS INCLUDED FOR .* ATOM PAIRS/{ calcNMRpairs = $5 ; next ; }
/ *X *Y *Z *SHIELDING *ANISOTROPY *ASYMMETRY *SIGMA11 *SIGMA22 *SIGMA33/{ shieldingsline = NR + 1 ; next ; }
(NMRatomA != "") && (NR == shieldingsline + NMRatomA) { calcNMRtypeA = $2 ; calcNMRisoA = $6 ; calcNMRs1A = $9 ; calcNMRs2A = $10 ; calcNMRs3A = $11 ; next ;}
(NMRatomB != "") && (NR == shieldingsline + NMRatomB) { calcNMRtypeB = $2 ; calcNMRisoB = $6 ; calcNMRs1B = $9 ; calcNMRs2B = $10 ; calcNMRs3B = $11 ; next ;}
(NMRatomC != "") && (NR == shieldingsline + NMRatomC) { calcNMRtypeC = $2 ; calcNMRisoC = $6 ; calcNMRs1C = $9 ; calcNMRs2C = $10 ; calcNMRs3C = $11 ; next ;}
(NMRatomD != "") && (NR == shieldingsline + NMRatomD) { calcNMRtypeD = $2 ; calcNMRisoD = $6 ; calcNMRs1D = $9 ; calcNMRs2D = $10 ; calcNMRs3D = $11 ; next ;}
(NMRatomE != "") && (NR == shieldingsline + NMRatomE) { calcNMRtypeE = $2 ; calcNMRisoE = $6 ; calcNMRs1E = $9 ; calcNMRs2E = $10 ; calcNMRs3E = $11 ; next ;}
(NMRatomF != "") && (NR == shieldingsline + NMRatomF) { calcNMRtypeF = $2 ; calcNMRisoF = $6 ; calcNMRs1F = $9 ; calcNMRs2F = $10 ; calcNMRs3F = $11 ; next ;}
(NMRatomG != "") && (NR == shieldingsline + NMRatomG) { calcNMRtypeG = $2 ; calcNMRisoG = $6 ; calcNMRs1G = $9 ; calcNMRs2G = $10 ; calcNMRs3G = $11 ; next ;}
(NMRatomH != "") && (NR == shieldingsline + NMRatomH) { calcNMRtypeH = $2 ; calcNMRisoH = $6 ; calcNMRs1H = $9 ; calcNMRs2H = $10 ; calcNMRs3H = $11 ; next ;}
/MEAN POLARIZABILITY/{ calcmeanpol[++imeanpol] = $4 ; next ; }
/FIRST HYPERPOLARIZABILITY/{ calcfirsthpol[++ifirsthpol] = $5 ; next ; }
/SECOND HYPERPOLARIZABILITY/{ calcsecondhpol[++isecondhpol] = $5 ; next ; }
END { 
  toltoten      = 0.000011  ; # kcal/mol
  tolheat       = 0.000011  ; # kcal/mol
  tolroomheat   = 0.0005    ; # kcal/mol
  tolzeroheat   = 0.0005    ; # kcal/mol
  tolfinaltoten = 0.0001    ; # kcal/mol
  tolfinalheat  = 0.0001    ; # kcal/mol
  tolcitoten    = 0.000011  ; # kcal/mol
  tolciheat     = 0.000011  ; # kcal/mol
  toldipole     = 0.0000011 ; # Debye
  toldeltae     = 0.0000011 ; # eV
  tolfr         = 0.0000011 ; # Debye
  tolfp         = 0.0000011 ; # Debye
  tolfrp        = 0.0000011 ; # Debye
  tolrot        = 0.0000011 ; # Debye
  tolpop        = 0.0000011 ; # au
  tolcalref     = 0.0001    ; # ?
  tolgrad       = 0.00001   ; # kcal/mol/angstrom
  tolgradnrm    = 0.000031  ; # kcal/mol/angstrom
  tolgradoptnrm = 0.02      ; # kcal/mol/angstrom
  tolcycl       = 1         ;
  tolhdlccycl   = 1         ;
  tolmem        = 0         ;
  tolrespvec    = 0.1       ;
  tolsymmdev    = 0.00003   ; # ?
  tolfreq       = 0.02      ; # cm^-1
  tolzpe        = 0.0005    ; # kcal/mol
  tolir         = 0.1       ; # kcal/mol
  tolCI         = 0.0001    ; # kcal/mol
  tolpairs      = 0         ;
  tolshift      = 0.005     ; # ppm
  tolmeanpol    = 0.001     ; # A**3
  tolfirsthpol  = 0.001     ; # 1.0D-30 ESU
  tolsecondhpol = 0.001     ; # 1.0D-36 ESU

  if( relaxfreq ){
    tolsymmdev *= 10 ;
    tolfreq    *= 10 ;
    tolzpe     *= 10 ;
    tolir      *= 10 ;
  }
  bad        = "" ;
  for ( i = 1 ; i <= ntoten ; i++ ) {
    if( atoten[i] != "" && (calctoten[i] < atoten[i] - toltoten ||
                           calctoten[i] > atoten[i] + toltoten ) )
    {
      if (i > 1)  bad = bad sprintf("Molecule %d: ", i);
      bad = bad sprintf("Wrong total energy: %7.5f instead of %7.5f\n",
                        calctoten[i],atoten[i]) ;
    }
  }
  for ( i = 1 ; i <= nheat ; i++ ) {
    if( aheat[i] != "" && (calcheat[i] < aheat[i] - tolheat ||
			   calcheat[i] > aheat[i] + tolheat ) )
    {
      if (i > 1)  bad = bad sprintf("Molecule %d: ", i);
      bad = bad sprintf("Wrong SCF heat of formation: %7.5f instead of %7.5f\n",
			calcheat[i],aheat[i]) ;
    }
  }
  for ( i = 1 ; i <= noptcyc ; i++ ) {
    if( aoptcyc[i] != "" && (calcoptcyc[i] < aoptcyc[i] - tolcyc ||
			     calcoptcyc[i] > aoptcyc[i] + tolcyc) ){
      bad = bad sprintf("Wrong number of optimization cycles: %d instead of %d\n",
			calcoptcyc[i],aoptcyc[i]) ;
    }
  }
  for ( i = 1 ; i <= nhdlccyc ; i++ ) {
    if( ahdlccyc[i] != "" && (calchdlccyc[i] < ahdlccyc[i] - tolhdlccyc ||
			      calchdlccyc[i] > ahdlccyc[i] + tolhdlccyc) ){
      bad = bad sprintf("Wrong number of hdlc optimization cycles: %d instead of %d\n",
			calchdlccyc[i],ahdlccyc[i]) ;
    }
  }
  for ( i = 1 ; i <= nfinaltoten ; i++ ) {
    if( afinaltoten[i] != "" && (calcfinaltoten[i] < afinaltoten[i] - tolfinaltoten ||
                                calcfinaltoten[i] > afinaltoten[i] + tolfinaltoten ) ){
      bad = bad sprintf("Wrong final total energy: %7.5f instead of %7.5f\n",
                        calcfinaltoten[i],afinaltoten[i]) ;
    }
  }
  for ( i = 1 ; i <= nfinalheat ; i++ ) {
    if( afinalheat[i] != "" && (calcfinalheat[i] < afinalheat[i] - tolfinalheat ||
				calcfinalheat[i] > afinalheat[i] + tolfinalheat ) ){
      bad = bad sprintf("Wrong final heat of formation: %7.5f instead of %7.5f\n",
			calcfinalheat[i],afinalheat[i]) ;
    }
  }
  for ( i = 1 ; i <= ncitoten ; i++ ) {
    if( acitoten[i] != "" && (calccitoten[i] < acitoten[i] - tolcitoten ||
                             calccitoten[i] > acitoten[i] + tolcitoten ) )
    {
      if (i > 1)  bad = bad sprintf("Molecule %d: ", i);
      bad = bad sprintf("Wrong CI total energy: %7.5f instead of %7.5f\n",
                        calccitoten[i],acitoten[i]) ;
    }
  }
  for ( i = 1 ; i <= nciheat ; i++ ) {
    if( aciheat[i] != "" && (calcciheat[i] < aciheat[i] - tolciheat ||
			     calcciheat[i] > aciheat[i] + tolciheat ) )
    {
      if (i > 1)  bad = bad sprintf("Molecule %d: ", i);
      bad = bad sprintf("Wrong CI heat of formation: %7.5f instead of %7.5f\n",
			calcciheat[i],aciheat[i]) ;
    }
  }
  for (i=1; i<=nstate; i++)
  {
    if (calcdeltae[i] == "")
    {
      bad = bad sprintf("Energy of state %d does not exist.\n", i);
    }
    else
    {
      if (calcdeltae[i] < adeltae[i] - toldeltae || calcdeltae[i] > adeltae[i] + toldeltae)
        bad = bad sprintf("Wrong energy of state %d: %8.6f instead of %8.6f\n",
                          i, calcdeltae[i], adeltae[i]);
    }
    if (adipole[i] != "")
    {
      if (calcdipole[i] == "")
      {
	bad = bad sprintf("State dipole moment %d does not exist.\n", i);
      }
      else
      {
	if (calcdipole[i] < adipole[i] - toldipole || calcdipole[i] > adipole[i] + toldipole)
	  bad = bad sprintf("Wrong state dipole moment %d: %8.6f instead of %8.6f\n",
			    i, calcdipole[i], adipole[i]);
      }
    }
  }
  for (i=1; i<=nstate2; i++)
  {
    if (calcdeltae2[i] == "")
    {
      bad = bad sprintf("Molecule 2: Energy of state %d does not exist.\n", i);
    }
    else
    {
      if (calcdeltae2[i] < adeltae2[i] - toldeltae || calcdeltae2[i] > adeltae2[i] + toldeltae)
        bad = bad sprintf("Molecule 2: Wrong energy of state %d: %8.6f instead of %8.6f\n",
                          i, calcdeltae2[i], adeltae2[i]);
    }
    if (adipole2[i] != "")
    {
      if (calcdipole2[i] == "")
      {
	bad = bad sprintf("Molecule 2: State dipole moment %d does not exist.\n", i);
      }
      else
      {
	if (calcdipole2[i] < adipole2[i] - toldipole || calcdipole2[i] > adipole2[i] + toldipole)
	  bad = bad sprintf("Molecule 2: Wrong state dipole moment %d: %8.6f instead of %8.6f\n",
			    i, calcdipole2[i], adipole2[i]);
      }
    }
  }
  for (i=1; i<=nstate3; i++)
  {
    if (calcdeltae3[i] == "")
    {
      bad = bad sprintf("Molecule 3: Energy of state %d does not exist.\n", i);
    }
    else
    {
      if (calcdeltae3[i] < adeltae3[i] - toldeltae || calcdeltae3[i] > adeltae3[i] + toldeltae)
        bad = bad sprintf("Molecule 3: Wrong energy of state %d: %8.6f instead of %8.6f\n",
                          i, calcdeltae3[i], adeltae3[i]);
    }
    if (adipole3[i] != "")
    {
      if (calcdipole3[i] == "")
      {
	bad = bad sprintf("Molecule 3: State dipole moment %d does not exist.\n", i);
      }
      else
      {
	if (calcdipole3[i] < adipole3[i] - toldipole || calcdipole3[i] > adipole3[i] + toldipole)
	  bad = bad sprintf("Molecule 3: Wrong state dipole moment %d: %8.6f instead of %8.6f\n",
			    i, calcdipole3[i], adipole3[i]);
      }
    }
  }
  for (i=1; i<=ntrans; i++)
  {
    if (calcfr[i] == "")
    {
      bad = bad sprintf("Dipole-length oscillator strength 1->%d does not exist.\n", i+1);
    }
    else
    {
      if (calcfr[i] < afr[i] - tolfr || calcfr[i] > afr[i] + tolfr)
        bad = bad sprintf("Wrong dipole-length oscillator strength 1->%d: %8.6f instead of %8.6f\n",
                          i+1, calcfr[i], afr[i]);
    }
    if (calcfp[i] == "")
    {
      bad = bad sprintf("Dipole-velocity oscillator strength 1->%d does not exist.\n", i+1);
    }
    else
    {
      if (calcfp[i] < afp[i] - tolfp || calcfp[i] > afp[i] + tolfp)
        bad = bad sprintf("Wrong dipole-velocity oscillator strength 1->%d: %8.6f instead of %8.6f\n",
                          i+1, calcfp[i], afp[i]);
    }
    if (calcfrp[i] == "")
    {
      bad = bad sprintf("Mixed oscillator strength 1->%d does not exist.\n", i+1);
    }
    else
    {
      if (calcfrp[i] < afrp[i] - tolfrp || calcfrp[i] > afrp[i] + tolfrp)
        bad = bad sprintf("Wrong mixed oscillator strength 1->%d: %8.6f instead of %8.6f\n",
                          i+1, calcfrp[i], afrp[i]);
    }
    if (calcrot[i] == "")
    {
      bad = bad sprintf("Rotational strength 1->%d does not exist.\n", i+1);
    }
    else
    {
      if (calcrot[i] < arot[i] - tolrot || calcrot[i] > arot[i] + tolrot)
        bad = bad sprintf("Wrong rotational strength 1->%d: %8.6f instead of %8.6f\n",
                          i+1, calcrot[i], arot[i]);
    }
  }
  for (i=1; i<=ntrans2; i++)
  {
    if (calcfr2[i] == "")
    {
      bad = bad sprintf("Molecule 2: Dipole-length oscillator strength 1->%d does not exist.\n", i+1);
    }
    else
    {
      if (calcfr2[i] < afr2[i] - tolfr || calcfr2[i] > afr2[i] + tolfr)
        bad = bad sprintf("Molecule 2: Wrong dipole-length oscillator strength 1->%d: %8.6f instead of %8.6f\n",
                          i+1, calcfr2[i], afr2[i]);
    }
    if (calcfp2[i] == "")
    {
      bad = bad sprintf("Molecule 2: Dipole-velocity oscillator strength 1->%d does not exist.\n", i+1);
    }
    else
    {
      if (calcfp2[i] < afp2[i] - tolfp || calcfp2[i] > afp2[i] + tolfp)
        bad = bad sprintf("Molecule 2: Wrong dipole-velocity oscillator strength 1->%d: %8.6f instead of %8.6f\n",
                          i+1, calcfp2[i], afp2[i]);
    }
    if (calcfrp2[i] == "")
    {
      bad = bad sprintf("Molecule 2: Mixed oscillator strength 1->%d does not exist.\n", i+1);
    }
    else
    {
      if (calcfrp2[i] < afrp2[i] - tolfrp || calcfrp2[i] > afrp2[i] + tolfrp)
        bad = bad sprintf("Molecule 2: Wrong mixed oscillator strength 1->%d: %8.6f instead of %8.6f\n",
                          i+1, calcfrp2[i], afrp2[i]);
    }
    if (calcrot2[i] == "")
    {
      bad = bad sprintf("Molecule 2: Rotational strength 1->%d does not exist.\n", i+1);
    }
    else
    {
      if (calcrot2[i] < arot2[i] - tolrot || calcrot2[i] > arot2[i] + tolrot)
        bad = bad sprintf("Molecule 2: Wrong rotational strength 1->%d: %8.6f instead of %8.6f\n",
                          i+1, calcrot2[i], arot2[i]);
    }
  }
  for (i=1; i<=ntrans3; i++)
  {
    if (calcfr3[i] == "")
    {
      bad = bad sprintf("Molecule 3: Dipole-length oscillator strength 1->%d does not exist.\n", i+1);
    }
    else
    {
      if (calcfr3[i] < afr3[i] - tolfr || calcfr3[i] > afr3[i] + tolfr)
        bad = bad sprintf("Molecule 3: Wrong dipole-length oscillator strength 1->%d: %8.6f instead of %8.6f\n",
                          i+1, calcfr3[i], afr3[i]);
    }
    if (calcfp3[i] == "")
    {
      bad = bad sprintf("Molecule 3: Dipole-velocity oscillator strength 1->%d does not exist.\n", i+1);
    }
    else
    {
      if (calcfp3[i] < afp3[i] - tolfp || calcfp3[i] > afp3[i] + tolfp)
        bad = bad sprintf("Molecule 3: Wrong dipole-velocity oscillator strength 1->%d: %8.6f instead of %8.6f\n",
                          i+1, calcfp3[i], afp3[i]);
    }
    if (calcfrp3[i] == "")
    {
      bad = bad sprintf("Molecule 3: Mixed oscillator strength 1->%d does not exist.\n", i+1);
    }
    else
    {
      if (calcfrp3[i] < afrp3[i] - tolfrp || calcfrp3[i] > afrp3[i] + tolfrp)
        bad = bad sprintf("Molecule 3: Wrong mixed oscillator strength 1->%d: %8.6f instead of %8.6f\n",
                          i+1, calcfrp3[i], afrp3[i]);
    }
    if (calcrot3[i] == "")
    {
      bad = bad sprintf("Molecule 3: Rotational strength 1->%d does not exist.\n", i+1);
    }
    else
    {
      if (calcrot3[i] < arot3[i] - tolrot || calcrot3[i] > arot3[i] + tolrot)
        bad = bad sprintf("Molecule 3: Wrong rotational strength 1->%d: %8.6f instead of %8.6f\n",
                          i+1, calcrot3[i], arot3[i]);
    }
  }
  for (i=1; i<=npop; i++)
  {
    if (calcpop[i] == "")
    {
      bad = bad sprintf("Charge of atom %d does not exist.\n", i);
    }
    else
    {
      if (calcpop[i] < apop[i] - tolpop || calcpop[i] > apop[i] + tolpop)
        bad = bad sprintf("Wrong charge of atom %d: %7.4f instead of %7.4f\n",
                          i, calcpop[i], apop[i]);
    }
  }
  for (i=1; i<=npop2; i++)
  {
    if (calcpop2[i] == "")
    {
      bad = bad sprintf("Molecule 2: Charge of atom %d does not exist.\n", i);
    }
    else
    {
      if (calcpop2[i] < apop2[i] - tolpop || calcpop2[i] > apop2[i] + tolpop)
        bad = bad sprintf("Molecule 2: Wrong charge of atom %d: %7.4f instead of %7.4f\n",
                          i, calcpop2[i], apop2[i]);
    }
  }
  for (i=1; i<=npop3; i++)
  {
    if (calcpop3[i] == "")
    {
      bad = bad sprintf("Molecule 3: Charge of atom %d does not exist.\n", i);
    }
    else
    {
      if (calcpop3[i] < apop3[i] - tolpop || calcpop3[i] > apop3[i] + tolpop)
        bad = bad sprintf("Molecule 3: Wrong charge of atom %d: %7.4f instead of %7.4f\n",
                          i, calcpop3[i], apop3[i]);
    }
  }
  for (i=1; i<=ncalref; i++)
  {
    if (calccalref[i] == "")
    {
      bad = bad sprintf("Reference value #%d does not exist.\n", i);
    }
    else
    {
      if (calccalref[i] < acalref[i] - tolcalref || calccalref[i] > acalref[i] + tolcalref)
        bad = bad sprintf("Wrong reference value #%d: %7.5f instead of %7.5f\n",
                          i, calccalref[i], acalref[i]);
    }
  }
  for (i=1; i<=ncalref2; i++)
  {
    if (calccalref2[i] == "")
    {
      bad = bad sprintf("Molecule 2: Reference value #%d does not exist.\n", i);
    }
    else
    {
      if (calccalref2[i] < acalref2[i] - tolcalref || calccalref2[i] > acalref2[i] + tolcalref)
        bad = bad sprintf("Molecule 2: Wrong reference value #%d: %7.5f instead of %7.5f\n",
                          i, calccalref2[i], acalref2[i]);
    }
  }
  for (i=1; i<=ncalref3; i++)
  {
    if (calccalref3[i] == "")
    {
      bad = bad sprintf("Molecule 3: Reference value #%d does not exist.\n", i);
    }
    else
    {
      if (calccalref3[i] < acalref3[i] - tolcalref || calccalref3[i] > acalref3[i] + tolcalref)
        bad = bad sprintf("Molecule 3: Wrong reference value #%d: %7.5f instead of %7.5f\n",
                          i, calccalref3[i], acalref3[i]);
    }
  }
  for ( i = 1 ; i <= nroomheat ; i++ ) {
    if( aroomheat[i] != "" && (calcroomheat[i] < aroomheat[i] - tolroomheat ||
                               calcroomheat[i] > aroomheat[i] + tolroomheat ) ){
      bad = bad sprintf("Wrong heat of formation at 298 K: %g instead of %g\n",
                        calcroomheat[i],aroomheat[i]) ;
    }
  }  
  for ( i = 1 ; i <= nzeroheat ; i++ ) {
    if( azeroheat[i] != "" && (calczeroheat[i] < azeroheat[i] - tolzeroheat ||
			       calczeroheat[i] > azeroheat[i] + tolzeroheat ) ){
      bad = bad sprintf("Wrong heat of formation at 0 K: %g instead of %g\n",
			calczeroheat[i],azeroheat[i]) ;
    }
  }
  for ( i = 1 ; i <= ncartgrd ; i++ )
  {
    if( acartgrd[i] != "" && (calccartgrd[i] < acartgrd[i] - tolgradnrm ||
			      calccartgrd[i] > acartgrd[i] + tolgradnrm) )
    {
      if (i > 1)  bad = bad sprintf("Molecule %d: ", i);
      bad = bad sprintf("Wrong cartesian gradient norm: %7.5f instead of %7.5f\n",
			calccartgrd[i],acartgrd[i]) ;
    }
  }
  for ( i = 1 ; i <= nintgrd ; i++ )
  {
    if( aintgrd[i] != "" && (calcintgrd[i] < aintgrd[i] - tolgradnrm ||
			     calcintgrd[i] > aintgrd[i] + tolgradnrm) )
    {
      if (i > 1)  bad = bad sprintf("Molecule %d: ", i);
      bad = bad sprintf("Wrong gradient norm: %7.5f instead of %7.5f\n",
			calcintgrd[i],aintgrd[i]) ;
    }
  }
  for ( i = 1 ; i <= nqmmmgrd ; i++ )
  {
    if( aqmmmgrd[i] != "" && (calcqmmmgrd[i] < aqmmmgrd[i] - tolgradnrm ||
			      calcqmmmgrd[i] > aqmmmgrd[i] + tolgradnrm) )
    {
      if (i > 1)  bad = bad sprintf("Molecule %d: ", i);
      bad = bad sprintf("Wrong cartesian qm+mm gradient norm: %7.5f instead of %7.5f\n",
			calcqmmmgrd[i],aqmmmgrd[i]) ;
    }
  }
  for ( i = 1 ; i <= nintgrdopt ; i++ ) {
    if( aintgrdopt[i] != "" && (calcintgrdopt[i] < aintgrdopt[i] - tolgradoptnrm ||
				calcintgrdopt[i] > aintgrdopt[i] + tolgradoptnrm) ){
      bad = bad sprintf("Wrong final internal gradient norm after optimization: %7.5f instead of %7.5f\n",
			calcintgrdopt[i],aintgrdopt[i]) ;
    }
  }
  for ( i = 1 ; i <= ncartgrdopt ; i++ ) {
    if( acartgrdopt[i] != "" && (calccartgrdopt[i] < acartgrdopt[i] - tolgradoptnrm ||
				 calccartgrdopt[i] > acartgrdopt[i] + tolgradoptnrm) ){
      bad = bad sprintf("Wrong final cartesian gradient norm after optimization: %7.5f instead of %7.5f\n",
			calccartgrdopt[i],acartgrdopt[i]) ;
    }
  }
  if( freqmem != "" && (calcfreqmem < freqmem - tolmem || calcfreqmem > freqmem + tolmem) ){
    bad = bad sprintf("Wrong memory requirements in CPHF: %d instead of %d\n",calcfreqmem,freqmem) ;
  }
  if( freqdisk != "" && (calcfreqdisk < freqdisk - tolmem || calcfreqdisk > freqdisk + tolmem) ){
    bad = bad sprintf("Wrong disk space requirements in CPHF: %d instead of %d\n",calcfreqdisk,freqdisk) ;
  }
  if( respvec != "" && (calcrespvec < respvec - tolrespvec || calcrespvec > respvec + tolrespvec) ){
    bad = bad sprintf("Wrong average number of responce vectors: %g instead of %g\n",calcrespvec,respvec) ;
  }
  if( cartgrda != "" && (calccartgrdaX < cartgrdaX - tolgrad || calccartgrdaX > cartgrdaX + tolgrad || \
			 calccartgrdaY < cartgrdaY - tolgrad || calccartgrdaY > cartgrdaY + tolgrad || \
			 calccartgrdaZ < cartgrdaZ - tolgrad || calccartgrdaZ > cartgrdaZ + tolgrad) ){
    bad = bad sprintf("Wrong cartesian gradient on atom %d: (%g,%g,%g) instead of (%g,%g,%g)\n", \
		      cartgrda,calccartgrdaX,calccartgrdaY,calccartgrdaZ,cartgrdaX,cartgrdaY,cartgrdaZ) ;
  }
  if( cartgrdb != "" && (calccartgrdbX < cartgrdbX - tolgrad || calccartgrdbX > cartgrdbX + tolgrad || \
			 calccartgrdbY < cartgrdbY - tolgrad || calccartgrdbY > cartgrdbY + tolgrad || \
			 calccartgrdbZ < cartgrdbZ - tolgrad || calccartgrdbZ > cartgrdbZ + tolgrad) ){
    bad = bad sprintf("Wrong cartesian gradient on atom %d: (%g,%g,%g) instead of (%g,%g,%g)\n", \
		      cartgrdb,calccartgrdbX,calccartgrdbY,calccartgrdbZ,cartgrdbX,cartgrdbY,cartgrdbZ) ;
  }
  if( cartgrdc != "" && (calccartgrdcX < cartgrdcX - tolgrad || calccartgrdcX > cartgrdcX + tolgrad || \
			 calccartgrdcY < cartgrdcY - tolgrad || calccartgrdcY > cartgrdcY + tolgrad || \
			 calccartgrdcZ < cartgrdcZ - tolgrad || calccartgrdcZ > cartgrdcZ + tolgrad) ){
    bad = bad sprintf("Wrong cartesian gradient on atom %d: (%g,%g,%g) instead of (%g,%g,%g)\n", \
		      cartgrdc,calccartgrdcX,calccartgrdcY,calccartgrdcZ,cartgrdcX,cartgrdcY,cartgrdcZ) ;
  }
  if( cartgrdd != "" && (calccartgrddX < cartgrddX - tolgrad || calccartgrddX > cartgrddX + tolgrad || \
			 calccartgrddY < cartgrddY - tolgrad || calccartgrddY > cartgrddY + tolgrad || \
			 calccartgrddZ < cartgrddZ - tolgrad || calccartgrddZ > cartgrddZ + tolgrad) ){
    bad = bad sprintf("Wrong cartesian gradient on atom %d: (%g,%g,%g) instead of (%g,%g,%g)\n", \
		      cartgrdd,calccartgrddX,calccartgrddY,calccartgrddZ,cartgrddX,cartgrddY,cartgrddZ) ;
  }
  if( intgrda != "" && (calcintgrdaX < intgrdaX - tolgrad || calcintgrdaX > intgrdaX + tolgrad) ){
    bad = bad sprintf("Wrong internal gradient %d: %g instead of %g\n",intgrda,calcintgrdaX,intgrdaX) ;
  }
  if( intgrdb != "" && (calcintgrdbX < intgrdbX - tolgrad || calcintgrdbX > intgrdbX + tolgrad) ){
    bad = bad sprintf("Wrong internal gradient %d: %g instead of %g\n",intgrdb,calcintgrdbX,intgrdbX) ;
  }
  if( intgrdc != "" && (calcintgrdcX < intgrdcX - tolgrad || calcintgrdcX > intgrdcX + tolgrad) ){
    bad = bad sprintf("Wrong internal gradient %d: %g instead of %g\n",intgrdc,calcintgrdcX,intgrdcX) ;
  }
  if( intgrdd != "" && (calcintgrddX < intgrddX - tolgrad || calcintgrddX > intgrddX + tolgrad) ){
    bad = bad sprintf("Wrong internal gradient %d: %g instead of %g\n",intgrdd,calcintgrddX,intgrddX) ;
  }
  for ( i = 1 ; i <= nCIstate1 ; i++ ) {
    if( aCIstate1[i] != "" && (calcCIstate1[i] < aCIstate1[i] - tolCI || 
			       calcCIstate1[i] > aCIstate1[i] + tolCI) ){
      bad = bad sprintf("Wrong CI state 1 energy: %g instead of %g\n",
			calcCIstate1[i],aCIstate1[i]) ;
    }
  }
  for ( i = 1 ; i <= nCIstate2 ; i++ ) {
    if( aCIstate2[i] != "" && (calcCIstate2[i] < aCIstate2[i] - tolCI || 
			       calcCIstate2[i] > aCIstate2[i] + tolCI) ){
      bad = bad sprintf("Wrong CI state 2 energy: %g instead of %g\n",
			calcCIstate2[i],aCIstate2[i]) ;
    }
  }
  for ( i = 1 ; i <= nCIstate3 ; i++ ) {
    if( aCIstate3[i] != "" && (calcCIstate3[i] < aCIstate3[i] - tolCI || 
			       calcCIstate3[i] > aCIstate3[i] + tolCI) ){
      bad = bad sprintf("Wrong CI state 3 energy: %g instead of %g\n",
			calcCIstate3[i],aCIstate3[i]) ;
    }
  }
  if( (calcsymmdev+0) > tolsymmdev ){
    bad = bad sprintf("Force constants symmetrization deviation too large: %g\n",calcsymmdev) ;
  }
  if( freq1 != "" && (calcfreq1 < freq1 - tolfreq || calcfreq1 > freq1 + tolfreq) ){
    bad = bad sprintf("Vibrational frequency 1 is wrong: %7.5f instead of %7.5f\n",calcfreq1,freq1) ;
  }
  if( freq2 != "" && (calcfreq2 < freq2 - tolfreq || calcfreq2 > freq2 + tolfreq) ){
    bad = bad sprintf("Vibrational frequency 2 is wrong: %7.5f instead of %7.5f\n",calcfreq2,freq2) ;
  }
  if( freq3 != "" && (calcfreq3 < freq3 - tolfreq || calcfreq3 > freq3 + tolfreq) ){
    bad = bad sprintf("Vibrational frequency 3 is wrong: %7.5f instead of %7.5f\n",calcfreq3,freq3) ;
  }
  if( freq4 != "" && (calcfreq4 < freq4 - tolfreq || calcfreq4 > freq4 + tolfreq) ){
    bad = bad sprintf("Vibrational frequency 4 is wrong: %7.5f instead of %7.5f\n",calcfreq4,freq4) ;
  }
  if( freq5 != "" && (calcfreq5 < freq5 - tolfreq || calcfreq5 > freq5 + tolfreq) ){
    bad = bad sprintf("Vibrational frequency 5 is wrong: %7.5f instead of %7.5f\n",calcfreq5,freq5) ;
  }
  if( freq6 != "" && (calcfreq6 < freq6 - tolfreq || calcfreq6 > freq6 + tolfreq) ){
    bad = bad sprintf("Vibrational frequency 6 is wrong: %7.5f instead of %7.5f\n",calcfreq6,freq6) ;
  }
  if( freq7 != "" && (calcfreq7 < freq7 - tolfreq || calcfreq7 > freq7 + tolfreq) ){
    bad = bad sprintf("Vibrational frequency 7 is wrong: %7.5f instead of %7.5f\n",calcfreq7,freq7) ;
  }
  if( freq8 != "" && (calcfreq8 < freq8 - tolfreq || calcfreq8 > freq8 + tolfreq) ){
    bad = bad sprintf("Vibrational frequency 8 is wrong: %7.5f instead of %7.5f\n",calcfreq8,freq8) ;
  }
  if( freq9 != "" && (calcfreq9 < freq9 - tolfreq || calcfreq9 > freq9 + tolfreq) ){
    bad = bad sprintf("Vibrational frequency 9 is wrong: %7.5f instead of %7.5f\n",calcfreq9,freq9) ;
  }
  if( freq10 != "" && (calcfreq10 < freq10 - tolfreq || calcfreq10 > freq10 + tolfreq) ){
    bad = bad sprintf("Vibrational frequency 10 is wrong: %7.5f instead of %7.5f\n",calcfreq10,freq10) ;
  }
  for ( i = 1 ; i <= nzpe ; i++ ) {
    if( azpe[i] != "" && (calczpe[i] < azpe[i] - tolzpe ||
			  calczpe[i] > azpe[i] + tolzpe ) ){
      bad = bad sprintf("ZPE is wrong: %g instead of %g\n",
			calczpe[i],azpe[i]) ;
    }
  }
  if( modea != "" && (calcwavea < wavea - tolfreq || calcwavea > wavea + tolfreq || \
		      calcira   < ira   - tolir   || calcira   > ira   + tolir ) ){
    bad = bad sprintf("IR band %d is wrong: %g(%g) instead of %g(%g)\n",modea,calcwavea,calcira,wavea,ira) ;
  }
  if( modeb != "" && (calcwaveb < waveb - tolfreq || calcwaveb > waveb + tolfreq || \
		      calcirb   < irb   - tolir   || calcirb   > irb   + tolir ) ){
    bad = bad sprintf("IR band %d is wrong: %g(%g) instead of %g(%g)\n",modeb,calcwaveb,calcirb,waveb,irb) ;
  }
  if( modec != "" && (calcwavec < wavec - tolfreq || calcwavec > wavec + tolfreq || \
		      calcirc   < irc   - tolir   || calcirc   > irc   + tolir ) ){
    bad = bad sprintf("IR band %d is wrong: %g(%g) instead of %g(%g)\n",modec,calcwavec,calcirc,wavec,irc) ;
  }
  if( moded != "" && (calcwaved < waved - tolfreq || calcwaved > waved + tolfreq || \
		      calcird   < ird   - tolir   || calcird   > ird   + tolir ) ){
    bad = bad sprintf("IR band %d is wrong: %g(%g) instead of %g(%g)\n",moded,calcwaved,calcird,waved,ird) ;
  }
  if( NMRpairs != "" && (calcNMRpairs < NMRpairs - tolpairs || calcNMRpairs > NMRpairs + tolpairs) ){
    bad = bad sprintf("Number of atom pairs in NMR calculation is wrong: %d instead of %d\n",calcNMRpairs,NMRpairs) ;
  }
  if( NMRatomA != "" && (calcNMRtypeA != NMRtypeA || \
			 calcNMRisoA < NMRisoA - tolshift || calcNMRisoA > NMRisoA + tolshift || \
			 calcNMRs1A  < NMRs1A  - tolshift || calcNMRs1A  > NMRs1A  + tolshift || \
			 calcNMRs2A  < NMRs2A  - tolshift || calcNMRs2A  > NMRs2A  + tolshift || \
			 calcNMRs3A  < NMRs3A  - tolshift || calcNMRs3A  > NMRs3A  + tolshift ) ){
    bad = bad sprintf("Wrong NMR results: type = %s, iso = %g, s11 = %g, s22 = %g, s33 = %g instead of %s, %s, %g, %g, %g\n", \
		      calcNMRtypeA, calcNMRisoA, calcNMRs1A, calcNMRs2A, calcNMRs3A, NMRtypeA, NMRisoA, NMRs1A, NMRs2A, NMRs3A) ;
  }
  if( NMRatomB != "" && (calcNMRtypeB != NMRtypeB || \
			 calcNMRisoB < NMRisoB - tolshift || calcNMRisoB > NMRisoB + tolshift || \
			 calcNMRs1B  < NMRs1B  - tolshift || calcNMRs1B  > NMRs1B  + tolshift || \
			 calcNMRs2B  < NMRs2B  - tolshift || calcNMRs2B  > NMRs2B  + tolshift || \
			 calcNMRs3B  < NMRs3B  - tolshift || calcNMRs3B  > NMRs3B  + tolshift ) ){
    bad = bad sprintf("Wrong NMR results: type = %s, iso = %g, s11 = %g, s22 = %g, s33 = %g instead of %s, %s, %g, %g, %g\n", \
		      calcNMRtypeB, calcNMRisoB, calcNMRs1B, calcNMRs2B, calcNMRs3B, NMRtypeB, NMRisoB, NMRs1B, NMRs2B, NMRs3B) ;
  }
  if( NMRatomC != "" && (calcNMRtypeC != NMRtypeC || \
			 calcNMRisoC < NMRisoC - tolshift || calcNMRisoC > NMRisoC + tolshift || \
			 calcNMRs1C  < NMRs1C  - tolshift || calcNMRs1C  > NMRs1C  + tolshift || \
			 calcNMRs2C  < NMRs2C  - tolshift || calcNMRs2C  > NMRs2C  + tolshift || \
			 calcNMRs3C  < NMRs3C  - tolshift || calcNMRs3C  > NMRs3C  + tolshift ) ){
    bad = bad sprintf("Wrong NMR results: type = %s, iso = %g, s11 = %g, s22 = %g, s33 = %g instead of %s, %s, %g, %g, %g\n", \
		      calcNMRtypeC, calcNMRisoC, calcNMRs1C, calcNMRs2C, calcNMRs3C, NMRtypeC, NMRisoC, NMRs1C, NMRs2C, NMRs3C) ;
  }
  if( NMRatomD != "" && (calcNMRtypeD != NMRtypeD || \
			 calcNMRisoD < NMRisoD - tolshift || calcNMRisoD > NMRisoD + tolshift || \
			 calcNMRs1D  < NMRs1D  - tolshift || calcNMRs1D  > NMRs1D  + tolshift || \
			 calcNMRs2D  < NMRs2D  - tolshift || calcNMRs2D  > NMRs2D  + tolshift || \
			 calcNMRs3D  < NMRs3D  - tolshift || calcNMRs3D  > NMRs3D  + tolshift ) ){
    bad = bad sprintf("Wrong NMR results: type = %s, iso = %g, s11 = %g, s22 = %g, s33 = %g instead of %s, %s, %g, %g, %g\n", \
		      calcNMRtypeD, calcNMRisoD, calcNMRs1D, calcNMRs2D, calcNMRs3D, NMRtypeD, NMRisoD, NMRs1D, NMRs2D, NMRs3D) ;
  }
  if( NMRatomE != "" && (calcNMRtypeE != NMRtypeE || \
			 calcNMRisoE < NMRisoE - tolshift || calcNMRisoE > NMRisoE + tolshift || \
			 calcNMRs1E  < NMRs1E  - tolshift || calcNMRs1E  > NMRs1E  + tolshift || \
			 calcNMRs2E  < NMRs2E  - tolshift || calcNMRs2E  > NMRs2E  + tolshift || \
			 calcNMRs3E  < NMRs3E  - tolshift || calcNMRs3E  > NMRs3E  + tolshift ) ){
    bad = bad sprintf("Wrong NMR results: type = %s, iso = %g, s11 = %g, s22 = %g, s33 = %g instead of %s, %s, %g, %g, %g\n", \
		      calcNMRtypeE, calcNMRisoE, calcNMRs1E, calcNMRs2E, calcNMRs3E, NMRtypeE, NMRisoE, NMRs1E, NMRs2E, NMRs3E) ;
  }
  if( NMRatomF != "" && (calcNMRtypeF != NMRtypeF || \
			 calcNMRisoF < NMRisoF - tolshift || calcNMRisoF > NMRisoF + tolshift || \
			 calcNMRs1F  < NMRs1F  - tolshift || calcNMRs1F  > NMRs1F  + tolshift || \
			 calcNMRs2F  < NMRs2F  - tolshift || calcNMRs2F  > NMRs2F  + tolshift || \
			 calcNMRs3F  < NMRs3F  - tolshift || calcNMRs3F  > NMRs3F  + tolshift ) ){
    bad = bad sprintf("Wrong NMR results: type = %s, iso = %g, s11 = %g, s22 = %g, s33 = %g instead of %s, %s, %g, %g, %g\n", \
		      calcNMRtypeF, calcNMRisoF, calcNMRs1F, calcNMRs2F, calcNMRs3F, NMRtypeF, NMRisoF, NMRs1F, NMRs2F, NMRs3F) ;
  }
  if( NMRatomG != "" && (calcNMRtypeG != NMRtypeG || \
			 calcNMRisoG < NMRisoG - tolshift || calcNMRisoG > NMRisoG + tolshift || \
			 calcNMRs1G  < NMRs1G  - tolshift || calcNMRs1G  > NMRs1G  + tolshift || \
			 calcNMRs2G  < NMRs2G  - tolshift || calcNMRs2G  > NMRs2G  + tolshift || \
			 calcNMRs3G  < NMRs3G  - tolshift || calcNMRs3G  > NMRs3G  + tolshift ) ){
    bad = bad sprintf("Wrong NMR results: type = %s, iso = %g, s11 = %g, s22 = %g, s33 = %g instead of %s, %s, %g, %g, %g\n", \
		      calcNMRtypeG, calcNMRisoG, calcNMRs1G, calcNMRs2G, calcNMRs3G, NMRtypeG, NMRisoG, NMRs1G, NMRs2G, NMRs3G) ;
  }
  if( NMRatomH != "" && (calcNMRtypeH != NMRtypeH || \
			 calcNMRisoH < NMRisoH - tolshift || calcNMRisoH > NMRisoH + tolshift || \
			 calcNMRs1H  < NMRs1H  - tolshift || calcNMRs1H  > NMRs1H  + tolshift || \
			 calcNMRs2H  < NMRs2H  - tolshift || calcNMRs2H  > NMRs2H  + tolshift || \
			 calcNMRs3H  < NMRs3H  - tolshift || calcNMRs3H  > NMRs3H  + tolshift ) ){
    bad = bad sprintf("Wrong NMR results: type = %s, iso = %g, s11 = %g, s22 = %g, s33 = %g instead of %s, %s, %g, %g, %g\n", \
		      calcNMRtypeH, calcNMRisoH, calcNMRs1H, calcNMRs2H, calcNMRs3H, NMRtypeH, NMRisoH, NMRs1H, NMRs2H, NMRs3H) ;
  }
  for ( i = 1 ; i <= nmeanpol ; i++ ) {
    if ( ameanpol[i] != "" && ( calcmeanpol[i] < ameanpol[i] - tolmeanpol ||
				calcmeanpol[i] > ameanpol[i] + tolmeanpol ) ){
      bad = bad sprintf("Wrong mean polarizability: %g instead of %g\n",
			calcmeanpol[i],ameanpol[i]) ;
    }
  } 
  for ( i = 1 ; i <= nfirsthpol ; i++ ) {
    if ( afirsthpol[i] != "" && ( calcfirsthpol[i] < afirsthpol[i] - tolfirsthpol ||
				  calcfirsthpol[i] > afirsthpol[i] + tolfirsthpol ) ){
      bad = bad sprintf("Wrong first hyperpolarizability: %g instead of %g\n",
			calcfirsthpol[i],afirsthpol[i]) ;
    }
  } 
  for ( i = 1 ; i <= nsecondhpol ; i++ ) {
    if ( asecondhpol[i] != "" && ( calcsecondhpol[i] < asecondhpol[i] - tolsecondhpol ||
				   calcsecondhpol[i] > asecondhpol[i] + tolsecondhpol ) ){
      bad = bad sprintf("Wrong second hyperpolarizability: %g instead of %g\n",
			calcsecondhpol[i],asecondhpol[i]) ;
    }
  } 
  if( bad != "" ){ 
    printf "%s", bad ; 
    print bad > (name ".fail") ;
  } 
  else { 
    print "Passed" ;
    printf "" > (name ".pass") ;
  }
}
