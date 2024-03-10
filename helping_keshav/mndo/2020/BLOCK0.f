      BLOCKDATA BLOCK0
C     *
C     INITIALIZATION FOR CONSTANTS, LABELS, AND SYMMETRY DATA.
C     *
C     COMMON BLOCKS:  ATSYM , CONSTF, CONSTN, ELEMTS, INOPT1, INOPT2,
C     ISTOPE, MULTAB, NBFILE, PAR43,  QMNAME, SYMALL, SYMLAB, SYMLIM,
C     TXPROP.
C     *
C     SPECIAL CONVENTION.
C     ELEMENT 86 IS A CONNECTION ATOM FOR QM/MM TREATMENTS.
C     EXPERIMENTAL HEAT OF FORMATION AND MASS FROM METHYL (CH3).
C     *
      USE LIMIT, ONLY: LMZ, LMPR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*2  ELEMNT
      CHARACTER*3  NGROUP
      CHARACTER*3  NPTGRP
      CHARACTER*4  IRREP
      CHARACTER*7  METHOD
      CHARACTER*30 PROP1
      CHARACTER*31 PROP2
      CHARACTER*37 PROP3
      CHARACTER*10 UNITS
      COMMON
     ./ATSYM / IELSYM(21)
     ./CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./ELEMTS/ ELEMNT(107)
     ./INOPT1/ IN1(300)
     ./INOPT2/ IN2(300)
     ./ISTOPE/ CMS(LMZ),BMS(LMZ)
     ./MULTAB/ IPROD(8,8)
     ./NBFILE/ NBF(20)
     ./PAR43 / EHEATP,IPHEXP
     ./QMNAME/ METHOD(-23:6)
     ./SYMALL/ NPTGRP(59)
     ./SYMLAB/ NGROUP(7),IRREP(26)
     ./SYMLIM/ IRREPA(7),IRREPN(7)
     ./TXPROP/ PROP1(LMPR),PROP2(LMPR),PROP3(LMPR),UNITS(LMPR)
C *** DEFAULT OPTIONS.
      DATA IN1 / 300*0 /
      DATA IN2 / 300*0 /
C *** DEFAULT FILE NUMBERS.
      DATA NBF/ 1, 2, 3, 4, 5, 6, 7, 8, 9,10,
     1         11,12,13,14,15,16,17,18,19,20/
C *** SYMMETRY LABELS.
C     BEHAVIOR OF CARTESIAN COORDINATES UNDER SYMMETRY OPERATIONS.
      DATA IELSYM/-1, 1, 1, 1,-1, 1, 1, 1,-1,
     1             1,-1,-1,-1, 1,-1,-1,-1, 1,-1,-1,-1/
C *** CONSTANT CONVERSION FACTORS AND CUTOFFS.
C     SOME OF THESE FACTORS DO NOT REPRESENT THE MOST ACCURATE VALUES
C     AVAILABLE NOWADAYS. HOWEVER, THESE VALUES WERE USED WHEN DERIVING
C     THE SEMIEMPIRICAL PARAMETERS. THEREFORE THESE VALUES MUST BE KEPT
C     UNLESS ONE WANTS TO REDEFINE THE PARAMETERS.
      DATA A0    /  0.529167D0/
      DATA AFACT / 57.29577951308232D0/
      DATA EV    / 27.21D0/
      DATA EVCAL / 23.061D0/
      DATA PI    /  3.141592653589793D0/
      DATA W1    / 14.399D0/
      DATA W2    /  7.1995D0/
      DATA BIGEXP/ 50.0D0/
C *** NUMERICAL CONSTANTS.
      DATA ZERO  /  0.0D0/
      DATA ONE   /  1.0D0/
      DATA TWO   /  2.0D0/
      DATA THREE /  3.0D0/
      DATA FOUR  /  4.0D0/
      DATA PT5   /  0.5D0/
      DATA PT25  /  0.25D0/
C *** THE PROTON'S HEAT OF FORMATION
      DATA EHEATP/367.171D0/
      DATA IPHEXP/313.5873D0/
C *** SYMBOLS FOR CHEMICAL ELEMENTS.
      DATA ELEMNT/'H','He',
     1 'Li','Be','B ','C ','N ','O ','F ','Ne',
     2 'Na','Mg','Al','Si','P ','S ','Cl','Ar',
     3 'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu',
     4 'Zn','Ga','Ge','As','Se','Br','Kr',
     5 'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag',
     6 'Cd','In','Sn','Sb','Te','I ','Xe',
     7 'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy',
     8 'Ho','Er','Tm','Yb','Lu','Hf','Ta','W ','Re','Os','Ir','Pt',
     9 'Au','Hg','Tl','Pb','Bi','Po','At','**',
     A 'Fr','Ra','Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','XX',
     B 'Fm','Md','Cb','++','+ ','--','- ','TV'/
C *** NAMES OF THE AVAILABLE QUANTUM-CHEMICAL METHODS.
      DATA METHOD/ 'ODM3   ','ODM2   ','       ','       ','       ',
     1   '       ','PM3/tm ','       ','       ','       ','MNDO/dH',
     2   'AM1/d  ','       ','MNDO/d ','OM4    ','OM3    ','PM3    ',
     3   'OM2    ','OM1    ','MNDO   ','MNDO/H ','AM1    ','MNDOC  ',
     4   'MNDO   ','MINDO/3','CNDO/2 ','       ','       ','DFTB   ',
     5   'DFTB   '/
C *** STANDARD ATOMIC WEIGHTS TAKEN FROM THE HANDBOOK
C     OF CHEMISTRY AND PHYSICS, 58TH EDITION, 1978.
      DATA  CMS  /  1.00790D0,  4.00260D0,  6.94000D0,  9.01218D0,
     1 10.81000D0, 12.01100D0, 14.00670D0, 15.99940D0, 18.99840D0,
     2 20.17900D0, 22.98977D0, 24.30500D0, 26.98154D0, 28.08550D0,
     3 30.97376D0, 32.06000D0, 35.45300D0, 39.94800D0, 39.09830D0,
     4 40.08000D0, 44.95590D0, 47.90000D0, 50.94150D0, 51.99600D0,
     5 54.93800D0, 55.84700D0, 58.93320D0, 58.71000D0, 63.54600D0,
     6 65.38000D0, 69.73500D0, 72.59000D0, 74.92160D0, 78.96000D0,
     7 79.90400D0, 83.80000D0, 85.46780D0, 87.62000D0, 88.90590D0,
     8 91.22000D0, 92.90640D0, 95.94000D0, 98.90620D0, 101.0700D0,
     9 102.9055D0, 106.4000D0, 107.8680D0, 112.4100D0, 114.8200D0,
     A 118.6900D0, 121.7500D0, 127.6000D0, 126.9045D0, 131.3000D0,
     B 132.9054D0, 137.3300D0, 138.9100D0, 140.1200D0, 140.9070D0,
     C 144.2400D0, 145.0000D0, 150.3500D0, 151.9600D0, 157.2500D0,
     D 158.9240D0, 162.5000D0, 164.9300D0, 167.2600D0, 168.9340D0,
     E 173.0400D0, 174.9700D0, 178.4900D0, 180.9479D0, 183.8500D0,
     F 186.2070D0, 190.2000D0, 192.2200D0, 195.0900D0, 196.9665D0,
     G 200.5900D0, 204.3700D0, 207.2000D0, 208.9804D0, 208.9829D0,
     H 210.0000D0,  15.0347D0/
C *** MASSES OF THE PRINCIPAL ISOTOPES TAKEN FROM THE HANDBOOK
C     OF CHEMISTRY AND PHYSICS, 58TH EDITION, 1978.
      DATA  BMS  / 1.007825D0,  4.00260D0,  7.01600D0,  9.01218D0,
     1 11.00931D0, 12.00000D0, 14.00307D0, 15.99491D0, 18.99840D0,
     2 19.99244D0, 22.98980D0, 23.98504D0, 26.98153D0, 27.97693D0,
     3 30.99376D0, 31.97207D0, 34.96885D0, 39.96272D0, 38.96371D0,
     4 39.96259D0, 44.95592D0, 47.90000D0, 50.94400D0, 51.94050D0,
     5 54.93810D0, 55.93490D0, 58.93320D0, 57.93530D0, 62.92980D0,
     6 63.92910D0, 68.92570D0, 73.92190D0, 74.92160D0, 79.91650D0,
     7 78.91830D0, 83.80000D0, 84.91170D0, 87.90560D0, 88.90540D0,
     8 89.90430D0, 92.90600D0, 97.90550D0, 98.90620D0, 101.9037D0,
     9 102.9048D0, 105.9032D0, 106.9041D0, 113.9036D0, 114.9041D0,
     A 117.9018D0, 120.9038D0, 129.9067D0, 126.9004D0, 131.9042D0,
     B 133.9051D0, 137.9050D0, 138.9061D0, 139.9053D0, 140.9070D0,
     C 141.9075D0, 145.0000D0, 151.9195D0, 152.9209D0, 157.9241D0,
     D 159.9250D0, 163.9265D0, 164.9303D0, 165.9304D0, 168.9344D0,
     E 173.9390D0, 174.9409D0, 179.9468D0, 180.9480D0, 183.9510D0,
     F 186.9560D0, 189.9586D0, 192.9633D0, 194.9648D0, 196.9666D0,
     G 201.9706D0, 204.9745D0, 207.9766D0, 208.9804D0, 208.9829D0,
     H 210.0000D0,  15.023475D0/
C *** MULTIPLICATION TABLE FOR ABELIAN POINT GROUPS.
C     GROUP       ORDER OF IRREDUCIBLE REPRESENTATIONS.
C     Cs          A'  A''
C     C2          A   B
C     C2v         A1  A2  B1  B2
C     D2h         Ag  Au  B1g B1u B2g B2u B3g B3u
C     C2h         Ag  Bg  Au  Bu
C     D2          A   B1  B2  B3
C     Ci          Ag  Au
      DATA IPROD/ 1,  2,  3,  4,  5,  6,  7,  8,
     1            2,  1,  4,  3,  6,  5,  8,  7,
     2            3,  4,  1,  2,  7,  8,  5,  6,
     3            4,  3,  2,  1,  8,  7,  6,  5,
     4            5,  6,  7,  8,  1,  2,  3,  4,
     5            6,  5,  8,  7,  2,  1,  4,  3,
     6            7,  8,  5,  6,  3,  4,  1,  2,
     7            8,  7,  6,  5,  4,  3,  2,  1/
C *** ABELIAN POINT GROUPS HANDLED BY THE PROGRAM.
      DATA NGROUP/'Cs ','C2 ','C2v','D2h','C2h','D2 ','Ci '/
C *** LABELS FOR THE IRREDUCIBLE REPRESENTATIONS OF THESE GROUPS.
      DATA IRREP/'A''  ','A'''' ','A   ','B   ','A1  ','A2  ','B1  ',
     1       'B2  ','Ag  ','Au  ','B1g ','B1u ','B2g ','B2u ','B3g ',
     2       'B3u ','Ag  ','Bg  ','Au  ','Bu  ','A   ','B1  ','B2  ',
     3       'B3  ','Ag  ','Au  '/
C *** LOCATION OF THE FIRST LABEL FOR THE I-TH GROUP IS IRREPA(I)+1.
      DATA IRREPA/0,2,4,8,16,20,24/
C *** ORDER OF THE I-TH GROUP IS IRREPN(I).
      DATA IRREPN/2,2,4,8,4,4,2/
C *** LABELS FOR ALL COMMON POINT GROUPS.
C      1     POINT GROUP C1 (NO SYMMETRY).
C      2- 8  ABELIAN POINT GROUPS IN STANDARD ORDER.
C      9-10  LINEAR MOLECULES (NOTATION C0v,D0h).
C     11-12  LINEAR MOLECULES (NOTATION C v,D h).
C     13-52  SYMMETRIC TOPS (UP TO 8-FOLD AXIS).
C     53-59  SPHERICAL TOPS.
      DATA NPTGRP/'C1 ','Cs ','C2 ','C2v','D2h','C2h',
     1            'D2 ','Ci ','C0v','D0h','C v','D h',
     2            'C3 ','C4 ','C5 ','C6 ','C7 ','C8 ',
     3            'D3 ','D4 ','D5 ','D6 ','D7 ','D8 ',
     4            'C3v','C4v','C5v','C6v','C7v','C8v',
     5            'C3h','C4h','C5h','C6h','C7h','C8h',
     6            'D3h','D4h','D5h','D6h','D7h','D8h',
     7            'D2d','D3d','D4d','D5d','D6d','D7d',
     8            'D8d','S4 ','S6 ','S8 ','Td ','Th ',
     9            'T  ','Oh ','O  ','Ih ','I  '/
C *** TEXT FOR PROPERTIES EVALUATED (SINGULAR).
      DATA PROP1 ( 1)/'Heat of formation'/
      DATA PROP1 ( 2)/'Bond length'/
      DATA PROP1 ( 3)/'Bond angle'/
      DATA PROP1 ( 4)/'Dihedral angle'/
      DATA PROP1 ( 5)/'Ionization potential'/
      DATA PROP1 ( 6)/'Excitation energy'/
      DATA PROP1 ( 7)/'Dipole moment'/
      DATA PROP1 ( 8)/'Polarisability: alpha'/
      DATA PROP1 ( 9)/'Polarisability: beta'/
      DATA PROP1 (10)/'Polarisability: gamma'/
      DATA PROP1 (11)/'Relative energy'/
      DATA PROP1 (12)/'Dissociation energy'/
      DATA PROP1 (13)/'Activation enthalpy at 298 K'/
      DATA PROP1 (14)/'Ionization energy'/
      DATA PROP1 (15)/'Electron affinity'/
      DATA PROP1 (16)/'Vibrational wavenumber'/
      DATA PROP1 (17)/'Atomic charge'/
      DATA PROP1 (18)/'Population'/
      DATA PROP1 (19)/'S**2 expectation value'/
      DATA PROP1 (20)/'Gradient norm'/
      DATA PROP1 (21)/'Internal gradient norm'/
      DATA PROP1 (22)/'Cartesian gradient norm'/
      DATA PROP1 (23)/'IPs: Higher ionization'/
      DATA PROP1 (24)/'IPs: Difference'/
      DATA PROP1 (25)/'Dipole moment: non-ZDO'/
      DATA PROP1 (26)/'Imaginary frequencies'/
      DATA PROP1 (27)/'Wavenumber difference'/
      DATA PROP1 (28)/'Solvation energy'/
      DATA PROP1 (29)/'Solvation energy: vert'/
      DATA PROP1 (30)/'NMR shift: gas phase'/
      DATA PROP1 (31)/'NMR chemical shift'/
      DATA PROP1 (32)/'NICS chemical shift'/
      DATA PROP1 (33)/'CI excitation energy'/
      DATA PROP1 (34)/'CI rotational strength'/
      DATA PROP1 (35)/'CI oscillator strength'/
      DATA PROP1 (36)/'CI dipole or transition moment'/
      DATA PROP1 (37)/'CI dip/tr moment angle: phi'/
      DATA PROP1 (38)/'CI dip/tr moment angle: theta'/
      DATA PROP1 (39)/'CI dipole moment component'/
      DATA PROP1 (40)/'Orbital energy'/
      DATA PROP1 (41)/'Interaction energy'/
      DATA PROP1 (42)/'Atomization energy at 0 K'/
      DATA PROP1 (43)/'Proton affinity'/
      DATA PROP1 (44)/'Reaction energy'/
      DATA PROP1 (45)/'Adiabatic excitation energy'/
      DATA PROP1 (46)/'Distance RMSD'/
      DATA PROP1 (47)/'Atomization enthalpy at 298 K'/
      DATA PROP1 (48)/'Relative energy at 0 K'/
      DATA PROP1 (49)/'Barrier at 0 K'/
C *** TEXT FOR PROPERTIES EVALUATED (PLURAL).
      DATA PROP2 ( 1)/'Heats of formation'/
      DATA PROP2 ( 2)/'Bond lengths'/
      DATA PROP2 ( 3)/'Bond angles'/
      DATA PROP2 ( 4)/'Dihedral angles'/
      DATA PROP2 ( 5)/'Ionization potentials'/
      DATA PROP2 ( 6)/'Excitation energies'/
      DATA PROP2 ( 7)/'Dipole moments'/
      DATA PROP2 ( 8)/'Polarisabilities: alpha'/
      DATA PROP2 ( 9)/'Polarisabilities: beta'/
      DATA PROP2 (10)/'Polarisabilities: gamma'/
      DATA PROP2 (11)/'Enthalpy change at 298 K'/
      DATA PROP2 (12)/'Dissociation energies'/
      DATA PROP2 (13)/'Activation enthalpies at 298 K'/
      DATA PROP2 (14)/'Ionization energies'/
      DATA PROP2 (15)/'Electron affinities'/
      DATA PROP2 (16)/'Vibrational wavenumbers'/
      DATA PROP2 (17)/'Atomic charges'/
      DATA PROP2 (18)/'Populations'/
      DATA PROP2 (19)/'S**2 expectation values'/
      DATA PROP2 (20)/'Gradient norms'/
      DATA PROP2 (21)/'Internal gradient norms'/
      DATA PROP2 (22)/'Cartesian gradient norms'/
      DATA PROP2 (23)/'IPs: Higher ionizations'/
      DATA PROP2 (24)/'IPs: Differences'/
      DATA PROP2 (25)/'Dipole moments: non-ZDO'/
      DATA PROP2 (26)/'Imaginary frequencies'/
      DATA PROP2 (27)/'Wavenumber differences'/
      DATA PROP2 (28)/'Solvation energies'/
      DATA PROP2 (29)/'Solvation energies: vert'/
      DATA PROP2 (30)/'NMR shifts: gas phase'/
      DATA PROP2 (31)/'NMR chemical shifts'/
      DATA PROP2 (32)/'NICS chemical shifts'/
      DATA PROP2 (33)/'CI excitation energies'/
      DATA PROP2 (34)/'CI rotational strengths'/
      DATA PROP2 (35)/'CI oscillator strengths'/
      DATA PROP2 (36)/'CI dipole or transition moments'/
      DATA PROP2 (37)/'CI dip/tr moment angles: phi'/
      DATA PROP2 (38)/'CI dip/tr moment angles: theta'/
      DATA PROP2 (39)/'CI dipole moment components'/
      DATA PROP2 (40)/'Orbital energies'/
      DATA PROP2 (41)/'Interaction energies'/
      DATA PROP2 (42)/'Atomization energies at 0 K'/
      DATA PROP2 (43)/'Proton affinities'/
      DATA PROP2 (44)/'Reaction energies'/
      DATA PROP2 (45)/'Adiabatic excitation energies'/
      DATA PROP2 (46)/'Distance RMSDs'/
      DATA PROP2 (47)/'Atomization enthalpies at 298 K'/
      DATA PROP2 (48)/'Relative energies at 0 K'/
      DATA PROP2 (49)/'Barriers at 0 K'/
C *** TEXT FOR PROPERTIES EVALUATED (WITH UNITS).
      DATA PROP3 ( 1)/'Heats of formation (kcal/mol)'/
      DATA PROP3 ( 2)/'Bond lengths (A)'/
      DATA PROP3 ( 3)/'Bond angles (deg)'/
      DATA PROP3 ( 4)/'Dihedral angles (deg)'/
      DATA PROP3 ( 5)/'Ionization potentials (eV)'/
      DATA PROP3 ( 6)/'Excitation energies (eV)'/
      DATA PROP3 ( 7)/'Dipole moments (D)'/
      DATA PROP3 ( 8)/'Polarisabilities: alpha (A**3)'/
      DATA PROP3 ( 9)/'Polarisabilities: beta (10D-30 esu)'/
      DATA PROP3 (10)/'Polarisabilities: gamma (10D-36 esu)'/
      DATA PROP3 (11)/'Enthalpy changes at 298 K (kcal/mol)'/
      DATA PROP3 (12)/'Dissociation energies (kcal/mol)'/
      DATA PROP3 (13)/'Activation enthalpies (kcal/mol)'/
      DATA PROP3 (14)/'Ionization energies (kcal/mol)'/
      DATA PROP3 (15)/'Electron affinities (kcal/mol)'/
      DATA PROP3 (16)/'Vibrational wavenumbers (cm-1)'/
      DATA PROP3 (17)/'Atomic charges (e)'/
      DATA PROP3 (18)/'Populations (e)'/
      DATA PROP3 (19)/'S**2 expectation values'/
      DATA PROP3 (20)/'Gradient norms (kcal/mol*A)'/
      DATA PROP3 (21)/'Internal gradient norms (kcal/mol*A)'/
      DATA PROP3 (22)/'Cartesian gradient norms (kcal/mol*A)'/
      DATA PROP3 (23)/'IPs: Higher ionizations (eV)'/
      DATA PROP3 (24)/'IPs: Differences (eV)'/
      DATA PROP3 (25)/'Dipole moments: non-ZDO (D)'/
      DATA PROP3 (26)/'Imaginary frequencies'/
      DATA PROP3 (27)/'Wavenumber differences (cm-1)'/
      DATA PROP3 (28)/'Solvation energies (kcal/mol)'/
      DATA PROP3 (29)/'Solvation energies: vert (kcal/mol)'/
      DATA PROP3 (30)/'NMR shifts: gas phase (ppm)'/
      DATA PROP3 (31)/'NMR chemical shifts (ppm)'/
      DATA PROP3 (32)/'NICS chemical shifts (ppm)'/
      DATA PROP3 (33)/'CI excitation energies (eV)'/
      DATA PROP3 (34)/'CI rotational strengths (DBM)'/
      DATA PROP3 (35)/'CI oscillator strengths (au)'/
      DATA PROP3 (36)/'CI dipole or transition moment (D)'/
      DATA PROP3 (37)/'CI dip/tr moment angles: phi (deg)'/
      DATA PROP3 (38)/'CI dip/tr moment angles: theta (deg)'/
      DATA PROP3 (39)/'CI dipole moment components (D)'/
      DATA PROP3 (40)/'Orbital energies (eV)'/
      DATA PROP3 (41)/'Interaction energies (kcal/mol)'/
      DATA PROP3 (42)/'Atomization energies (kcal/mol)'/
      DATA PROP3 (43)/'Proton affinities (kcal/mol)'/
      DATA PROP3 (44)/'Reaction energies (kcal/mol)'/
      DATA PROP3 (45)/'Adiabatic excitation energies (eV)'/
      DATA PROP3 (46)/'Distance RMSD (A)'/
      DATA PROP3 (47)/'Atomization enthalpies (kcal/mol)'/
      DATA PROP3 (48)/'Relative energies at 0 K (kcal/mol)'/
      DATA PROP3 (49)/'Barriers at 0 K (kcal/mol)'/
C *** TEXT FOR UNITS OF PROPERTIES.
      DATA UNITS ( 1)/'kcal/mol'/
      DATA UNITS ( 2)/'Angstrom'/
      DATA UNITS ( 3)/'deg'/
      DATA UNITS ( 4)/'deg'/
      DATA UNITS ( 5)/'eV'/
      DATA UNITS ( 6)/'eV'/
      DATA UNITS ( 7)/'Debye'/
      DATA UNITS ( 8)/'A**3 '/
      DATA UNITS ( 9)/'10D-30 esu'/
      DATA UNITS (10)/'10D-36 esu'/
      DATA UNITS (11)/'kcal/mol'/
      DATA UNITS (12)/'kcal/mol'/
      DATA UNITS (13)/'kcal/mol'/
      DATA UNITS (14)/'kcal/mol'/
      DATA UNITS (15)/'kcal/mol'/
      DATA UNITS (16)/'cm-1'/
      DATA UNITS (17)/'e'/
      DATA UNITS (18)/'e'/
      DATA UNITS (19)/' '/
      DATA UNITS (20)/'kcal/mol*A'/
      DATA UNITS (21)/'kcal/mol*A'/
      DATA UNITS (22)/'kcal/mol*A'/
      DATA UNITS (23)/'eV'/
      DATA UNITS (24)/'eV'/
      DATA UNITS (25)/'Debye'/
      DATA UNITS (26)/' '/
      DATA UNITS (27)/'cm-1'/
      DATA UNITS (28)/'kcal/mol'/
      DATA UNITS (29)/'kcal/mol'/
      DATA UNITS (30)/'ppm'/
      DATA UNITS (31)/'ppm'/
      DATA UNITS (32)/'ppm'/
      DATA UNITS (33)/'eV'/
      DATA UNITS (34)/'DBM'/
      DATA UNITS (35)/'au'/
      DATA UNITS (36)/'D'/
      DATA UNITS (37)/'deg'/
      DATA UNITS (38)/'deg'/
      DATA UNITS (39)/'D'/
      DATA UNITS (40)/'eV'/
      DATA UNITS (41)/'kcal/mol'/
      DATA UNITS (42)/'kcal/mol'/
      DATA UNITS (43)/'kcal/mol'/
      DATA UNITS (44)/'kcal/mol'/
      DATA UNITS (45)/'eV'/
      DATA UNITS (46)/'Angstrom'/
      DATA UNITS (47)/'kcal/mol'/
      DATA UNITS (48)/'kcal/mol'/
      DATA UNITS (49)/'kcal/mol'/
      END
