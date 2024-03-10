noinst_LIBRARIES = libes.a

libes_a_SOURCES = esmain.f90 es.f90 cisenergy.f90 \
esglobal.f90 cisgrad.f90 cisnac.f90 cishess.f90 2e031.f90    \
vecutil.f90 dipole2.f90 2edrv.f 2e131.f 2e132.f cismain.f90   \
naccsf.f90 rpamain.f90 rpaenergy.f90 cissym.f90 cissolve.f90    \
sasfcisenergy.f90 sasfcisnac.f90 sasfcisz.f90 ciszsolve.f90   \
sasfcisutil.f90 cisz.f90 matutil.f90 esgradient.f90

CONST.o: ../CONST.f90
	$(F90) $(FREE)$(FFLAGS2_ES) -c ../CONST.f90

LIMIT.o: ../LIMIT.f90
	$(F90) $(FREE)$(FFLAGS2_ES) -c ../LIMIT.f90

esmain.o: esmain.f90 es.o
	$(F90) $(FREE)$(FFLAGS2_ES) -c esmain.f90

es.o: es.f90 esglobal.o LIMIT.o
	$(F90) $(FREE)$(FFLAGS2_ES) -c es.f90

cismain.o: cismain.f90 cisenergy.o cisgrad.o cishess.o cisnac.o \
	esglobal.o sasfcisenergy.o sasfcisnac.o LIMIT.o
	$(F90) $(FREE)$(FFLAGS2_ES) -c cismain.f90

cisenergy.o: cisenergy.f90 cissolve.o cissym.o esglobal.o LIMIT.o
	$(F90) $(FREE)$(FFLAGS2_ES) -c cisenergy.f90

esglobal.o: esglobal.f90 LIMIT.o
	$(F90) $(FREE)$(FFLAGS2_ES) -c esglobal.f90

vecutil.o: vecutil.f90
	$(F90) $(FREE)$(FFLAGS2_ES) -c vecutil.f90

dipole2.o: dipole2.f90 CONST.o LIMIT.o
	$(F90) $(FREE)$(FFLAGS2_ES) -c dipole2.f90

cisgrad.o: cisgrad.f90 cissolve.o ciszsolve.o esglobal.o LIMIT.o
	$(F90) $(FREE)$(FFLAGS2_ES) -c cisgrad.f90

cisnac.o: cisnac.f90 cisgrad.o cissolve.o cisz.o \
	esglobal.o esgradient.o naccsf.o LIMIT.o
	$(F90) $(FREE)$(FFLAGS2_ES) -c cisnac.f90

cishess.o: cishess.f90
	$(F90) $(FREE)$(FFLAGS2_ES) -c cishess.f90

2edrv.o: 2edrv.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_ES) -c 2edrv.f

2e031.o: 2e031.f90 LIMIT.o
	$(F90) $(FREE)$(FFLAGS2_ES) -c 2e031.f90

2e131.o: 2e131.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_ES) -c 2e131.f

2e132.o: 2e132.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_ES) -c 2e132.f

naccsf.o: naccsf.f90 LIMIT.o CONST.o
	$(F90) $(FIXED)$(FFLAGS2_ES) -c naccsf.f90

rpamain.o: rpamain.f90 esglobal.o rpaenergy.o LIMIT.o
	$(F90) $(FREE)$(FFLAGS2_ES) -c rpamain.f90

rpaenergy.o: rpaenergy.f90 cissolve.o esglobal.o LIMIT.o
	$(F90) $(FREE)$(FFLAGS2_ES) -c rpaenergy.f90

cissym.o: cissym.f90 esglobal.o
	$(F90) $(FREE)$(FFLAGS2_ES) -c cissym.f90

cissolve.o: cissolve.f90 esglobal.o
	$(F90) $(FREE)$(FFLAGS2_ES) -c cissolve.f90

sasfcisenergy.o: sasfcisenergy.f90 cissolve.o esglobal.o \
	sasfcisutil.o LIMIT.o
	$(F90) $(FREE)$(FFLAGS2_ES) -c sasfcisenergy.f90

sasfcisnac.o: sasfcisnac.f90 cissolve.o esglobal.o \
	naccsf.o sasfcisutil.o sasfcisz.o LIMIT.o
	$(F90) $(FREE)$(FFLAGS2_ES) -c sasfcisnac.f90

sasfcisz.o: sasfcisz.f90 cissolve.o ciszsolve.o esglobal.o \
	sasfcisutil.o
	$(F90) $(FREE)$(FFLAGS2_ES) -c sasfcisz.f90

sasfcisutil.o : sasfcisutil.f90 esglobal.o
	$(F90) $(FREE)$(FFLAGS2_ES) -c sasfcisutil.f90

ciszsolve.o: ciszsolve.f90 cissolve.o esglobal.o
	$(F90) $(FREE)$(FFLAGS2_ES) -c ciszsolve.f90
 
cisz.o: cisz.f90 cissolve.o ciszsolve.o esglobal.o esgradient.o
	$(F90) $(FREE)$(FFLAGS2_ES) -c cisz.f90

matutil.o: matutil.f90
	$(F90) $(FREE)$(FFLAGS2_ES) -c matutil.f90

esgradient.o: esgradient.f90 esglobal.o LIMIT.o
	$(F90) $(FREE)$(FFLAGS2_ES) -c esgradient.f90
