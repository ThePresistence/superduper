noinst_LIBRARIES = liblindir.a

liblindir_a_SOURCES = \
	SCFCG.f       \
	ITERCG.f      \
	CGDMS.f       \
	DMSGRD.f      \
	PRECO.f       \
	DMSCOF.f      \
	CGROOT.f      \
	PURIFY.f      \
	CGUPD.f       \
	PSCAL.f       \
	BORDCG.f      \
	COMMUT.f      \
	MATDEV.f      \
	MATSCL.f      \
	MATUPM.f      \
	TRACE1.f      \
	TRACE2.f      \
	TRACE3.f      \
	TRMSUM.f      \
	DIRCG.f       \
	MATCUT.f

LIMIT.o: ../LIMIT.f90
	$(F90) $(FREE)$(FFLAGS2_CGDMS) -c ../LIMIT.f90

SCFCG.o: SCFCG.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_CGDMS) -c SCFCG.f

ITERCG.o: ITERCG.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_CGDMS) -c ITERCG.f

CGDMS.o: CGDMS.f
	$(F90) $(FIXED)$(FFLAGS2_CGDMS) -c CGDMS.f

DMSGRD.o: DMSGRD.f
	$(F90) $(FIXED)$(FFLAGS2_CGDMS) -c DMSGRD.f

PRECO.o: PRECO.f
	$(F90) $(FIXED)$(FFLAGS2_CGDMS) -c PRECO.f

DMSCOF.o: DMSCOF.f
	$(F90) $(FIXED)$(FFLAGS2_CGDMS) -c DMSCOF.f

CGROOT.o: CGROOT.f
	$(F90) $(FIXED)$(FFLAGS2_CGDMS) -c CGROOT.f

PURIFY.o: PURIFY.f
	$(F90) $(FIXED)$(FFLAGS2_CGDMS) -c PURIFY.f

CGUPD.o: CGUPD.f
	$(F90) $(FIXED)$(FFLAGS2_CGDMS) -c CGUPD.f

PSCAL.o: PSCAL.f
	$(F90) $(FIXED)$(FFLAGS2_CGDMS) -c PSCAL.f

BORDCG.o: BORDCG.f
	$(F90) $(FIXED)$(FFLAGS2_CGDMS) -c BORDCG.f

COMMUT.o: COMMUT.f
	$(F90) $(FIXED)$(FFLAGS2_CGDMS) -c COMMUT.f

MATDEV.o: MATDEV.f
	$(F90) $(FIXED)$(FFLAGS2_CGDMS) -c MATDEV.f

MATSCL.o: MATSCL.f
	$(F90) $(FIXED)$(FFLAGS2_CGDMS) -c MATSCL.f

MATUPM.o: MATUPM.f
	$(F90) $(FIXED)$(FFLAGS2_CGDMS) -c MATUPM.f

TRACE1.o: TRACE1.f
	$(F90) $(FIXED)$(FFLAGS2_CGDMS) -c TRACE1.f

TRACE2.o: TRACE2.f
	$(F90) $(FIXED)$(FFLAGS2_CGDMS) -c TRACE2.f

TRACE3.o: TRACE3.f
	$(F90) $(FIXED)$(FFLAGS2_CGDMS) -c TRACE3.f

TRMSUM.o: TRMSUM.f
	$(F90) $(FIXED)$(FFLAGS2_CGDMS) -c TRMSUM.f

DIRCG.o: DIRCG.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_CGDMS) -c DIRCG.f

MATCUT.o: MATCUT.f
	$(F90) $(FIXED)$(FFLAGS2_CGDMS) -c MATCUT.f
