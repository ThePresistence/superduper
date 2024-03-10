bin_PROGRAMS = mndo2020
MNDO = ../../mndo2020

mndo2020_SOURCES = MNDO2020.f LIMIT.f90 CONST.f90

if COND_BLAS
MAYBE_BLAS = BLAS
endif

if COND_LAPACK
MAYBE_LAPACK = LAPACK
endif

if COND_MACHINE
MAYBE_MACHINE = MACHINE
endif

if COND_PS
MAYBE_PS    = $(OPTDIR_PS)
MAYBE_TSTPS = tstps
MAYBE_ADDPS = addps
PRINT_PS    = yes
endif

if COND_GUGACI
MAYBE_GUGACI    = $(OPTDIR_GUGACI)
MAYBE_TSTGUGACI = tstguga
MAYBE_ADDGUGACI = addguga
PRINT_GUGACI    = yes
endif

if COND_PSEX
MAYBE_PSEX    = $(OPTDIR_PSEX)
MAYBE_TSTPSEX = tstpsex
MAYBE_ADDPSEX = addpsex
PRINT_PSEX    = yes
endif

if COND_HDLCOPT
MAYBE_HDLCOPT    = $(OPTDIR_HDLCOPT)
MAYBE_TSTHDLCOPT = tsthdlc
MAYBE_ADDHDLCOPT = addhdlc
PRINT_HDLCOPT    = yes
endif

if COND_CGDMS
MAYBE_CGDMS    = $(OPTDIR_CGDMS)
MAYBE_TSTCGDMS = tstlindir
MAYBE_ADDCGDMS = addlindir
PRINT_CGDMS = yes
endif

if COND_SPARSE
MAYBE_SPARSE    = $(OPTDIR_SPARSE)
MAYBE_TSTSPARSE = tstsparse
MAYBE_ADDSPARSE = addsparse
PRINT_SPARSE    = yes
endif

if COND_SCOSMO
MAYBE_SCOSMO    = $(OPTDIR_SCOSMO)
MAYBE_TSTSCOSMO = tstscosmo
MAYBE_ADDSCOSMO = addscosmo
PRINT_SCOSMO    = yes
endif

if COND_DFTB
MAYBE_DFTB    = $(OPTDIR_DFTB)
MAYBE_TSTDFTB = tstdftb
MAYBE_ADDDFTB = adddftb
PRINT_DFTB    = yes
endif

if COND_VB
MAYBE_VB    = $(OPTDIR_VB)
MAYBE_TSTVB = tstvb
MAYBE_ADDVB = addvb
PRINT_VB    = yes
endif

if COND_ES
MAYBE_ES    = $(OPTDIR_ES)
##MAYBE_TSTES = tstes
MAYBE_TSTES = 
MAYBE_ADDES = addes 
PRINT_ES    = yes
endif

if COND_NOPS
MAYBE_RMPS = rmps
endif

if COND_NOGUGACI
MAYBE_RMGUGACI = rmguga
endif

if COND_NOPSEX
MAYBE_RMPSEX = rmpsex
endif

if COND_NOHDLCOPT
MAYBE_RMHDLCOPT = rmhdlc
endif

if COND_NOCGDMS
MAYBE_RMCGDMS = rmlindir
endif

if COND_NOSPARSE
MAYBE_RMSPARSE = rmsparse
endif

if COND_NOSCOSMO
MAYBE_RMSCOSMO = rmscosmo
endif

if COND_NODFTB
MAYBE_RMDFTB = rmdftb
endif

if COND_NOVB
MAYBE_RMVB = rmvb
endif

if COND_NOES
MAYBE_RMES = rmes
endif

SUBDIRS = 2020             \
          $(MAYBE_BLAS)    \
          $(MAYBE_LAPACK)  \
          $(MAYBE_MACHINE) \
          $(MAYBE_PS)      \
          $(MAYBE_GUGACI)  \
          $(MAYBE_PSEX)    \
          $(MAYBE_HDLCOPT) \
          $(MAYBE_CGDMS)   \
          $(MAYBE_SPARSE)  \
          $(MAYBE_SCOSMO)  \
          $(MAYBE_DFTB)    \
          $(MAYBE_VB)      \
	  $(MAYBE_ES)

DISTCLEANFILES = test*/*.out  \
                 test*/*.pass \
                 test*/*.fail \
                 test*/fort.*

all: mndo2020

LIMIT.o: LIMIT.f90
	$(F90) $(FREE)$(FFLAGS2_MNDO) -c LIMIT.f90

CONST.o: CONST.f90
	$(F90) $(FREE)$(FFLAGS2_MNDO) -c CONST.f90

MNDO2020.o: MNDO2020.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_MNDO) -c MNDO2020.f

mndo2020: MNDO2020.o CONST.o $(ARCHIVES)
	$(F90) -o mndo2020 MNDO2020.o CONST.o $(LDFLAGS) $(LIBS)

libmndo.a: 2020/libmndo.a
	cp 2020/libmndo.a .

libanalytic.a: PS/libanalytic.a
	cp PS/libanalytic.a .

libguga.a: GUGACI/libguga.a
	cp GUGACI/libguga.a .

libsurface.a: PSEX/libsurface.a
	cp PSEX/libsurface.a .

libhdlc.a: HDLCOPT/libhdlc.a
	cp HDLCOPT/libhdlc.a .

liblindir.a: CGDMS/liblindir.a
	cp CGDMS/liblindir.a .

libsparse.a: SPARSE/libsparse.a
	cp SPARSE/libsparse.a .

libscosmo.a: SCOSMO/libscosmo.a
	cp SCOSMO/libscosmo.a .

libdftb.a: DFTB/libdftb.a
	cp DFTB/libdftb.a .

libvb.a: VB/libvb.a
	cp VB/libvb.a .

libbl2011.a: BLAS/libbl2011.a
	cp BLAS/libbl2011.a .

libla350.a: LAPACK/libla350.a
	cp LAPACK/libla350.a .

libmachine.a: MACHINE/libmachine.a
	cp MACHINE/libmachine.a .

libes.a: ES/libes.a
	cp ES/libes.a .

dist-hook: $(MAYBE_RMPS)       \
           $(MAYBE_RMGUGACI)   \
           $(MAYBE_RMPSEX)     \
           $(MAYBE_RMHDLCOPT)  \
           $(MAYBE_RMCGDMS)    \
           $(MAYBE_RMSPARSE)   \
           $(MAYBE_RMSCOSMO)   \
           $(MAYBE_RMDFTB)     \
           $(MAYBE_RMVB)       \
           $(MAYBE_ADDPS)      \
           $(MAYBE_ADDHDLCOPT) \
           $(MAYBE_ADDPSEX)    \
           $(MAYBE_ADDGUGACI)  \
           $(MAYBE_ADDSPARSE)  \
           $(MAYBE_ADDSCOSMO)  \
           $(MAYBE_ADDDFTB)    \
           $(MAYBE_ADDCGDMS)   \
           $(MAYBE_ADDVB)      \
	   $(MAYBE_ADDES)
	tar -c -C $(srcdir) -f - mndo2020.doc mndo2020.txt FAQ \
          testjobs/Makefile testjobs/READ_ME testjobs/check.awk \
          testjobs/*.inp testjobs/*.opt \
          testkeys/Makefile testkeys/READ_ME testkeys/check.awk \
          testkeys/*.inp testkeys/*.opt \
          testmopac6/Makefile testmopac6/READ_ME \
          testmopac6/check.awk testmopac6/testdata.dat \
          testmopac6/*.inp testmopac6/*.opt | \
	tar -x -C $(distdir) -f -

rmps:
	rm -rf $(distdir)/PS

rmguga:
	rm -rf $(distdir)/GUGACI

rmpsex:
	rm -rf $(distdir)/PSEX

rmhdlc:
	rm -rf $(distdir)/HDLCOPT

rmlindir:
	rm -rf $(distdir)/CGDMS

rmsparse:
	rm -rf $(distdir)/SPARSE

rmscosmo:
	rm -rf $(distdir)/SCOSMO

rmdftb:
	rm -rf $(distdir)/DFTB

rmvb:
	rm -rf $(distdir)/VB

rmes:
	rm -rf $(distdir)/ES

addps:
	tar -c -C $(srcdir) -f - \
          testps/Makefile testps/READ_ME testps/check.awk testps/awk.out \
          testps/*.inp testps/*.opt | \
	tar -x -C $(distdir) -f -

addguga:
	tar -c -C $(srcdir) -f - \
          testguga/Makefile testguga/check.awk testguga/overview.txt \
          testguga/*.inp testguga/*.opt | \
	tar -x -C $(distdir) -f -

addpsex:
	tar -c -C $(srcdir) -f - \
          testpsex/Makefile testpsex/check.awk \
          testpsex/*.inp testpsex/*.opt | \
	tar -x -C $(distdir) -f -

addhdlc:
	tar -c -C $(srcdir) -f - \
          testhdlc/Makefile testhdlc/READ_ME testhdlc/check.awk \
          testhdlc/*.inp testhdlc/*.opt | \
	tar -x -C $(distdir) -f -

addlindir:
	tar -c -C $(srcdir) -f - \
          testlindir/Makefile testlindir/REMARKS testlindir/check.awk \
          testlindir/*.inp testlindir/*.opt | \
	tar -x -C $(distdir) -f -

addsparse:
	tar -c -C $(srcdir) -f - \
          testsparse/Makefile testsparse/REMARKS testsparse/check.awk \
          testsparse/*.inp testsparse/*.opt | \
	tar -x -C $(distdir) -f -

addscosmo:
	tar -c -C $(srcdir) -f - \
          testscosmo/Makefile testscosmo/check.awk \
          testscosmo/*.inp testscosmo/*.opt | \
	tar -x -C $(distdir) -f -

adddftb:
	tar -c -C $(srcdir) -f - \
          testdftb/Makefile testdftb/check.awk \
          testdftb/*.inp testdftb/*.opt testdftb/SLKO/*.spl | \
	tar -x -C $(distdir) -f -

addvb:
	tar -c -C $(srcdir) -f - \
          testvb/Makefile testvb/check.awk \
          testvb/*.inp testvb/*.opt | \
	tar -x -C $(distdir) -f -

addes:
	tar -c -C $(srcdir) -f - \
          testes/Makefile testes/check.awk \
          testes/*.inp testes/*.opt | \
	tar -x -C $(distdir) -f -

test:
	@$(MAKE) -s alltests

tests:
	@$(MAKE) -s alltests

tstjobs:
	(cd testjobs;    $(MAKE) MNDO=$(MNDO))

tstkeys:
	(cd testkeys;    $(MAKE) MNDO=$(MNDO))

tstmopac6:
	(cd testmopac6;  $(MAKE) MNDO=$(MNDO))

tstps:
	(cd testps;      $(MAKE) MNDO=$(MNDO))

tstguga:
	(cd testguga;    $(MAKE) MNDO=$(MNDO))

tstpsex:
	(cd testpsex;    $(MAKE) MNDO=$(MNDO))

tsthdlc:
	(cd testhdlc;    $(MAKE) MNDO=$(MNDO))

tstlindir:
	(cd testlindir;  $(MAKE) MNDO=$(MNDO))

tstsparse:
	(cd testsparse;  $(MAKE) MNDO=$(MNDO))

tstscosmo:
	(cd testscosmo;  $(MAKE) MNDO=$(MNDO))

tstdftb:
	(cd testdftb;    $(MAKE) MNDO=$(MNDO))

tstvb:
	(cd testvb;      $(MAKE) MNDO=$(MNDO))

tstes:
	(cd testes;      $(MAKE) MNDO=$(MNDO))

alltests: tstjobs             \
          tstkeys             \
          tstmopac6           \
          $(MAYBE_TSTPS)      \
          $(MAYBE_TSTGUGACI)  \
          $(MAYBE_TSTPSEX)    \
          $(MAYBE_TSTHDLCOPT) \
          $(MAYBE_TSTCGDMS)   \
          $(MAYBE_TSTSPARSE)  \
          $(MAYBE_TSTSCOSMO)  \
          $(MAYBE_TSTDFTB)    \
          $(MAYBE_TSTVB)      \
          $(MAYBE_TSTES)      
	echo ''
	echo ''
	echo ''
	echo 'Final results of test jobs:'
	echo ''
	echo ''
	echo 'Suite                    | Tests   Passed  Failed'
	echo '-------------------------------------------------'
	(cd testjobs; \
	  ls *inp  2>/dev/null | wc | $(AWK) '{printf "Standard - formatted     |   %2d", $$1}'; \
	  ls *pass 2>/dev/null | wc | $(AWK) '{printf "      %2d", $$1}'; \
	  ls *fail 2>/dev/null | wc | $(AWK) '{printf "      %2d\n", $$1}'; \
	)
	(cd testkeys; \
	  ls *inp  2>/dev/null | wc | $(AWK) '{printf "Standard - keywords      |   %2d", $$1}'; \
	  ls *pass 2>/dev/null | wc | $(AWK) '{printf "      %2d", $$1}'; \
	  ls *fail 2>/dev/null | wc | $(AWK) '{printf "      %2d\n", $$1}'; \
	)
	(cd testmopac6; \
	  ls *inp  2>/dev/null | wc | $(AWK) '{printf "MOPAC6 Compatibility     |   %2d", $$1}'; \
	  ls *pass 2>/dev/null | wc | $(AWK) '{printf "      %2d", $$1}'; \
	  ls *fail 2>/dev/null | wc | $(AWK) '{printf "      %2d\n", $$1}'; \
	)
	test "$(PRINT_PS)" = "yes"  &&  \
	(cd testps; \
	  ls *inp  2>/dev/null | wc | $(AWK) '{printf "Analytic Derivatives     |   %2d", $$1}'; \
	  ls *pass 2>/dev/null | wc | $(AWK) '{printf "      %2d", $$1}'; \
	  ls *fail 2>/dev/null | wc | $(AWK) '{printf "      %2d\n", $$1}'; \
	)  ||  true
	test "$(PRINT_GUGACI)" = "yes"  &&  \
	(cd testguga; \
	  ls *inp  2>/dev/null | wc | $(AWK) '{printf "GUGA-CI Module           |  %3d", $$1}'; \
	  ls *pass 2>/dev/null | wc | $(AWK) '{printf "     %3d", $$1}'; \
	  ls *fail 2>/dev/null | wc | $(AWK) '{printf "     %3d\n", $$1}'; \
	)  ||  true
	test "$(PRINT_PSEX)" = "yes"  &&  \
	(cd testpsex; \
	  ls *inp  2>/dev/null | wc | $(AWK) '{printf "Surface Crossing/Hopping |   %2d", $$1}'; \
	  ls *pass 2>/dev/null | wc | $(AWK) '{printf "      %2d", $$1}'; \
	  ls *fail 2>/dev/null | wc | $(AWK) '{printf "      %2d\n", $$1}'; \
	)  ||  true
	test "$(PRINT_HDLCOPT)" = "yes"  &&  \
	(cd testhdlc; \
	  ls *inp  2>/dev/null | wc | $(AWK) '{printf "HDLC Optimizer           |   %2d", $$1}'; \
	  ls *pass 2>/dev/null | wc | $(AWK) '{printf "      %2d", $$1}'; \
	  ls *fail 2>/dev/null | wc | $(AWK) '{printf "      %2d\n", $$1}'; \
	)  ||  true
	test "$(PRINT_CGDMS)" = "yes"  &&  \
	(cd testlindir; \
	  ls *inp  2>/dev/null | wc | $(AWK) '{printf "Full-matrix DMS          |  %3d", $$1}'; \
	  ls *pass 2>/dev/null | wc | $(AWK) '{printf "     %3d", $$1}'; \
	  ls *fail 2>/dev/null | wc | $(AWK) '{printf "     %3d\n", $$1}'; \
	)  ||  true
	test "$(PRINT_SPARSE)" = "yes"  &&  \
	(cd testsparse; \
	  ls *inp  2>/dev/null | wc | $(AWK) '{printf "Direct sparse DMS        |  %3d", $$1}'; \
	  ls *pass 2>/dev/null | wc | $(AWK) '{printf "     %3d", $$1}'; \
	  ls *fail 2>/dev/null | wc | $(AWK) '{printf "     %3d\n", $$1}'; \
	)  ||  true
	test "$(PRINT_SCOSMO)" = "yes"  &&  \
	(cd testscosmo; \
	  ls *inp  2>/dev/null | wc | $(AWK) '{printf "Smooth COSMO solv. model |  %3d", $$1}'; \
	  ls *pass 2>/dev/null | wc | $(AWK) '{printf "     %3d", $$1}'; \
	  ls *fail 2>/dev/null | wc | $(AWK) '{printf "     %3d\n", $$1}'; \
	)  ||  true
	test "$(PRINT_DFTB)" = "yes"  &&  \
	(cd testdftb; \
	  ls *inp  2>/dev/null | wc | $(AWK) '{printf "Semiemp. DFTB module     |  %3d", $$1}'; \
	  ls *pass 2>/dev/null | wc | $(AWK) '{printf "     %3d", $$1}'; \
	  ls *fail 2>/dev/null | wc | $(AWK) '{printf "     %3d\n", $$1}'; \
	)  ||  true
	test "$(PRINT_VB)" = "yes"  &&  \
	(cd testvb; \
	  ls *inp  2>/dev/null | wc | $(AWK) '{printf "Valence bond module      |  %3d", $$1}'; \
	  ls *pass 2>/dev/null | wc | $(AWK) '{printf "     %3d", $$1}'; \
	  ls *fail 2>/dev/null | wc | $(AWK) '{printf "     %3d\n", $$1}'; \
	)  ||  true
	test "$(PRINT_ES)" = "yes"  &&  \
	(cd testes; \
	  ls *inp  2>/dev/null | wc | $(AWK) '{printf "Excited-state module     |  %3d", $$1}'; \
	  ls *pass 2>/dev/null | wc | $(AWK) '{printf "     %3d", $$1}'; \
	  ls *fail 2>/dev/null | wc | $(AWK) '{printf "     %3d\n", $$1}'; \
	)  ||  true
	echo '-------------------------------------------------'
	echo ''
	echo ''

testclean:
				(cd testjobs;   $(MAKE) clean)
				(cd testkeys;   $(MAKE) clean)
				(cd testmopac6; $(MAKE) clean)
	test -d testps      &&  (cd testps;     $(MAKE) clean)  ||  true
	test -d testguga    &&  (cd testguga;   $(MAKE) clean)  ||  true
	test -d testpsex    &&  (cd testpsex;   $(MAKE) clean)  ||  true
	test -d testhdlc    &&  (cd testhdlc;   $(MAKE) clean)  ||  true
	test -d testlindir  &&  (cd testlindir; $(MAKE) clean)  ||  true
	test -d testsparse  &&  (cd testsparse; $(MAKE) clean)  ||  true
	test -d testscosmo  &&  (cd testscosmo; $(MAKE) clean)  ||  true
	test -d testdftb    &&  (cd testdftb;   $(MAKE) clean)  ||  true
	test -d testvb      &&  (cd testvb;     $(MAKE) clean)  ||  true
	test -d testes      &&  (cd testes;     $(MAKE) clean)  ||  true
