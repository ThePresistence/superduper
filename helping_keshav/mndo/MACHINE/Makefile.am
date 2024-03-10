noinst_LIBRARIES = libmachine.a

libmachine_a_SOURCES =

if COND_CUDA
libmachine_a_SOURCES += chkgpu.cu
else
libmachine_a_SOURCES += nogpu.f
endif

if COND_DEFMEM
libmachine_a_SOURCES += DEFMEM.f
endif

if COND_DMR
libmachine_a_SOURCES += dmr.f
endif

if COND_GEMMBL
libmachine_a_SOURCES += gemmbl.f
endif

if COND_MEMX
libmachine_a_SOURCES += memx.c
endif

if COND_NOPGPLOT
libmachine_a_SOURCES += nopgplot.f
endif

if COND_NOPS
libmachine_a_SOURCES += nops.f
endif

if COND_NOGUGACI
libmachine_a_SOURCES += noguga.f
endif

if COND_NOPSEX
libmachine_a_SOURCES += nopsex.f
endif

if COND_NOHDLCOPT
libmachine_a_SOURCES += nohdlc.f
endif

if COND_NOCGDMS
libmachine_a_SOURCES += nolindir.f
endif

if COND_NOSPARSE
libmachine_a_SOURCES += nosparse.f
endif

if COND_NOSCOSMO
libmachine_a_SOURCES += noscosmo.f
endif

if COND_NODFTB
libmachine_a_SOURCES += nodftb.f
endif

if COND_NOVB
libmachine_a_SOURCES += novb.f
endif

if COND_NOES
libmachine_a_SOURCES += noes.f
endif

.cu.o:
	$(CXX) $(NVFLAGS_MACHINE) -c $<

DEFMEM.o: DEFMEM.f
	$(F90) $(FIXED)$(FFLAGS2_MACHINE) -c DEFMEM.f

dmr.o: dmr.f
	$(F90) $(FIXED)$(FFLAGS2_MACHINE) -c dmr.f

gemmbl.o: gemmbl.f
	$(F90) $(FIXED)$(FFLAGS2_MACHINE) -c gemmbl.f

memx.o: memx.c
	$(CC) $(CFLAGS_MACHINE) -c memx.c

nopgplot.o: nopgplot.f
	$(F90) $(FIXED)$(FFLAGS2_MACHINE) -c nopgplot.f

nops.o: nops.f
	$(F90) $(FIXED)$(FFLAGS2_MACHINE) -c nops.f

noguga.o: noguga.f
	$(F90) $(FIXED)$(FFLAGS2_MACHINE) -c noguga.f

nopsex.o: nopsex.f
	$(F90) $(FIXED)$(FFLAGS2_MACHINE) -c nopsex.f

nohdlc.o: nohdlc.f
	$(F90) $(FIXED)$(FFLAGS2_MACHINE) -c nohdlc.f

nolindir.o: nolindir.f
	$(F90) $(FIXED)$(FFLAGS2_MACHINE) -c nolindir.f

nosparse.o: nosparse.f
	$(F90) $(FIXED)$(FFLAGS2_MACHINE) -c nosparse.f

noscosmo.o: noscosmo.f
	$(F90) $(FIXED)$(FFLAGS2_MACHINE) -c noscosmo.f

nodftb.o: nodftb.f
	$(F90) $(FIXED)$(FFLAGS2_MACHINE) -c nodftb.f

novb.o: novb.f
	$(F90) $(FIXED)$(FFLAGS2_MACHINE) -c novb.f

noes.o: noes.f
	$(F90) $(FIXED)$(FFLAGS2_MACHINE) -c noes.f
