noinst_LIBRARIES = libdftb.a

libdftb_a_SOURCES =       \
        atomen.f          \
        broyden.f         \
        dftb.f            \
        dispersion_egr.f  \
        eglcao.f          \
        ewevge.f          \
        externalchgrad.f  \
        externalshift.f   \
        fermi.f           \
        gamma.f           \
        gammamat.f        \
        gettab.f          \
        long_range.f      \
        my_gradient.f     \
        my_mulliken.f     \
        output.f          \
        repulsiv.f        \
        self.f            \
        shift.f           \
        short_range.f     \
        skpar.f           \
        slkode.f          \
        slktrafo.f        \
        dftbpar.f

LIMIT.o: ../LIMIT.f90
	$(F90) $(FREE)$(FFLAGS2_DFTB) -c ../LIMIT.f90

atomen.o: atomen.f
	$(F90) $(FIXED)$(FFLAGS2_DFTB) -c atomen.f

broyden.o: broyden.f
	$(F90) $(FIXED)$(FFLAGS2_DFTB) -c broyden.f

dftb.o: dftb.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_DFTB) -c dftb.f

dispersion_egr.o: dispersion_egr.f
	$(F90) $(FIXED)$(FFLAGS2_DFTB) -c dispersion_egr.f

eglcao.o: eglcao.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_DFTB) -c eglcao.f

ewevge.o: ewevge.f
	$(F90) $(FIXED)$(FFLAGS2_DFTB) -c ewevge.f

externalchgrad.o: externalchgrad.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_DFTB) -c externalchgrad.f

externalshift.o: externalshift.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_DFTB) -c externalshift.f

fermi.o: fermi.f
	$(F90) $(FIXED)$(FFLAGS2_DFTB) -c fermi.f

gamma.o: gamma.f
	$(F90) $(FIXED)$(FFLAGS2_DFTB) -c gamma.f

gammamat.o: gammamat.f
	$(F90) $(FIXED)$(FFLAGS2_DFTB) -c gammamat.f

gettab.o: gettab.f
	$(F90) $(FIXED)$(FFLAGS2_DFTB) -c gettab.f

long_range.o: long_range.f
	$(F90) $(FIXED)$(FFLAGS2_DFTB) -c long_range.f

my_gradient.o: my_gradient.f
	$(F90) $(FIXED)$(FFLAGS2_DFTB) -c my_gradient.f

my_mulliken.o: my_mulliken.f
	$(F90) $(FIXED)$(FFLAGS2_DFTB) -c my_mulliken.f

output.o: output.f
	$(F90) $(FIXED)$(FFLAGS2_DFTB) -c output.f

repulsiv.o: repulsiv.f
	$(F90) $(FIXED)$(FFLAGS2_DFTB) -c repulsiv.f

self.o: self.f
	$(F90) $(FIXED)$(FFLAGS2_DFTB) -c self.f

shift.o: shift.f
	$(F90) $(FIXED)$(FFLAGS2_DFTB) -c shift.f

short_range.o: short_range.f
	$(F90) $(FIXED)$(FFLAGS2_DFTB) -c short_range.f

slkode.o: slkode.f
	$(F90) $(FIXED)$(FFLAGS2_DFTB) -c slkode.f

skpar.o: skpar.f
	$(F90) $(FIXED)$(FFLAGS2_DFTB) -c skpar.f

slktrafo.o: slktrafo.f
	$(F90) $(FIXED)$(FFLAGS2_DFTB) -c slktrafo.f

dftbpar.o: dftbpar.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_DFTB) -c dftbpar.f
