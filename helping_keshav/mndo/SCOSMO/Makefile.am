noinst_LIBRARIES = libscosmo.a

libscosmo_a_SOURCES =          \
      BNDPOP.f                 \
      cosmo_helper.f90         \
      cosmo_io.f90             \
      cosmo_mod.f90            \
      cosmo_util.f90           \
      CPEFERMI.f               \
      CPEMNM.f                 \
      CPEMod.f                 \
      GENSW.f                  \
      MinimizeMod.f90          \
      PELLENQ.f                \
      QPOP.f                   \
      rbMod.f90                \
      SCosmoMNDOInterfaces.f90 \
      SCosmoMNDOMod.f90        \
      SIMPLEDISP.f             \
      SSCPE.f                  \
      SurfaceMod.f90           \
      VEPMod.f90               \
      Y77LIB.f                 \
      YLib.f90

LIMIT.o: ../LIMIT.f90
	$(F90) $(FREE)$(FFLAGS2_SCOSMO) -c ../LIMIT.f90

BNDPOP.o: BNDPOP.f
	$(F90) $(FIXED)$(FFLAGS2_SCOSMO) -c BNDPOP.f

cosmo_helper.o: cosmo_helper.f90 YLib.o
	$(F90) $(FREE)$(FFLAGS2_SCOSMO) -c cosmo_helper.f90

cosmo_io.o: cosmo_io.f90 cosmo_mod.o MinimizeMod.o SurfaceMod.o YLib.o
	$(F90) $(FREE)$(FFLAGS2_SCOSMO) -c cosmo_io.f90

cosmo_mod.o: cosmo_mod.f90 cosmo_util.o rbMod.o SurfaceMod.o YLib.o
	$(F90) $(FREE)$(FFLAGS2_SCOSMO) -c cosmo_mod.f90

cosmo_util.o: cosmo_util.f90 YLib.o
	$(F90) $(FREE)$(FFLAGS2_SCOSMO) -c cosmo_util.f90

CPEFERMI.o: CPEFERMI.f
	$(F90) $(FIXED)$(FFLAGS2_SCOSMO) -c CPEFERMI.f

CPEMNM.o: CPEMNM.f
	$(F90) $(FIXED)$(FFLAGS2_SCOSMO) -c CPEMNM.f

CPEMod.o: CPEMod.f
	$(F90) $(FIXED)$(FFLAGS2_SCOSMO) -c CPEMod.f

GENSW.o: GENSW.f
	$(F90) $(FIXED)$(FFLAGS2_SCOSMO) -c GENSW.f

MinimizeMod.o: MinimizeMod.f90 cosmo_helper.o cosmo_mod.o cosmo_util.o \
               rbMod.o SurfaceMod.o YLib.o
	$(F90) $(FREE)$(FFLAGS2_SCOSMO) -c MinimizeMod.f90

PELLENQ.o: PELLENQ.f
	$(F90) $(FIXED)$(FFLAGS2_SCOSMO) -c PELLENQ.f

QPOP.o: QPOP.f
	$(F90) $(FIXED)$(FFLAGS2_SCOSMO) -c QPOP.f

rbMod.o: rbMod.f90 YLib.o
	$(F90) $(FREE)$(FFLAGS2_SCOSMO) -c rbMod.f90

SCosmoMNDOInterfaces.o: SCosmoMNDOInterfaces.f90 cosmo_mod.o cosmo_util.o  \
                        MinimizeMod.o rbMod.o SurfaceMod.o SCosmoMNDOMod.o \
                        YLib.o
	$(F90) $(FREE)$(FFLAGS2_SCOSMO) -c SCosmoMNDOInterfaces.f90

SCosmoMNDOMod.o: SCosmoMNDOMod.f90 cosmo_io.o cosmo_mod.o MinimizeMod.o \
                 SurfaceMod.o YLib.o
	$(F90) $(FREE)$(FFLAGS2_SCOSMO) -c SCosmoMNDOMod.f90

SIMPLEDISP.o: SIMPLEDISP.f
	$(F90) $(FIXED)$(FFLAGS2_SCOSMO) -c SIMPLEDISP.f

SSCPE.o: SSCPE.f
	$(F90) $(FIXED)$(FFLAGS2_SCOSMO) -c SSCPE.f

SurfaceMod.o: SurfaceMod.f90 cosmo_util.o YLib.o
	$(F90) $(FREE)$(FFLAGS2_SCOSMO) -c SurfaceMod.f90

VEPMod.o: VEPMod.f90 cosmo_helper.o cosmo_util.o SurfaceMod.o YLib.o
	$(F90) $(FREE)$(FFLAGS2_SCOSMO) -c VEPMod.f90

Y77LIB.o: Y77LIB.f
	$(F90) $(FIXED)$(FFLAGS2_SCOSMO) -c Y77LIB.f

YLib.o: YLib.f90
	$(F90) $(FREE)$(FFLAGS2_SCOSMO) -c YLib.f90
